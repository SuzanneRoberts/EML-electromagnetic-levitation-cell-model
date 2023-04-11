# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 13:23:33 2015

@author: Suzanne
"""

import numpy as np, pylab as pl, unittest, time, copy

import defineDesign as dd
import discSample as ds
import plotting as pp
import plotLitData as pld
from energyBalance import sampleTemp
from samplePos import samplePos
from emlcSim import emlcSim
import scipy.optimize as op


# Set up matplotlib setting to format figures by setting mpl.rcParams, 
# see http://matplotlib.org/users/customizing.html
# and http://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
#pl.rcParams['font.size'] = 18 #22

# Dissertation settings
pl.rcParams['grid.linewidth'] = 2
pl.rcParams['lines.linewidth'] = 2
pl.rcParams['lines.markersize'] = 6
pl.rcParams['legend.fontsize'] = 12
pl.rcParams['axes.labelsize'] = 20
pl.rcParams['xtick.labelsize'] = 16
pl.rcParams['ytick.labelsize'] = 16

## Presentation settings
#pl.rcParams['grid.linewidth'] = 2
#pl.rcParams['lines.linewidth'] = 4
#pl.rcParams['lines.markersize'] = 10
#pl.rcParams['legend.fontsize'] = 22
#pl.rcParams['axes.labelsize'] = 24
#pl.rcParams['xtick.labelsize'] = 20
#pl.rcParams['ytick.labelsize'] = 20


def sumFzFunc(aSamplePos, nodes, layers, sections, slices, coil, sample, sourceSample, rPrimeSample, MAGrPRIMEsample, UNITrPRIMEsample, Verts):
    
    nodesPos, Jcomp, Jcart, Jmag, Bfield, forceField, F_lift_z, vertsPos = samplePos(aSamplePos, nodes, layers, sections, slices, coil, sample, sourceSample, rPrimeSample, MAGrPRIMEsample, UNITrPRIMEsample, Verts)
    
    # force balance
    sampleVolume  = (4/3)*np.pi*sample.R**3.0
    sampleMass    = sample.rho * sampleVolume
    sampleWeight  = sampleMass*9.81
    
    sumFz = F_lift_z - sampleWeight
    
    return sumFz, sampleWeight, nodesPos, Jcomp, Jcart, Jmag, Bfield, forceField, F_lift_z, vertsPos


# Levitation simulation function
def levSim(coil, sample, atmosphere, layers, sections, slices):
    
    tic = time.time()

    nodes, nodesSph, sourceSample, sourceSampleSph, rPrimeSample, MAGrPRIMEsample, UNITrPRIMEsample, Verts, zeniths, azimuths \
     = ds.discSample(sample.R, layers, sections, slices)
     
    # Sweep from the bottom of the coil to find the sample levitation position
    flag = 0
    aSamplePos = np.min(coil.z_i)
    sumFzOld = -1
    posVec = []
    forceVec = []
    powerVec = []
    tempVec = []
    while (flag == 0) & (aSamplePos < (np.max(coil.z_i) + 6.0*coil.R)):
        
        # wrapped function to give sumFz as the first output for optimization root-finding
        sumFz, sampleWeight, nodesPos, Jcomp, Jcart, Jmag, Bfield, forceField, F_lift_z, vertsPos = sumFzFunc(aSamplePos, nodes, layers, sections, slices, coil, sample, sourceSample, rPrimeSample, MAGrPRIMEsample, UNITrPRIMEsample, Verts)
        
        posVec.append(aSamplePos)        
        forceVec.append(F_lift_z)
        
        if (sumFzOld > 0) & (sumFz < 0):
            flag = 1
            
            # Determine exact root
            func = lambda x: (sumFzFunc(x, nodes, layers, sections, slices, coil, sample, sourceSample, rPrimeSample, MAGrPRIMEsample, UNITrPRIMEsample, Verts)[0])**2.0
            x0   = posVec[-2]
            bnds = [(posVec[-2], posVec[-1])]
            
            fmin = op.minimize(fun=func, x0=x0, method='SLSQP', bounds=bnds, options={'ftol': 1e-12})
            
            # save the levitation position
            aSamplePos = copy.copy(fmin.x) + sample.posVar
            aSamplePosGrad = copy.copy(fmin.jac[0])

#            # Jacobian check            
#            print('aSamplePosGrad',aSamplePosGrad)
#            inc = 1E-5
#            aSamplePosGradFD = (func(aSamplePos + inc) - func(aSamplePos - inc))/(2*inc)
#            print('aSamplePosGradFD',aSamplePosGradFD)
            
            nodesPos, Jcomp, Jcart, Jmag, Bfield, forceField, F_lift_z, vertsPos = samplePos(aSamplePos, nodes, layers, sections, slices, coil, sample, sourceSample, rPrimeSample, MAGrPRIMEsample, UNITrPRIMEsample, Verts)
           
        else:
            sumFzOld = copy.copy(sumFz)
            aSamplePos += (coil.R) # each coil loop contribute a "peak" to the force curve, so I do not expect force sign changes in increments smaller than the coil tube diameter
            
    if flag == 0:
        # Minimize difference when sample is too heavy to be levitated
        func = lambda x: (sumFzFunc(x, nodes, layers, sections, slices, coil, sample, sourceSample, rPrimeSample, MAGrPRIMEsample, UNITrPRIMEsample, Verts)[0])**2.0
        x0   = posVec[np.argmax(forceVec)] # start searching the position of the maximum leviation force encountered in the sweep. np.argmax returns the index of the maximum value in the argument array.
        bnds = [(np.min(coil.z_i), (np.max(coil.z_i) + 6.0*coil.R))]
        
        fmin = op.minimize(fun=func, x0=x0, method='SLSQP', bounds=bnds, options={'ftol': 1e-12})
        print('Sample not levitated')
        
        # save the levitation position
        aSamplePos = copy.copy(fmin.x) + sample.posVar
        aSamplePosGrad = copy.copy(fmin.jac[0])
        nodesPos, Jcomp, Jcart, Jmag, Bfield, forceField, F_lift_z, vertsPos = samplePos(aSamplePos, nodes, layers, sections, slices, coil, sample, sourceSample, rPrimeSample, MAGrPRIMEsample, UNITrPRIMEsample, Verts)
    
    # Plot sample levitation position inside coil
    plotLevPosWithCoil = 1
    if plotLevPosWithCoil == 1:    
        
        pl.figure()
        pl.plot(nodesPos[:,0]*1000, nodesPos[:,2]*1000, 'k.', markersize=2)
        pl.plot(coil.loops[:,0]*1000, coil.loops[:,2]*1000, 'ko', markerfacecolor='w', markersize=14, markeredgewidth=2)
        pl.plot(coil.x_i[np.where(coil.k_i == 1)]*1000, coil.z_i[np.where(coil.k_i == 1)]*1000, 'kx', markersize=10, markeredgewidth=2)
        pl.plot(coil.x_i[np.where(coil.k_i == 0)]*1000, coil.z_i[np.where(coil.k_i == 0)]*1000, 'k.', markersize=10)
        pl.plot([0.0,0.0],[np.min(coil.z_i)*1000,np.max(coil.z_i)*1000],'k--') # centerline loop2dir1
        #pl.plot([0.0,0.0],[np.min(coil.z_i), 0.015],'k--') # centerline loop1dir1
        pl.xlabel('x-axis [mm]')
        pl.ylabel('z-axis [mm]')
        pl.axis('equal')
        #pl.axis([-0.01, 0.02, -0.015, 0.015]) #loop2dir1
        #pl.axis([-0.002, 0.012, -0.002, 0.016]) #loop1dir1
        pl.show()    
        
    # Power calculation (cartesian coordinates)
    specPcomp = (Jcart[:(layers*sections),0]*np.conj(Jcart[:(layers*sections),0])) + (Jcart[:(layers*sections),1]*np.conj(Jcart[:(layers*sections),1])) + (Jcart[:(layers*sections),2]*np.conj(Jcart[:(layers*sections),2]))
    specP     = np.real(specPcomp)/sample.sigma
    P         = specP*(nodesPos[:,3]*2.0*np.pi*nodesPos[:,0]) # loop volume (not element volume)
    #P         = (nodesPos[:,3]*2.0*np.pi*nodesPos[:,0]) # loop volume - check volume
    Ptotal    = np.sum(P)
    
    powerVec.append(Ptotal)


    # Temperature calculation
    T, Qrad, Qconv = sampleTemp(Ptotal, sample, atmosphere)
    tempVec.append(T)
        
    # Plot fields in levitation position
    plotLevPosFields = 1
    if plotLevPosFields == 1:    
    #    pp.plotNodes(nodesPos[:,0], nodesPos[:,2])
    #    pp.plotSampleField(nodesPos[:,0], nodesPos[:,2], Jmag*1E-6, 'Current density, J', 'A/mm$^2$')
    
    #    pp.plot_polar_contour(Jmag*1E-6, azimuths, zeniths)
    #    pp.patchSampleField(vertsPos, Jmag*1E-6, 'Current density, J', 'A/mm$^2$')
    #    pp.plotSampleField(nodesPos[:,0], nodesPos[:,2], Jmag*1E-6, 'Current density, J', 'A/mm$^2$')
        
        pp.patchSampleField(vertsPos, specP*1E-9, 'Specific power absorbed, p', 'W/mm$^3$')
    #    pp.patchSampleField(vertsPos, np.arctan2(np.imag(Jcomp),np.real(Jcomp)), 'Phase, theta (correct quadrant)', 'Radians')
        pp.patchSampleField(vertsPos, np.sqrt(np.imag(Jcomp)**2.0+np.real(Jcomp)**2.0)*1e-6, 'Current density magnitude, |J|', 'A/mm^2')
    #    pp.plot_polar_contour(np.arctan(np.imag(Jcomp)/np.real(Jcomp)), azimuths, zeniths)
    #    pp.patchSampleField(vertsPos, np.arctan(np.imag(Jcomp)/np.real(Jcomp)), 'Phase, theta (first quadrant)', 'Radians')
    #    pp.patchSampleField(vertsPos, np.real(Bfield[:,0]), 'Real part of B_x', 'Wb/m^2')
    #    pp.patchSampleField(vertsPos, np.real(Bfield[:,1]), 'Real part of B_y', 'Wb/m^2')
        pp.patchSampleField(vertsPos, np.real(Bfield[:,2]), 'Real part of the magnetic flux in the z-direction, Re(B$_{z})$', 'Wb/m^2')
    #    pp.patchSampleField(vertsPos, np.imag(Bfield[:,0]), 'Imag part of B_x', 'Wb/m^2')
    #    pp.patchSampleField(vertsPos, np.imag(Bfield[:,1]), 'Imag part of B_y', 'Wb/m^2')
    #    pp.patchSampleField(vertsPos, np.imag(Bfield[:,2]), 'Imag part of B_z', 'Wb/m^2')
    #    pp.patchSampleField(vertsPos, forceField[:,0]*1e-9, 'Force_x', 'N/mm^3')
    #    pp.patchSampleField(vertsPos, forceField[:,0]/np.abs(forceField[:,0]), 'Sign of force_x', 'N/m^3')   
    #    pp.patchSampleField(vertsPos, forceField[:,1]*1e-9, 'Force_y', 'N/mm^3')
        pp.patchSampleField(vertsPos, forceField[:,2]*1e-9, 'Force in the z-direction, F$_z$', 'N/mm^3')
        pp.patchSampleField(vertsPos, forceField[:,2], 'Force in the z-direction, F$_z$', 'N/m^3')
    #    pp.plot_polar_contour(forceField[:,2], azimuths, zeniths)
    #    pp.patchSampleField(vertsPos, (np.arctan2(forceField[:,1],forceField[:,0])), 'Force_\{phi}', 'N/m^3')
    #    pp.patchSampleField(vertsPos, (np.arctan2(np.real(Bfield[:,1]),np.real(Bfield[:,0]))), 'B_\{phi}', 'N/m^3')
        
        pl.show()

    
    toc = time.time()
    simTime = toc - tic
   
    return (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime, aSamplePosGrad)
    
    
def meshIndependent(indCoil, indSample, indAtmosphere):
    # mesh independence study of the levitation position and sample temperature in this position
    #meshDensity = [10, 15, 20, 25, 30, 35]
    #meshDensity = [8, 9, 10, 11, 12, 13]
    #meshDensity = [8, 10, 13, 15, 17, 20]
    meshDensity = [15, 17, 20, 25, 30, 35]#, 40]

    aSamplePosVec = []
    Tvec = []
    simTimeVec = []
                   
    for idx, n in enumerate(meshDensity):
        
        levSimTuple = levSim(indCoil, indSample, indAtmosphere, n, n, n)
        # levSimTuple = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime, aSamplePosGrad)

        aSamplePosVec.append(levSimTuple[2])
        Tvec.append(levSimTuple[7])
        simTimeVec.append(levSimTuple[8])

    # convert lists to np.arrays
    meshDensity   = np.array(meshDensity)
    aSamplePosVec = np.array(aSamplePosVec)
    Tvec          = np.array(Tvec)
    simTimeVec    = np.array(simTimeVec)

    # levitation position mesh independence plots
    pl.figure()
    pl.plot(meshDensity**3.0, aSamplePosVec, '-ks', markersize = 10, linewidth = 2)
    pl.xlabel('Number of cells', fontsize = 22)
    pl.ylabel('Sample levitation position [m]', fontsize = 22)
    
    pl.figure()
    pl.plot(meshDensity**3.0, np.abs(aSamplePosVec - aSamplePosVec[-1])*100.0/aSamplePosVec[-1], '-k*', markersize = 10, linewidth = 2)
    pl.xlabel('Number of cells', fontsize = 22)
    pl.ylabel('Absolute difference in sample levitation position \n relative to case with the most cells [%]', fontsize = 22)
    
    # sample temperature mesh independence plots
    pl.figure()
    pl.plot(meshDensity**3.0, Tvec, '-ks', markersize = 10, linewidth = 2)
    pl.xlabel('Number of cells', fontsize = 22)
    pl.ylabel('Sample temperature [K]', fontsize = 22)
    
    pl.figure()
    pl.plot(meshDensity**3.0, np.abs(Tvec - Tvec[-1])*100.0/Tvec[-1], '-k*', markersize = 10, linewidth = 2)
    pl.xlabel('Number of cells', fontsize = 22)
    pl.ylabel('Absolute difference in sample temperature \n relative to case with the most cells [%]', fontsize = 22)   
    
    # simulation time mesh independence plots
    pl.figure()
    pl.plot(meshDensity**3.0, simTimeVec, '-ks', markersize = 10, linewidth = 2)
    pl.xlabel('Number of cells', fontsize = 22)
    pl.ylabel('Time required to run simulation [s]', fontsize = 22)
    
    pl.figure()
    pl.plot(meshDensity**3.0, np.abs(simTimeVec - simTimeVec[-1])*100.0/simTimeVec[-1], '-k*', markersize = 10, linewidth = 2)
    pl.xlabel('Number of cells', fontsize = 22)
    pl.ylabel('Absolute difference in simulation time \n relative to case with the most cells [%]', fontsize = 22) 
    
    pl.show()    
    

class testEMLC(unittest.TestCase):
    def testJsymmetry(self): # set function above to only compute on symmetry plane!
        # TEST - J symmetry: is the J field symmetric for a sample in the middle of a symmetric coil?
        layers   = 27
        sections = 27
        slices   = 27        
        nForce   = 1
        minForce = 0.0
        maxForce = 0.0
        
        thePosVec, theForceVec, Jcomp, powerVec = emlcSim(dd.sz, dd.Fe, dd.Ar, layers, sections, slices, nForce, minForce, maxForce)
        
        maxComparison = 0     # initialization
        for i in np.arange(sections):
            testLayer    = Jcomp[(i*sections):(i+1)*sections]
            testMod      = np.mod(len(testLayer),2)
            sampleTop    = testLayer[0:int((len(testLayer)-testMod)/2)]
            sampleBottom = testLayer[int((len(testLayer)+testMod)/2):len(testLayer)]
            comparison   = sampleTop + np.flipud(sampleBottom)
            thisMax      = np.max(np.abs(comparison))
            
            if thisMax > maxComparison:
                maxComparison = thisMax
            
#            pl.figure()
#            pl.plot(comparison)
            
        self.assertTrue(maxComparison < 1E-3)
        


pl.ion # turns interactive mode in pylab on - prevents figures from vanishing in ubuntu
#%matplotlib inline        

# main
# meshIndependent(dd.sz, dd.Fe, dd.Ar)
# emlcSim(dd.sz, dd.Fe, dd.Ar, 20, 20, 20, 1, 0.0, 0.0)

# switches
#valWithFjF    = 0
#plotWithFjF   = 0
#plotWithSzF   = 0
#plotWithMoF   = 0
#plotWithKeF   = 0
#plotRootsFigs = 0
#plotSymmGeo   = 0

#sensIplot = 0
#sensPropPlot = 0
sensTemp = 0

#varyMatProps = 0
#meshIndependence = 0
#szCenter = 0
#fjCase = 0
runUnittests = 0

#ringVerify = 0
#ring4handCompare = 0

levPos = 1
#compareRoyerTemp = 0
#labGeoVarCurrent = 0
#labGeoMoveLid    = 0


if levPos == True:
#levPos.levPos()
#def levPos():
    Layers   = 15
    Sections = 15
    Slices   = 15
    
    myCoil       = dd.sz # dd.loop2dir1
    mySample     = dd.Cu
    myAtmosphere = dd.Ar

#    myCoil = dd.ro
#    mySample = dd.Al
#    myAtmosphere = dd.Ar
#    
#    myCoil.I = 130
#    dd.ro.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
    
    print('powerVec', tupOut[5])
    print('tempVec', tupOut[6])
    print('T', tupOut[7])
    print('simTime', tupOut[8])
    
    pl.figure()
    pl.plot(tupOut[0], tupOut[1], '-bo')
    pl.plot(tupOut[2], tupOut[3], 'rx')
    pl.plot(tupOut[0], tupOut[4]*np.ones(len(tupOut[0])), 'k')
    pl.show()


if sensTemp == True:
    Layers   = 10 #30
    Sections = 10 #30
    Slices   = 10 #30
    
    myCoil = dd.roOpt
    mySample = dd.Al
    myAtmosphere = dd.Ar
    
    sensParameter = 'emissivity'
#    sensParameter = 'electricalConductivity'
#    sensParameter = 'magneticPermeability'
#    sensParameter = 'sampleRadius'
#    sensParameter = 'thermalConductivity'
#    sensParameter = 'kinematicViscosity'
#    sensParameter = 'position'
    
    # Vary the parameter
    epsilonBaseValue = 0.1
    sigmaBaseValue   = 1E7
    muBaseValue      = 4E-7*np.pi
    rBaseValue       = 0.0031
    kBaseValue       = 152E-3
    etaBaseValue      = 199E-7
    percentChangeS = 10
    percentChangeL = 50
    positionChangeS = 0.004
    positionChangeL = 0.008

    if sensParameter == 'emissivity':
        mySample.epsilon = (1 - percentChangeL/100)*epsilonBaseValue
    if sensParameter == 'electricalConductivity':
        mySample.sigma = (1 - percentChangeL/100)*sigmaBaseValue
    if sensParameter == 'magneticPermeability':
        mySample.mu0 = (1 - percentChangeL/100)*muBaseValue
    if sensParameter == 'sampleRadius':
        mySample.R = (1 - percentChangeL/100)*rBaseValue
    if sensParameter == 'thermalConductivity':
        myAtmosphere.k = (1 - percentChangeL/100)*kBaseValue
    if sensParameter == 'kinematicViscosity':
        myAtmosphere.eta = (1 - percentChangeL/100)*etaBaseValue
    if sensParameter == 'position':
        mySample.posVar = - positionChangeL
    Ivec = (np.arange(20)+17)*10
    TvecA = []
    zvecA = []
    for myI in Ivec:
        myCoil.I = myI
        dd.roOpt.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
        tupOut = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
        # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
        TvecA.append(tupOut[7])
        zvecA.append(tupOut[2])
    
    if sensParameter == 'emissivity':    
        mySample.epsilon = (1 - percentChangeS/100)*epsilonBaseValue
    if sensParameter == 'electricalConductivity':    
        mySample.sigma = (1 - percentChangeS/100)*sigmaBaseValue
    if sensParameter == 'magneticPermeability':
        mySample.mu0 = (1 - percentChangeS/100)*muBaseValue
    if sensParameter == 'sampleRadius':
        mySample.R = (1 - percentChangeS/100)*rBaseValue
    if sensParameter == 'thermalConductivity':
        myAtmosphere.k = (1 - percentChangeS/100)*kBaseValue
    if sensParameter == 'kinematicViscosity':
        myAtmosphere.eta = (1 - percentChangeS/100)*etaBaseValue
    if sensParameter == 'position':
        mySample.posVar = - positionChangeS
    Ivec = (np.arange(20)+17)*10
    TvecB = []
    zvecB = []
    for myI in Ivec:
        myCoil.I = myI
        dd.roOpt.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
        tupOut = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
        # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
        TvecB.append(tupOut[7])
        zvecB.append(tupOut[2])

    if sensParameter == 'emissivity':        
        mySample.epsilon = epsilonBaseValue
    if sensParameter == 'electricalConductivity':        
        mySample.sigma = sigmaBaseValue
    if sensParameter == 'magneticPermeability':
        mySample.mu0 = muBaseValue
    if sensParameter == 'sampleRadius':
        mySample.R = rBaseValue
    if sensParameter == 'thermalConductivity':
        myAtmosphere.k = kBaseValue
    if sensParameter == 'kinematicViscosity':
        myAtmosphere.eta = etaBaseValue
    if sensParameter == 'position':
        mySample.posVar = 0.0
    Ivec = (np.arange(20)+17)*10
    TvecC = []
    zvecC = []
    for myI in Ivec:
        myCoil.I = myI
        dd.roOpt.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
        tupOut = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
        # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
        TvecC.append(tupOut[7])
        zvecC.append(tupOut[2])

    if sensParameter == 'emissivity':
        mySample.epsilon = (1 + percentChangeS/100)*epsilonBaseValue
    if sensParameter == 'electricalConductivity':
        mySample.sigma = (1 + percentChangeS/100)*sigmaBaseValue
    if sensParameter == 'magneticPermeability':
        mySample.mu0 = (1 + percentChangeS/100)*muBaseValue
    if sensParameter == 'sampleRadius':
        mySample.R = (1 + percentChangeS/100)*rBaseValue
    if sensParameter == 'thermalConductivity':
        myAtmosphere.k = (1 + percentChangeS/100)*kBaseValue
    if sensParameter == 'kinematicViscosity':
        myAtmosphere.eta = (1 + percentChangeS/100)*etaBaseValue
    if sensParameter == 'position':
        mySample.posVar = positionChangeS
    #Ivec = (np.arange(20)+17)*10
    Ivec = (np.arange(20)+17)*10
    TvecD = []
    zvecD = []
    for myI in Ivec:
        myCoil.I = myI
        dd.roOpt.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
        tupOut = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
        # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
        TvecD.append(tupOut[7])
        zvecD.append(tupOut[2])
        
    if sensParameter == 'emissivity':
        mySample.epsilon = (1 + percentChangeL/100)*epsilonBaseValue
    if sensParameter == 'electricalConductivity':
        mySample.sigma = (1 + percentChangeL/100)*sigmaBaseValue
    if sensParameter == 'magneticPermeability':
        mySample.mu0 = (1 + percentChangeL/100)*muBaseValue
    if sensParameter == 'sampleRadius':
        mySample.R = (1 + percentChangeL/100)*rBaseValue
    if sensParameter == 'thermalConductivity':
        myAtmosphere.k = (1 + percentChangeL/100)*kBaseValue
    if sensParameter == 'kinematicViscosity':
        myAtmosphere.eta = (1 + percentChangeL/100)*etaBaseValue
    if sensParameter == 'position':
        mySample.posVar = positionChangeL
    Ivec = (np.arange(20)+17)*10
    TvecE = []
    zvecE = []
    for myI in Ivec:
        myCoil.I = myI
        dd.roOpt.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
        tupOut = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
        # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
        TvecE.append(tupOut[7])
        zvecE.append(tupOut[2])
  
  
    pl.figure()
    
#    pl.plot(Ivec, np.array(TvecA), '--mo')
#    pl.plot(Ivec, np.array(TvecB), '--r^')
#    pl.plot(Ivec, np.array(TvecC), '-g*')
#    pl.plot(Ivec, np.array(TvecD), '--bv')
#    pl.plot(Ivec, np.array(TvecE), '--cs')
    
    pl.plot(Ivec, np.array(TvecA), '--ko')
    pl.plot(Ivec, np.array(TvecB), '--k^')
    pl.plot(Ivec, np.array(TvecC), '-k*')
    pl.plot(Ivec, np.array(TvecD), '--kv')
    pl.plot(Ivec, np.array(TvecE), '--ks')

    pl.xlabel('Coil current  [A]')
    pl.ylabel('Temperature [$^\circ$C]')
    
    if sensParameter == 'emissivity':
        pl.legend(['50% decrease in sample emissivity',
                   '10% decrease in sample emissivity',
                   'Emissivity = 0.1',
                   '10% increase in sample emissivity',
                   '50% increase in sample emissivity'], loc='upper right')    
    if sensParameter == 'electricalConductivity':
        pl.legend([str(percentChangeL) + '% decrease in electrical conductivity',
                   str(percentChangeS) + '% decrease in electrical conductivity',
                   'Sample conductivity = %.1E S/m' %sigmaBaseValue,
                   str(percentChangeS) + '% increase in electrical conductivity',
                   str(percentChangeL) + '% increase in electrical conductivity'], loc='upper right') 
        pl.ylim([800,1600])    
    if sensParameter == 'magneticPermeability':
        pl.legend([str(percentChangeL) + '% decrease in magnetic permeability',
                   str(percentChangeS) + '% decrease in magnetic permeability',
                   'Magnetic permeability = %.1E H/m' %muBaseValue,
                   str(percentChangeS) + '% increase in magnetic permeability',
                   str(percentChangeL) + '% increase in magnetic permeability'], loc='upper right')
        pl.ylim([1000,2200])
    if sensParameter == 'sampleRadius':    
        pl.legend([str(percentChangeL) + '% decrease in sample radius',
                   str(percentChangeS) + '% decrease in sample radius',
                   'Sample radius = ' + str(rBaseValue) + ' m',
                   str(percentChangeS) + '% increase in sample radius',
                   str(percentChangeL) + '% increase in sample radius'], loc='upper right')  
        pl.ylim([600,1600])
    if sensParameter == 'thermalConductivity':
        pl.legend([str(percentChangeL) + '% decrease in thermal conductivity',
                   str(percentChangeS) + '% decrease in thermal conductivity',
                   'Atmosphere thermal conductivity = %.1E W/mK' %sigmaBaseValue,
                   str(percentChangeS) + '% increase in thermal conductivity',
                   str(percentChangeL) + '% increase in thermal conductivity'], loc='upper right')
        pl.ylim([1000,1900])
    if sensParameter == 'kinematicViscosity':
        pl.legend([str(percentChangeL) + '% decrease in kinematic viscosity',
                   str(percentChangeS) + '% decrease in kinematic viscosity',
                   'Atmosphere kinematic viscosity = %.1E m$^2$/s' %etaBaseValue,
                   str(percentChangeS) + '% increase in kinematic viscosity',
                   str(percentChangeL) + '% increase in kinematic viscosity'], loc='upper right')        
    if sensParameter == 'position':
        pl.legend([str(positionChangeL*1000) + 'mm decrease in levitation position',
                   str(positionChangeS*1000) + 'mm decrease in levitation position',
                   'Actual levitation position',
                   str(positionChangeS*1000) + 'mm increase in levitation position',
                   str(positionChangeL*1000) + 'mm increase in levitation position'], loc='upper left')
        pl.ylim([0,5000])
    
    pl.grid(True)               
    pl.tick_params(labelsize=20)
    
    pl.show()
    

        
# Run unittests when module is called as a script
if runUnittests == True: #__name__ == '__main__':
    
    try:
        unittest.main()
    except SystemExit as inst:
        if inst.args[0] is True: # raised by sys.exit(True) when tests failed
            raise
            
