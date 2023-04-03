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

#import sys
#sys.path.append('..')
#import valWithFjF
#import plotWithFjF
#import plotWithSzF
#import plotWithMoF
#import plotWithKeF
#import plotRootsFigs


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
#valWithFjF    = 1
#plotWithFjF   = 0
#plotWithSzF   = 0
#plotWithMoF   = 0
#plotWithKeF   = 0
#plotRootsFigs = 0
plotSymmGeo   = 1

sensIplot = 0
sensPropPlot = 0
sensTemp = 0

varyMatProps = 0
meshIndependence = 0
szCenter = 0
fjCase = 0
runUnittests = 0

ringVerify = 0
ring4handCompare = 0

levPos = 1
compareRoyerTemp = 0
labGeoVarCurrent = 0
labGeoMoveLid    = 0


if labGeoMoveLid == True:
    ticAll = time.time()    
    
    Layers   = 25
    Sections = 25
    Slices   = 25
    
    myCoil = copy.copy(dd.lb)
    mySample = dd.Al
    myAtmosphere = dd.Ar
    
    nForce   = 20
    minForce = np.min(myCoil.z_i) - 2*myCoil.R
    maxForce = np.max(myCoil.z_i) + 2*myCoil.R
    
    lidInc = 0.01
    
    mySample.R = 0.005
    
    tic = time.time()
    myCoil.z_i[9:] += lidInc*np.ones((3))
    myCoil.loops[:,2]      = myCoil.z_i
    myCoil.loopsSph[:,0]   = np.sqrt(myCoil.x_i**2.0 + myCoil.z_i**2.0)
    myCoil.loopsSph[:,1]   = np.arctan2(myCoil.x_i,myCoil.z_i)
    tupOut1 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
    thePosVec1, theForceVec1, Jcomp1, powerVec1 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    toc = time.time()
    
    myCoil.z_i[9:] += lidInc*np.ones((3))
    myCoil.loops[:,2]      = myCoil.z_i    
    myCoil.loopsSph[:,0]   = np.sqrt(myCoil.x_i**2.0 + myCoil.z_i**2.0)
    myCoil.loopsSph[:,1]   = np.arctan2(myCoil.x_i,myCoil.z_i)
    tupOut2 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec2, theForceVec2, Jcomp2, powerVec2 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    myCoil.z_i[9:] += lidInc*np.ones((3))
    myCoil.loops[:,2]      = myCoil.z_i    
    myCoil.loopsSph[:,0]   = np.sqrt(myCoil.x_i**2.0 + myCoil.z_i**2.0)
    myCoil.loopsSph[:,1]   = np.arctan2(myCoil.x_i,myCoil.z_i)
    tupOut3 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec3, theForceVec3, Jcomp3, powerVec3 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    myCoil.z_i[9:] += lidInc*np.ones((3))
    myCoil.loops[:,2]      = myCoil.z_i
    myCoil.loopsSph[:,0]   = np.sqrt(myCoil.x_i**2.0 + myCoil.z_i**2.0)
    myCoil.loopsSph[:,1]   = np.arctan2(myCoil.x_i,myCoil.z_i)
    tupOut4 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec4, theForceVec4, Jcomp4, powerVec4 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    myCoil.z_i[9:] += lidInc*np.ones((3))
    myCoil.loops[:,2]      = myCoil.z_i
    myCoil.loopsSph[:,0]   = np.sqrt(myCoil.x_i**2.0 + myCoil.z_i**2.0)
    myCoil.loopsSph[:,1]   = np.arctan2(myCoil.x_i,myCoil.z_i)
    tupOut5 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec5, theForceVec5, Jcomp5, powerVec5 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    

    mySample.R = 0.007
        
    myCoil.z_i[9:] -= 5*lidInc*np.ones((3))
    myCoil.loops[:,2]      = myCoil.z_i
    myCoil.loopsSph[:,0]   = np.sqrt(myCoil.x_i**2.0 + myCoil.z_i**2.0)
    myCoil.loopsSph[:,1]   = np.arctan2(myCoil.x_i,myCoil.z_i)
    tupOut6 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
    thePosVec6, theForceVec6, Jcomp6, powerVec6 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    toc = time.time()
    
    myCoil.z_i[9:] += lidInc*np.ones((3))
    myCoil.loops[:,2]      = myCoil.z_i
    myCoil.loopsSph[:,0]   = np.sqrt(myCoil.x_i**2.0 + myCoil.z_i**2.0)
    myCoil.loopsSph[:,1]   = np.arctan2(myCoil.x_i,myCoil.z_i)
    tupOut7 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec7, theForceVec7, Jcomp7, powerVec7 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
   
    myCoil.z_i[9:] += lidInc*np.ones((3))
    myCoil.loops[:,2]      = myCoil.z_i
    myCoil.loopsSph[:,0]   = np.sqrt(myCoil.x_i**2.0 + myCoil.z_i**2.0)
    myCoil.loopsSph[:,1]   = np.arctan2(myCoil.x_i,myCoil.z_i)
    tupOut8 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec8, theForceVec8, Jcomp8, powerVec8 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    myCoil.z_i[9:] += lidInc*np.ones((3))
    myCoil.loops[:,2]      = myCoil.z_i
    myCoil.loopsSph[:,0]   = np.sqrt(myCoil.x_i**2.0 + myCoil.z_i**2.0)
    myCoil.loopsSph[:,1]   = np.arctan2(myCoil.x_i,myCoil.z_i)
    tupOut9 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec9, theForceVec9, Jcomp9, powerVec9 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
   
    myCoil.z_i[9:] += lidInc*np.ones((3))
    myCoil.loops[:,2]      = myCoil.z_i
    myCoil.loopsSph[:,0]   = np.sqrt(myCoil.x_i**2.0 + myCoil.z_i**2.0)
    myCoil.loopsSph[:,1]   = np.arctan2(myCoil.x_i,myCoil.z_i)
    tupOut10 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec10, theForceVec10, Jcomp10, powerVec10 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    
#    pl.figure()
#    pl.plot(tupOut1[0], tupOut1[1], '-bo')
#    pl.plot(tupOut1[2], tupOut1[3], 'rx')
#    pl.plot(tupOut1[0], tupOut1[4]*np.ones(len(tupOut1[0])), 'k')
    
    pp.plotDissertation((np.array(thePosVec1)*1000, theForceVec1, 
#                         np.array(tupOut1[0])*1000, tupOut1[1],
#                         tupOut1[2]*1000, tupOut1[3],
                         np.array(thePosVec2)*1000, theForceVec2,
#                         np.array(tupOut2[0])*1000, tupOut2[1],
#                         tupOut2[2]*1000, tupOut2[3], 
                         np.array(thePosVec3)*1000, theForceVec3,
#                         np.array(tupOut3[0])*1000, tupOut3[1],
#                         tupOut3[2]*1000, tupOut3[3],
                         np.array(thePosVec4)*1000, theForceVec4,
#                         np.array(tupOut4[0])*1000, tupOut4[1],
#                         tupOut4[2]*1000, tupOut4[3], 
                         np.array(thePosVec5)*1000, theForceVec5,
#                         np.array(tupOut5[0])*1000, tupOut5[1],
#                         tupOut5[2]*1000, tupOut5[3],
                         np.array(thePosVec5)*1000, tupOut5[4]*np.ones(len(thePosVec5)) ), 
                        'forcePlotMoveLid5mmAl', 
                        xString='Sample position along the centerline of the coil [mm]', 
                        yString='Force [N]', 
                        legendList=['Levitation force, '+str(1*lidInc*1000)+'mm, '+'Sample @ '+str(round(tupOut1[2][0]*1000))+' mm', 
                                    'Levitation force, '+str(2*lidInc*1000)+'mm, '+'Sample @ '+str(round(tupOut2[2][0]*1000))+' mm', 
                                    'Levitation force, '+str(3*lidInc*1000)+'mm, '+'Sample @ '+str(round(tupOut3[2][0]*1000))+' mm',
                                    'Levitation force, '+str(4*lidInc*1000)+'mm, '+'Sample @ '+str(round(tupOut4[2][0]*1000))+' mm',
                                    'Levitation force, '+str(5*lidInc*1000)+'mm, '+'Sample @ '+str(round(tupOut5[2][0]*1000))+' mm',
                                    'Sample weight (radius = 5 mm)'], 
                        locString='lower right')
                        
    pp.plotDissertation((np.array(thePosVec6)*1000, theForceVec6, 
#                         np.array(tupOut6[0])*1000, tupOut6[1],
#                         tupOut6[2]*1000, tupOut6[3],
                         np.array(thePosVec7)*1000, theForceVec7,
#                         np.array(tupOut7[0])*1000, tupOut7[1],
#                         tupOut7[2]*1000, tupOut7[3], 
                         np.array(thePosVec8)*1000, theForceVec8,  
#                         np.array(tupOut8[0])*1000, tupOut8[1],
#                         tupOut8[2]*1000, tupOut8[3],
                         np.array(thePosVec9)*1000, theForceVec9,
#                         np.array(tupOut9[0])*1000, tupOut9[1],
#                         tupOut9[2]*1000, tupOut9[3], 
                         np.array(thePosVec10)*1000, theForceVec10, 
#                         np.array(tupOut10[0])*1000, tupOut10[1],
#                         tupOut10[2]*1000, tupOut10[3],
                         np.array(thePosVec10)*1000, tupOut10[4]*np.ones(len(thePosVec10)) ), 
                        'forcePlotMoveLid7mmAl', 
                        xString='Sample position along the centerline of the coil [mm]', 
                        yString='Force [N]', 
                        legendList=['Levitation force, '+str(1*lidInc*1000)+'mm, '+'Sample @ '+str(round(tupOut6[2][0]*1000))+' mm', 
                                    'Levitation force, '+str(2*lidInc*1000)+'mm, '+'Sample @ '+str(round(tupOut7[2][0]*1000))+' mm', 
                                    'Levitation force, '+str(3*lidInc*1000)+'mm, '+'Sample @ '+str(round(tupOut8[2][0]*1000))+' mm',
                                    'Levitation force, '+str(4*lidInc*1000)+'mm, '+'Sample @ '+str(round(tupOut9[2][0]*1000))+' mm',
                                    'Levitation force, '+str(5*lidInc*1000)+'mm, '+'Sample @ '+str(round(tupOut10[2][0]*1000))+' mm',
                                    'Sample weight (radius = 7 mm)',], 
                        locString='lower right') 
                        
    pp.plotDissertation((1000*lidInc*(np.arange(5)+1), np.array([tupOut1[7], tupOut2[7], tupOut3[7], tupOut4[7], tupOut5[7] ])-273.15,
                         1000*lidInc*(np.arange(5)+1), np.array([tupOut6[7], tupOut7[7], tupOut8[7], tupOut9[7], tupOut10[7] ])-273.15 ), 
                        'tempPlotMoveLidAl', 
                        xString='Lid height [mm]', 
                        yString='Temperature [$^{\circ}$C]', 
                        legendList=['5mm sample radius', 
                                    '7mm sample radius'], 
                        locString='upper left')

    pl.show()
    
    tocAll = time.time()
    
    print('Time to run one simmulation', toc-tic)
    print('Time to run the full labGeo', tocAll-ticAll)
    
    #meshIndependent(myCoil, mySample, myAtmosphere)


if labGeoVarCurrent == True:
    ticAll = time.time()    
    
    Layers   = 25
    Sections = 25
    Slices   = 25
    
    myCoil = dd.lb
    mySample = dd.Fe
    myAtmosphere = dd.Ar
    
    nForce   = 20
    minForce = np.min(myCoil.z_i) - 2*myCoil.R
    maxForce = np.max(myCoil.z_i) + 2*myCoil.R
    
    currentVec = np.array([100, 150, 200, 250, 300])
    
    mySample.R = 0.005
    
    tic = time.time()
    myCoil.I = currentVec[0]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut1 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
    thePosVec1, theForceVec1, Jcomp1, powerVec1 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    toc = time.time()
    
    myCoil.I = currentVec[1]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut2 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec2, theForceVec2, Jcomp2, powerVec2 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    myCoil.I = currentVec[2]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut3 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec3, theForceVec3, Jcomp3, powerVec3 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    myCoil.I = currentVec[3]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut4 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec4, theForceVec4, Jcomp4, powerVec4 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    myCoil.I = currentVec[4]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut5 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec5, theForceVec5, Jcomp5, powerVec5 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    

    mySample.R = 0.007
    
    myCoil.I = currentVec[0]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut6 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
    thePosVec6, theForceVec6, Jcomp6, powerVec6 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    toc = time.time()
    
    myCoil.I = currentVec[1]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut7 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec7, theForceVec7, Jcomp7, powerVec7 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
   
    myCoil.I = currentVec[2]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut8 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec8, theForceVec8, Jcomp8, powerVec8 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    myCoil.I = currentVec[3]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut9 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec9, theForceVec9, Jcomp9, powerVec9 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
   
    myCoil.I = currentVec[4]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut10 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec10, theForceVec10, Jcomp10, powerVec10 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    
#    pl.figure()
#    pl.plot(tupOut1[0], tupOut1[1], '-bo')
#    pl.plot(tupOut1[2], tupOut1[3], 'rx')
#    pl.plot(tupOut1[0], tupOut1[4]*np.ones(len(tupOut1[0])), 'k')
    
    pp.plotDissertation((np.array(thePosVec1)*1000, theForceVec1, 
#                         np.array(tupOut1[0])*1000, tupOut1[1],
#                         tupOut1[2]*1000, tupOut1[3],
                         np.array(thePosVec2)*1000, theForceVec2,
#                         np.array(tupOut2[0])*1000, tupOut2[1],
#                         tupOut2[2]*1000, tupOut2[3], 
                         np.array(thePosVec3)*1000, theForceVec3,
#                         np.array(tupOut3[0])*1000, tupOut3[1],
#                         tupOut3[2]*1000, tupOut3[3],
                         np.array(thePosVec4)*1000, theForceVec4,
#                         np.array(tupOut4[0])*1000, tupOut4[1],
#                         tupOut4[2]*1000, tupOut4[3], 
                         np.array(thePosVec5)*1000, theForceVec5,
#                         np.array(tupOut5[0])*1000, tupOut5[1],
#                         tupOut5[2]*1000, tupOut5[3],
                         np.array(thePosVec5)*1000, tupOut5[4]*np.ones(len(thePosVec5)) ), 
                        'forcePlotVarCurrent5mmFe', 
                        xString='Sample position along the centerline of the coil [mm]', 
                        yString='Force [N]', 
                        legendList=['Levitation force, '+str(currentVec[0])+'A, '+'Not levitated', 
                                    'Levitation force, '+str(currentVec[1])+'A, '+'Sample @ '+str(round(tupOut2[2][0]*1000,1))+' mm', 
                                    'Levitation force, '+str(currentVec[2])+'A, '+'Sample @ '+str(round(tupOut3[2][0]*1000,1))+' mm',
                                    'Levitation force, '+str(currentVec[3])+'A, '+'Sample @ '+str(round(tupOut4[2][0]*1000,1))+' mm',
                                    'Levitation force, '+str(currentVec[4])+'A, '+'Sample @ '+str(round(tupOut5[2][0]*1000,1))+' mm',
                                    'Sample weight (radius = 5 mm)'], 
                        locString='lower right')
                        
    pp.plotDissertation((np.array(thePosVec6)*1000, theForceVec6, 
#                         np.array(tupOut6[0])*1000, tupOut6[1],
#                         tupOut6[2]*1000, tupOut6[3],
                         np.array(thePosVec7)*1000, theForceVec7,
#                         np.array(tupOut7[0])*1000, tupOut7[1],
#                         tupOut7[2]*1000, tupOut7[3], 
                         np.array(thePosVec8)*1000, theForceVec8,  
#                         np.array(tupOut8[0])*1000, tupOut8[1],
#                         tupOut8[2]*1000, tupOut8[3],
                         np.array(thePosVec9)*1000, theForceVec9,
#                         np.array(tupOut9[0])*1000, tupOut9[1],
#                         tupOut9[2]*1000, tupOut9[3], 
                         np.array(thePosVec10)*1000, theForceVec10, 
#                         np.array(tupOut10[0])*1000, tupOut10[1],
#                         tupOut10[2]*1000, tupOut10[3],
                         np.array(thePosVec10)*1000, tupOut10[4]*np.ones(len(thePosVec10)) ), 
                        'forcePlotVarCurrent7mmFe', 
                        xString='Sample position along the centerline of the coil [mm]', 
                        yString='Force [N]', 
                        legendList=['Levitation force, '+str(currentVec[0])+'A, '+'Not levitated', 
                                    'Levitation force, '+str(currentVec[1])+'A, '+'Sample @ '+str(round(tupOut7[2][0]*1000,1))+' mm', 
                                    'Levitation force, '+str(currentVec[2])+'A, '+'Sample @ '+str(round(tupOut8[2][0]*1000,1))+' mm',
                                    'Levitation force, '+str(currentVec[3])+'A, '+'Sample @ '+str(round(tupOut9[2][0]*1000,1))+' mm',
                                    'Levitation force, '+str(currentVec[4])+'A, '+'Sample @ '+str(round(tupOut10[2][0]*1000,1))+' mm',
                                    'Sample weight (radius = 7 mm)',], 
                        locString='lower right') 
                        
    pp.plotDissertation((currentVec[1:], np.array([ tupOut2[7], tupOut3[7], tupOut4[7], tupOut5[7] ])-273.15,
                         currentVec[1:], np.array([ tupOut7[7], tupOut8[7], tupOut9[7], tupOut10[7] ])-273.15 ), 
                        'tempPlotVarCurrentFe', 
                        xString='Coil Current [A]', 
                        yString='Temperature [$^{\circ}$C]', 
                        legendList=['5mm sample radius', 
                                    '7mm sample radius'], 
                        locString='upper left')

    pl.show()
    
    tocAll = time.time()
    
    print('Time to run one simmulation', toc-tic)
    print('Time to run the full labGeo', tocAll-ticAll)
    
    #meshIndependent(myCoil, mySample, myAtmosphere)


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


if compareRoyerTemp == True:
    Layers   = 25
    Sections = 25
    Slices   = 25
    
    myCoil = dd.roOpt
    mySample = dd.Al
    myAtmosphere = dd.Ar
    
    # Vary the sample emissivity
    baseValue = 0.1
    percentChangeS = 70
    percentChangeL = 110

    mySample.epsilon = (1 + percentChangeL/100)*baseValue
    
    #Ivec = (np.arange(20)+17)*10
    Ivec = (np.arange(20)+17)*10
    Tvec = []
    zvec = []
    for myI in Ivec:
        myCoil.I = myI
        dd.roOpt.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
        tupOut = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
        # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
        Tvec.append(tupOut[7])
        zvec.append(tupOut[2])
  
    pld.plotWithRoyerOptTemp(Ivec, np.array(Tvec)-273.15) 
    
    pl.show()
      
##    myCoil = dd.ro
##    
##    TvecSeed = []
##    for myI in Ivec:
##        myCoil.I = myI
##        dd.roOpt.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
##        tupOut = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
##        # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
##        TvecSeed.append(tupOut[7])
#    
#    pl.figure()
#    pl.plot(Ivec, np.array(Tvec)-273.15, '-bo')
##    pl.plot(Ivec, np.array(TvecSeed)-273.15, '-ro')
#    pl.show()
#    
#    pl.figure()
#    pl.plot(Ivec, np.array(zvec), '-go')
##    pl.plot(Ivec, np.array(TvecSeed)-273.15, '-ro')
#    pl.show()


#valWithFjF.valWithFjF() 

#plotWithFjF.plotWithFjF()
        
#plotWithSzF.plotWithSzF()

#plotWithMoF.plotWithMoF()

#plotWithKeF.plotWithKeF() 

#plotRootsFigs.plotRootsFigs()



if plotSymmGeo == True:

    Layers   = 15
    Sections = 15
    Slices   = 15
    nForce   = 50
    minForce = -0.01
    maxForce =  0.01
    
    myCoil       = dd.sz
    mySample     = dd.Fe
    myAtmosphere = dd.Ar
    
    thePosVec, theForceVec, Jcomp, powerVec = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    # figure
    fig = pl.figure()
    ax = fig.add_subplot(111)
    
    pl.plot(thePosVec, theForceVec, 'k-')
    
    pl.xlabel('Sample position along the z-axis (centerline of the coil) [m]')
    pl.ylabel('Force [N]')
    
   # for loop1dir1
    ax.annotate('coil \nloops \n1 & 2', xy=(myCoil.z_i[0], -0.039), xytext=(myCoil.z_i[0] + 0.002, -0.029),
            arrowprops=dict(facecolor='black', shrink=0.08), fontsize = 16)   
    ax.annotate('coil \nloops \n3 & 4', xy=(myCoil.z_i[2], -0.039), xytext=(myCoil.z_i[2] + 0.002, -0.029),
            arrowprops=dict(facecolor='black', shrink=0.08), fontsize = 16) 
    ax.annotate('coil \nloops \n5 & 6', xy=(myCoil.z_i[4], -0.039), xytext=(myCoil.z_i[4] + 0.002, -0.029),
            arrowprops=dict(facecolor='black', shrink=0.08), fontsize = 16)   
    ax.annotate('coil \nloops \n7 & 8', xy=(myCoil.z_i[6], -0.039), xytext=(myCoil.z_i[6] + 0.002, -0.029),
            arrowprops=dict(facecolor='black', shrink=0.08), fontsize = 16)             
           
    pl.grid(True)
    pl.show()
    
    meshIndependent(myCoil, mySample, myAtmosphere)    



if sensPropPlot == True:
    
    Layers   = 15
    Sections = 15
    Slices   = 15
    nForce   = 25
    minForce = -0.005 # -0.005 # -0.045 # -0.0023   #
    maxForce = 0.005 # 0.005  #  0.030 # 0.017 #
    
    coilDesign = dd.sz #mo #sz
    material   = dd.Fe #Ni #Fe
    atmosphere = dd.Ar  
    
    #varProp = material.sigma - does not work to change the parameter in one place only -> Python pointer error?
    baseValue = 4E-7*np.pi #1E7 #7870 #0.0031 #4E-7*np.pi
    percentChangeS = 10
    percentChangeL = 50

    material.mu0 = (1 - percentChangeL/100)*baseValue
    thePosVecA, theForceVecA, Jcomp, powerVec = emlcSim(coilDesign, material, atmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    material.mu0 = (1 - percentChangeS/100)*baseValue
    thePosVecB, theForceVecB, Jcomp, powerVec = emlcSim(coilDesign, material, atmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    material.mu0 = baseValue
    thePosVecC, theForceVecC, Jcomp, powerVec = emlcSim(coilDesign, material, atmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)

    material.mu0 = (1 + percentChangeS/100)*baseValue
    thePosVecD, theForceVecD, Jcomp, powerVec = emlcSim(coilDesign, material, atmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)

    material.mu0 = (1 + percentChangeL/100)*baseValue
    thePosVecE, theForceVecE, Jcomp, powerVec = emlcSim(coilDesign, material, atmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    pl.figure()   
#    pl.plot(thePosVecA, (np.array(theForceVecA) - dd.Zn.weight), '-r*')
#    pl.plot(thePosVecB, (np.array(theForceVecB) - dd.Zn.weight), '-b^')
#    pl.plot(thePosVecC, (np.array(theForceVecC) - dd.Zn.weight), '-go')
#    pl.plot(thePosVecD, (np.array(theForceVecD) - dd.Zn.weight), '-mv')
#    pl.plot(thePosVecE, (np.array(theForceVecE) - dd.Zn.weight), '-cs')
    
    pl.plot(thePosVecA, np.array(theForceVecA), '--mo')
    pl.plot(thePosVecB, np.array(theForceVecB), '--r^')
    pl.plot(thePosVecC, np.array(theForceVecC), '-g*')
    pl.plot(thePosVecD, np.array(theForceVecD), '--bv')
    pl.plot(thePosVecE, np.array(theForceVecE), '--cs')
    
#    pl.plot(thePosVecA, np.array(theForceVecA), '--ko')
#    pl.plot(thePosVecB, np.array(theForceVecB), '--k^')
#    pl.plot(thePosVecC, np.array(theForceVecC), '-k*')
#    pl.plot(thePosVecD, np.array(theForceVecD), '--kv')
#    pl.plot(thePosVecE, np.array(theForceVecE), '--ks')
    
    
#    # Compare with Fromm & Jehn's experimental data
#    fj10cmDynes = np.loadtxt('FrommJehnExperimentalData/FrommJehnFig4-10mmDiam.txt')
#    
#    fj10 = np.zeros(np.shape(fj10cmDynes))
#    fj10[:,0] = fj10cmDynes[:,0]*1E-2
#    fj10[:,1] = fj10cmDynes[:,1]*1E-5
#    
#    pl.plot(fj10[:,0],fj10[:,1],'kd', markersize = 10)
    
    
#    # Compare with El-Kaddah and Szekely model results
#    sz200 = np.loadtxt('SzekelyModel/szFig4force.txt')
# 
#    pl.plot(sz200[:,0],sz200[:,1],'kd', markersize = 10)
    
    
#    # Compare with Moghimi model
#    moghimiModel = np.loadtxt('MoghimiModel/MoghimiFig8.txt')
#    moghimiZ     = moghimiModel[:,0]
#    moghimiF     = moghimiModel[:,1]  
#    
#    pl.plot(moghimiZ, moghimiF, 'kd', markersize = 10)
    
    
#    # Compare with Kermanpur model case
#    KermanpurModel = np.loadtxt('KermanpurModel/KermanpurFig8.txt')
#    KermanpurZ     = KermanpurModel[:,0]
#    KermanpurF     = KermanpurModel[:,1]  
#    
#    pl.plot(KermanpurZ, KermanpurF, 'kd', markersize = 10)
    
    
    pl.xlabel('Sample position along the centerline of the coil (symmetric coil symmetry plane at zero)  [m]')
    #pl.xlabel('Sample position along the centerline of the coil  [m]')
    pl.ylabel('Lifting force acting on the sample [N]')
#    pl.legend(['25% decrease in coil current',
#               '10% decrease in coil current',
#               'Coil current = 250 A',
#               '10% increase in coil current',
#               '25% increase in coil current'], loc='upper right')    
    
    pl.legend([str(percentChangeL) + '% decrease in magnetic permeability',
               str(percentChangeS) + '% decrease in magnetic permeability',
               'Magnetic permeability = %.1E H/m' %baseValue,
               str(percentChangeS) + '% increase in magnetic permeability',
               str(percentChangeL) + '% increase in magnetic permeability'], loc='upper right')
    
#    pl.legend([str(percentChangeL) + '% decrease in electrical conductivity',
#               str(percentChangeS) + '% decrease in electrical conductivity',
#               'Sample conductivity = %.1E S/m' %baseValue,
#               str(percentChangeS) + '% increase in electrical conductivity',
#               str(percentChangeL) + '% increase in electrical conductivity'], loc='upper right')    
               
#    pl.legend([str(percentChangeL) + '% decrease in sample radius',
#               str(percentChangeS) + '% decrease in sample radius',
#               'Sample radius = ' + str(baseValue) + ' m',
#               str(percentChangeS) + '% increase in sample radius',
#               str(percentChangeL) + '% increase in sample radius'], loc='upper right')                 

#    pl.legend(['-25% A model with Kermanpur et al. coil design',
#               '-10% A model with Kermanpur et al. coil design',
#               '250 A model with Kermanpur et al. coil design',
#               '+10% A model with Kermanpur et al. coil design',
#               '+25% A model with Kermanpur et al. coil design',
#               'Kermanpur et al. model results with 500 A'], loc='lower right')
    
    pl.grid(True)               
    pl.tick_params(labelsize=20)
        
    pl.show()


if sensIplot == True:
    
    Layers   = 15
    Sections = 15
    Slices   = 15
    nForce   = 25
    minForce = -0.005 # -0.005 # -0.045 # -0.0023   #
    maxForce =  0.005 # 0.005  #  0.030 # 0.017 #
    
    coilDesign = dd.sz #sz
    material   = dd.Fe #Fe
    atmosphere = dd.Ar
    
    coilDesign.I = (1 - 50/100)*250
    coilDesign.Ivec = (coilDesign.loops[:,4]*coilDesign.I) + (np.abs(coilDesign.loops[:,4]-1)*(-coilDesign.I))
    thePosVecA, theForceVecA, Jcomp, powerVec = emlcSim(coilDesign, material, atmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    coilDesign.I = (1 - 10/100)*250
    coilDesign.Ivec = (coilDesign.loops[:,4]*coilDesign.I) + (np.abs(coilDesign.loops[:,4]-1)*(-coilDesign.I))
    thePosVecB, theForceVecB, Jcomp, powerVec = emlcSim(coilDesign, material, atmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    coilDesign.I = 250
    coilDesign.Ivec = (coilDesign.loops[:,4]*coilDesign.I) + (np.abs(coilDesign.loops[:,4]-1)*(-coilDesign.I))
    thePosVecC, theForceVecC, Jcomp, powerVec = emlcSim(coilDesign, material, atmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)

    coilDesign.I = (1 + 10/100)*250
    coilDesign.Ivec = (coilDesign.loops[:,4]*coilDesign.I) + (np.abs(coilDesign.loops[:,4]-1)*(-coilDesign.I))
    thePosVecD, theForceVecD, Jcomp, powerVec = emlcSim(coilDesign, material, atmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)

    coilDesign.I = (1 + 50/100)*250
    coilDesign.Ivec = (coilDesign.loops[:,4]*coilDesign.I) + (np.abs(coilDesign.loops[:,4]-1)*(-coilDesign.I))
    thePosVecE, theForceVecE, Jcomp, powerVec = emlcSim(coilDesign, material, atmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    pl.figure()
   
    
#    pl.plot(thePosVecA, (np.array(theForceVecA) - dd.Zn.weight), '-r*')
#    pl.plot(thePosVecB, (np.array(theForceVecB) - dd.Zn.weight), '-b^')
#    pl.plot(thePosVecC, (np.array(theForceVecC) - dd.Zn.weight), '-go')
#    pl.plot(thePosVecD, (np.array(theForceVecD) - dd.Zn.weight), '-mv')
#    pl.plot(thePosVecE, (np.array(theForceVecE) - dd.Zn.weight), '-cs')
    
    pl.plot(thePosVecA, np.array(theForceVecA), '--mo')
    pl.plot(thePosVecB, np.array(theForceVecB), '--r^')
    pl.plot(thePosVecC, np.array(theForceVecC), '-g*')
    pl.plot(thePosVecD, np.array(theForceVecD), '--bv')
    pl.plot(thePosVecE, np.array(theForceVecE), '--cs')
    
#    pl.plot(thePosVecA, np.array(theForceVecA), '--ko')
#    pl.plot(thePosVecB, np.array(theForceVecB), '--k^')
#    pl.plot(thePosVecC, np.array(theForceVecC), '-k*')
#    pl.plot(thePosVecD, np.array(theForceVecD), '--kv')
#    pl.plot(thePosVecE, np.array(theForceVecE), '--ks')    
    
#    # Plot Szekely coil loop positions
#    pl.plot(  0.005*np.ones(2), np.array([-0.015, 0.015]), '-k')
#    pl.plot( 0.0081*np.ones(2), np.array([-0.015, 0.015]), '-k')
#    pl.plot( -0.005*np.ones(2), np.array([-0.015, 0.015]), '-k')
#    pl.plot(-0.0081*np.ones(2), np.array([-0.015, 0.015]), '-k')
    
    
#    # Compare with Fromm & Jehn's experimental data
#    fj10cmDynes = np.loadtxt('FrommJehnExperimentalData/FrommJehnFig4-10mmDiam.txt')
#    
#    fj10 = np.zeros(np.shape(fj10cmDynes))
#    fj10[:,0] = fj10cmDynes[:,0]*1E-2
#    fj10[:,1] = fj10cmDynes[:,1]*1E-5
#    
#    pl.plot(fj10[:,0],fj10[:,1],'kd', markersize = 10)
    
    
#    # Compare with El-Kaddah and Szekely model results
#    sz200 = np.loadtxt('SzekelyModel/szFig4force.txt')
# 
#    pl.plot(sz200[:,0],sz200[:,1],'kd', markersize = 10)
    
    
#    # Compare with Moghimi model
#    moghimiModel = np.loadtxt('MoghimiModel/MoghimiFig8.txt')
#    moghimiZ     = moghimiModel[:,0]
#    moghimiF     = moghimiModel[:,1]  
#    
#    pl.plot(moghimiZ, moghimiF, 'kd', markersize = 10)
    
    
#    # Compare with Kermanpur model case
#    KermanpurModel = np.loadtxt('KermanpurModel/KermanpurFig8.txt')
#    KermanpurZ     = KermanpurModel[:,0]
#    KermanpurF     = KermanpurModel[:,1]  
#    
#    pl.plot(KermanpurZ, KermanpurF, 'kd', markersize = 10)
    
    
    pl.xlabel('Sample position along the centerline of the coil (symmetric coil symmetry plane at zero)  [m]')
#    pl.xlabel('Sample position along the z-axis (centerline of the coil) [m]. \n Symmetric coil loops are at z = - 0.008, z = - 0.005, z = 0.005, z = 0.008,  with the symmetry plane at z = 0.')
#    pl.xlabel('Sample position along the centerline of the coil  [m]')
    pl.ylabel('Lifting force acting on the sample [N]')
    pl.legend(['50% decrease in coil current',
               '10% decrease in coil current',
               'Coil current = 250 A',
               '10% increase in coil current',
               '50% increase in coil current'], loc='upper right')    
#    pl.legend(['-25% A model with Kermanpur et al. coil design',
#               '-10% A model with Kermanpur et al. coil design',
#               '250 A model with Kermanpur et al. coil design',
#               '+10% A model with Kermanpur et al. coil design',
#               '+25% A model with Kermanpur et al. coil design',
#               'Kermanpur et al. model results with 500 A'], loc='lower right', fontsize=18)
    
    pl.grid(True)               
    pl.tick_params(labelsize=20)
        
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
    
   

if varyMatProps == True:
    
    Layers   = 20
    Sections = 20
    Slices   = 20
    nForce   = 1
    minForce =  0.010 # -0.045
    maxForce =  0.010 # 0.030

    scaleFactor = np.array([0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 4.0, 6.0, 8.0, 10])
    force10s = []
    
    for n in scaleFactor:
        
#        #dd.Fe.mu0 = n*4E-7*np.pi
#        dd.Cu.sigma = n*5.96E7      
#        thePosVec10, theForceVec10, Jcomp, powerVec = emlcSim(dd.fj, dd.Cu, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
#        force10s.append(theForceVec10[0])
        
#        dd.Cu.sigma = n*5.96E7
#        dd.Cu.R = 0.006
#        thePosVec12, theForceVec12, Jcomp, powerVec = emlcSim(dd.fj, dd.Cu, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
#        
#        dd.Cu.sigma = n*5.96E7
#        dd.Cu.R = 0.0075
#        thePosVec15, theForceVec15, Jcomp, powerVec = emlcSim(dd.fj, dd.Cu, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
#        
        dd.Cu.sigma = n*5.96E7
        dd.Cu.R = 0.010
        thePosVec20, theForceVec20, Jcomp, powerVec = emlcSim(dd.fj, dd.Cu, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
        force10s.append(theForceVec20[0])
        
#        pl.figure(101)
#        #pl.plot(thePosVec, theForceVec)
#       
#        #pld.plotWithMoghimiCaseForce(thePosVec, theForceVec)
#        #pld.plotWithKermanpurCaseForce(thePosVec, (theForceVec - dd.Zn.weight)) #/20
#        #pl.legend(['Model (current implementation) divided by 20','Kermanpur et al. model reported results','Kermanpur et al. reported data points'], loc='lower right', fontsize=22)
#        pld.plotWithFrommJehnCaseForce(thePosVec10, theForceVec10, thePosVec12, theForceVec12, thePosVec15, theForceVec15, thePosVec20, theForceVec20)
#        pl.legend(['Fromm & Jehn experimental data, d=20mm',
#               'Model (current implementation), d=20mm',
#               'Fromm & Jehn experimental data, d=15mm',
#               'Model (current implementation), d=15mm',
#               'Fromm & Jehn experimental data, d=12mm',
#               'Model (current implementation), d=12mm',
#               'Fromm & Jehn experimental data, d=10mm',
#               'Model (current implementation), d=10mm'], loc='upper right', fontsize=18)
    pl.figure(102)
    pl.plot(scaleFactor, np.array(force10s), '-bo')
        
    pl.show()
        


if meshIndependence == True:
    # mesh independence study of the force curve
    meshDensity = [10, 15, 20, 25, 30, 35]
    #meshDensity = [8, 9, 10, 11, 12, 13]
    
    forceIndVecFj = []
    forceIndVecSz = []
    forceIndVecMo = []
    # forceIndVecKe = []
    
    powerAbsVecFj = []
    powerAbsVecSz = []
    powerAbsVecMo = []
    # powerAbsVecKe = []
    
    # initialize old vectors
    theForceVecFj = []
    theForceVecSz = []
    theForceVecMo = []
    
    powerVecFj = []
    powerVecSz = []
    powerVecMo = []
    
    # initialize old vectors
    oldForceVecFj = []
    oldForceVecSz = []
    oldForceVecMo = []
    
    oldpowerVecFj = []
    oldpowerVecSz = []
    oldpowerVecMo = []
    
    for idx, n in enumerate(meshDensity):
                
        oldForceVecFj = theForceVecFj
        oldForceVecSz = theForceVecSz
        oldForceVecMo = theForceVecMo
        
        oldpowerVecFj = powerVecFj
        oldpowerVecSz = powerVecSz
        oldpowerVecMo = powerVecMo        
        
        Layers   = n
        Sections = n
        Slices   = n
        nForce   = 17
        minForce = -0.04
        maxForce =  0.03
        
        thePosVecFj, theForceVecFj, JcompFj, powerVecFj = emlcSim(dd.fj, dd.Cu, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
        thePosVecSz, theForceVecSz, JcompSz, powerVecSz = emlcSim(dd.sz, dd.Fe, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
        thePosVecMo, theForceVecMo, JcompMo, powerVecMo = emlcSim(dd.mo, dd.Ni, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
        # thePosVecKe, theForceVecKe, JcompKe, powerVecKe = emlcSim(dd.ke, dd.Zn, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
        
        forceIndVecFj.append(np.max(theForceVecFj))
        forceIndVecSz.append(np.max(theForceVecSz))
        forceIndVecMo.append(np.max(theForceVecMo))
        # forceIndVecKe.append(np.max(theForceVecKe))
        
        powerAbsVecFj.append(np.max(powerVecFj))
        powerAbsVecSz.append(np.max(powerVecSz))
        powerAbsVecMo.append(np.max(powerVecMo))
        # powerAbsVecKe.append(powerVecKe[8])
        
        pl.figure(100)
        pl.plot(thePosVecFj, theForceVecFj, 'b', linewidth = 1)
        pl.plot(thePosVecSz, theForceVecSz, 'r', linewidth = 1)
        pl.plot(thePosVecMo, theForceVecMo, 'g', linewidth = 1)
        # pl.plot(thePosVecKe, theForceVecKe)
        pl.hold(True)
        
        pl.figure(101)
        pl.plot(thePosVecFj, powerVecFj, 'b', linewidth = 1)
        pl.plot(thePosVecSz, powerVecSz, 'r', linewidth = 1)
        pl.plot(thePosVecMo, powerVecMo, 'g', linewidth = 1)
        # pl.plot(thePosVecKe, powerVecKe)
        pl.hold(True)

        
        if idx > 0:
            pl.figure(102)
            pl.plot(thePosVecFj, np.abs(theForceVecFj - oldForceVecFj), 'b', linewidth = 1)
            pl.plot(thePosVecSz, np.abs(theForceVecSz - oldForceVecSz), 'r', linewidth = 1)
            pl.plot(thePosVecMo, np.abs(theForceVecMo - oldForceVecMo), 'g', linewidth = 1)
            # pl.plot(thePosVecKe, theForceVecKe)
            pl.hold(True)
            
            pl.figure(103)
            pl.plot(thePosVecFj, np.abs(powerVecFj - oldpowerVecFj), 'b', linewidth = 1)
            pl.plot(thePosVecSz, np.abs(powerVecSz - oldpowerVecSz), 'r', linewidth = 1)
            pl.plot(thePosVecMo, np.abs(powerVecMo - oldpowerVecMo), 'g', linewidth = 1)
            # pl.plot(thePosVecKe, powerVecKe)
            pl.hold(True)

    # convert lists to np.arrays
    meshDensity   = np.array(meshDensity)
    forceIndVecFj = np.array(forceIndVecFj)
    forceIndVecSz = np.array(forceIndVecSz)
    forceIndVecMo = np.array(forceIndVecMo)
    # forceIndVecKe = np.array(forceIndVecKe)

    pl.figure(100)
    pl.xlabel('Sample position along the centerline of the coil [m]', fontsize=24)
    pl.ylabel('Lifting force acting on the sample [N]', fontsize=24)
    pl.legend(['Fromm & Jehn','El-Kaddah & Szekely','Moghimi et al.'], loc='upper right', fontsize = 22)
    # pl.legend([str(meshDensity[0]), str(meshDensity[1]), str(meshDensity[2]), str(meshDensity[3]), str(meshDensity[4]), str(meshDensity[5])])

    
    pl.figure(101)
    pl.xlabel('Sample position along the centerline of the coil [m]', fontsize=24)
    pl.ylabel('Power absorbed by the sample [W]', fontsize=24)
    pl.legend(['Fromm & Jehn','El-Kaddah & Szekely','Moghimi et al.'], loc='upper right', fontsize = 22)
    # pl.legend([str(meshDensity[0]), str(meshDensity[1]), str(meshDensity[2]), str(meshDensity[3]), str(meshDensity[4]), str(meshDensity[5])])
    pl.show()
    
    
    # lifting force mesh independence plots
    pl.figure()
    pl.plot(meshDensity**3.0, forceIndVecFj, '-b*', markersize = 10, linewidth = 2)
    pl.plot(meshDensity**3.0, forceIndVecSz, '-ro', markersize = 10, linewidth = 2)
    pl.plot(meshDensity**3.0, forceIndVecMo, '-g^', markersize = 10, linewidth = 2)
    # pl.plot(meshDensity**3.0, forceIndVecKe, '-ks', markersize = 10)
    pl.xlabel('Number of cells', fontsize = 22)
    pl.ylabel('Maximum lifting force [N]', fontsize = 22)
    pl.legend(['Fromm & Jehn','El-Kaddah & Szekely','Moghimi et al.'], loc='right', fontsize = 22)
    
    pl.figure()
    pl.plot(meshDensity**3.0, np.abs(forceIndVecFj - forceIndVecFj[-1])*100.0/forceIndVecFj[-1], '-b*', markersize = 10, linewidth = 2)
    pl.plot(meshDensity**3.0, np.abs(forceIndVecSz - forceIndVecSz[-1])*100.0/forceIndVecSz[-1], '-ro', markersize = 10, linewidth = 2)
    pl.plot(meshDensity**3.0, np.abs(forceIndVecMo - forceIndVecMo[-1])*100.0/forceIndVecMo[-1], '-g^', markersize = 10, linewidth = 2)
    # pl.plot(meshDensity**3.0, np.abs(forceIndVecKe - forceIndVecKe[-1])*100.0/forceIndVecKe[-1], '-ks', markersize = 10)
    pl.xlabel('Number of cells', fontsize = 22)
    pl.ylabel('Absolute difference in maximum lifting force \n relative to case with the most cells [%]', fontsize = 22)
    pl.legend(['Fromm & Jehn','El-Kaddah & Szekely','Moghimi et al.'], loc='upper right', fontsize = 22)
    
    # power absorbed mesh independence plots
    pl.figure()
    pl.plot(meshDensity**3.0, powerAbsVecFj, '-b*', markersize = 10, linewidth = 2)
    pl.plot(meshDensity**3.0, powerAbsVecSz, '-ro', markersize = 10, linewidth = 2)
    pl.plot(meshDensity**3.0, powerAbsVecMo, '-g^', markersize = 10, linewidth = 2)
    # pl.plot(meshDensity**3.0, powerAbsVecKe, '-ks', markersize = 10)
    pl.xlabel('Number of cells', fontsize = 22)
    pl.ylabel('Maximum power absorbed [W]', fontsize = 22)
    pl.legend(['Fromm & Jehn','El-Kaddah & Szekely','Moghimi et al.'], loc='right', fontsize = 22)
    
    pl.figure()
    pl.plot(meshDensity**3.0, np.abs(powerAbsVecFj - powerAbsVecFj[-1])*100.0/powerAbsVecFj[-1], '-b*', markersize = 10, linewidth = 2)
    pl.plot(meshDensity**3.0, np.abs(powerAbsVecSz - powerAbsVecSz[-1])*100.0/powerAbsVecSz[-1], '-ro', markersize = 10, linewidth = 2)
    pl.plot(meshDensity**3.0, np.abs(powerAbsVecMo - powerAbsVecMo[-1])*100.0/powerAbsVecMo[-1], '-g^', markersize = 10, linewidth = 2)
    # pl.plot(meshDensity**3.0, np.abs(powerAbsVecKe - powerAbsVecKe[-1])*100.0/powerAbsVecKe[-1], '-ks', markersize = 10)
    pl.xlabel('Number of cells', fontsize = 22)
    pl.ylabel('Absolute difference in maximum power absorbed \n relative to case with the most cells [%]', fontsize = 22)
    pl.legend(['Fromm & Jehn','El-Kaddah & Szekely','Moghimi et al.'], loc='upper right', fontsize = 22)
    
#    # induced current mesh independence plots
#    pl.figure()
#    pl.plot(meshDensity**3.0, np.sum(np.real(JcompFj)), '-b*')
#    pl.plot(meshDensity**3.0, np.sum(np.real(JcompSz)), '-ro')
#    pl.plot(meshDensity**3.0, np.sum(np.real(JcompMo)), '-g^')
#    pl.plot(meshDensity**3.0, np.sum(np.real(JcompKe)), '-ks')
#    pl.xlabel('Number of cells', fontsize = 22)
#    pl.ylabel('Total real current induced [A]', fontsize = 22)
#    pl.legend(['Fromm & Jehn','El-Kaddah & Szekely','Moghimi','Kermanpur'], loc='upper right', fontsize = 22)
   
    pl.show()    
    
   
if szCenter == True:
    
    meshDensity = [25]
    
    for n in meshDensity:
        Layers   = n
        Sections = n
        Slices   = n
        nForce   = 1
        minForce = 0.005
        maxForce = 0.005        
        
        thePosVec, theForceVec, Jcomp, powerVec = emlcSim(dd.sz, dd.Fe, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
        
        pl.figure(100)
        pl.plot(thePosVec, theForceVec)
        pl.hold(True)
        
#    pl.figure(101)
#    pl.legend([str(meshDensity[0]), str(meshDensity[1]), str(meshDensity[2])])   
   
if fjCase == True:
    layers   = 25
    sections = 25
    slices   = 25
    nForce   = 5
    minForce = 0.0
    maxForce = 0.01
    
    thePosVec, theForceVec, Jcomp, powerVec = emlcSim(dd.fj, dd.Cu, dd.Ar, layers, sections, slices, nForce, minForce, maxForce)
    
    pl.figure()
    pl.plot(thePosVec, theForceVec)
    
    
if ring4handCompare == True:
    layers   = 1
    sections = 1
    slices   = 2
    nForce   = 1
    minForce = 0.0
    maxForce = 0.0
    
    thePosVec, theForceVec, Jcomp, powerVec = emlcSim(dd.ring, dd.CuRing, dd.Ar, layers, sections, slices, nForce, minForce, maxForce)
    
    
if ringVerify == True:
    layers   = 1
    sections = 1
    slices   = 400
    nForce   = 1
    minForce = 0.0
    maxForce = 0.0
    
  
    #ringRadii = np.array([0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1])
    #ringRadii = np.array([0.005, 0.01, 0.02, 0.04, 0.05, 0.06, 0.08, 0.1, 0.2, 0.5, 1.0, 1.5, 5.0, 10.0])
    #ringRadii = np.array([0.5, 1.0, 1.5, 5.0, 10.0])
    ringRadii = (np.arange(40) + 1)/2.0
    Fvec      = []
    FanaVec   = []
    indIVec   = []
    indImagVec = []
    
    for radius in ringRadii:
        # change loop radii
        dd.ring.x_i = np.array([radius])
        dd.ring.loops[:,0] = np.array([radius])
        dd.ring.loopsSph[:,0]   = np.sqrt(dd.ring.loops[:,0]**2.0 + dd.ring.loops[:,2]**2.0)
        dd.ring.loopsSph[:,1]   = np.arctan2(dd.ring.loops[:,0], dd.ring.loops[:,2])
        dd.CuRing.R = radius
        
       
        # run model
        thePosVec, theForceVec, Jcomp, powerVec = emlcSim(dd.ring, dd.CuRing, dd.Ar, layers, sections, slices, nForce, minForce, maxForce)
        Fvec.append(theForceVec[0]*np.pi*0.226**2.0)
        
        # analytical calulation
        print(4E-7*np.pi)
        print(dd.ring.I)
        print(np.real(Jcomp)*np.pi*0.226**2.0)
        print(np.sqrt(np.real(Jcomp)**2.0 + np.imag(Jcomp)**2.0)*np.pi*0.226**2.0)
        print(2.0*np.pi*radius)
        print(dd.ring.z_i - maxForce)
        #Fana = (4E-7*np.pi)*dd.ring.I*(np.real(Jcomp)*np.pi*0.005**2.0)*(2.0*np.pi*radius)/(2.0*np.pi*(dd.ring.z_i - maxForce))
        #FanaVec.append(Fana)
        #print(Fana)
        
        indIVec.append(np.real(Jcomp)*np.pi*0.226**2.0)
        indImagVec.append(np.sqrt(np.real(Jcomp)**2.0 + np.imag(Jcomp)**2.0)*np.pi*0.226**2.0)
        
        # analytical calulation per m
        Fana = (4E-7*np.pi)*dd.ring.I*(np.real(Jcomp)*np.pi*0.226**2.0)/(2.0*np.pi*(dd.ring.z_i - maxForce))
        FanaVec.append(Fana)
        print('Force per meter',Fana)
    
    pl.figure()
    pl.plot(ringRadii, np.array(Fvec),'-bo')
#    pl.plot(ringRadii, np.array(Fvec),'-ko', markerfacecolor='None', markeredgewidth=2)
    pl.hold(True)
    pl.plot(ringRadii, np.array(FanaVec),'-r*')
#    pl.plot(ringRadii, np.array(FanaVec),'-rx', markerfacecolor='None', markeredgewidth=2)
    pl.xlabel('Ring radii [m]')
    pl.ylabel('Magnetic force [N/m]')
    pl.legend(['Model (coaxial rings)','Analytical (parallel wires)'], loc = 'upper right')
    
    Fvec = np.array(Fvec)
    FanaVec = np.array(FanaVec)
    pl.figure()
    pl.semilogy(ringRadii, np.abs(Fvec - FanaVec[:,0]),'-bo')
    pl.xlabel('Ring radii [m]')
    pl.ylabel('Absolute difference between the model and analytical solutions')
    
    pl.figure()
    pl.plot(ringRadii, np.array(indIVec),'-bo')
    pl.hold(True)
    pl.plot(ringRadii, np.array(indImagVec),'-r*')
    pl.xlabel('Ring radii')
    pl.ylabel('Induced current [A]')
    pl.legend(['Re(J) x area','mag(J) x area'], loc = 'upper right')
    
   
    # mesh independence
    meshVec = np.linspace(50,500,11)
    forceResultVec = []
    for n in meshVec:
        thePosVec, theForceVec, Jcomp, powerVec = emlcSim(dd.ring, dd.CuRing, dd.Ar, layers, sections, n, nForce, minForce, maxForce)
        forceResultVec.append(theForceVec[0]*np.pi*0.226**2.0) 
    
    forceResultVec = np.array(forceResultVec)
    pl.figure()
    pl.plot(meshVec,forceResultVec,'-ko')
    pl.xlabel('Number of cells')
    pl.ylabel('Magnetic force [N/m]')

    pl.show()    

        
# Run unittests when module is called as a script
if runUnittests == True: #__name__ == '__main__':
    
    try:
        unittest.main()
    except SystemExit as inst:
        if inst.args[0] is True: # raised by sys.exit(True) when tests failed
            raise
            
