import numpy as np
import time
from model import discSample as ds
from model.samplePos import samplePos
from model.energyBalance import sampleTemp
from model import plotLitData as pld


# EMLC simmulation function
def emlcSim(coil, sample, atmosphere, layers, sections, slices, nForce, minForce, maxForce):
    
    tic = time.time()

    nodes, nodesSph, sourceSample, sourceSampleSph, rPrimeSample, MAGrPRIMEsample, UNITrPRIMEsample, Verts, zeniths, azimuths \
     = ds.discSample(sample.R, layers, sections, slices)
     
#    nodes, nodesSph, sourceSample, sourceSampleSph, rPrimeSample, MAGrPRIMEsample, UNITrPRIMEsample, Verts, zeniths, azimuths \
#     = ds.discRing(sample.R, 0.226, slices)

    
    # Find sample levitation position
    ####
    posVec   = np.linspace(minForce, maxForce, nForce)
    forceVec = np.zeros((nForce,1))
    powerVec = np.zeros((nForce,1))
    
    for p, aSamplePos in enumerate(posVec):
        nodesPos, Jcomp, Jcart, Jmag, Bfield, forceField, F_lift_z, vertsPos = samplePos(aSamplePos, nodes, layers, sections, slices, coil, sample, sourceSample, rPrimeSample, MAGrPRIMEsample, UNITrPRIMEsample, Verts)
        forceVec[p,0] = F_lift_z
                
        # Power calculation (cartesian coordinates)
        specPcomp = (Jcart[:(layers*sections),0]*np.conj(Jcart[:(layers*sections),0])) + (Jcart[:(layers*sections),1]*np.conj(Jcart[:(layers*sections),1])) + (Jcart[:(layers*sections),2]*np.conj(Jcart[:(layers*sections),2]))
#        specP     = 0.5*np.real(specPcomp)/sample.sigma
        specP     = np.real(specPcomp)/sample.sigma
        P         = specP*(nodesPos[:,3]*2.0*np.pi*nodesPos[:,0]) # loop volume (not element volume)
        #P         = (nodesPos[:,3]*2.0*np.pi*nodesPos[:,0]) # loop volume - check volume
        Ptotal    = np.sum(P)
                
        powerVec[p,0] = Ptotal
    
        #    # power calculation in in the spherical coordinate system (can be used as a unittest later, spherical coordinate system result should be the same as the cartesian result)
        #    specPold  = 0.5*np.real(Jcomp*Jcomp.conj())/sample.sigma
        #    Pold      = specPold*(nodesPos[:,3]*2.0*np.pi*nodesPos[:,0])
        #    Ptotalold = np.sum(Pold)
        #    print('Total power absorbed',Ptotalold)
    ####
    
    #print('powerVec',powerVec)    
    
#    pl.figure()
#    pl.plot(posVec, forceVec, 'k')
#    pl.plot(posVec, forceVec, 'o')
#    
#    
#    pl.figure()
#    pl.plot(posVec, powerVec, 'k')
#    pl.plot(posVec, powerVec, 'o')
#    pl.xlabel('Distance from the center of the coil [m] (center of the coil is 0)')
#    pl.ylabel('Power absorbed [W]')
    
    
    # Temperature calculation
    T, Qrad, Qconv = sampleTemp(Ptotal, sample, atmosphere)
   
    toc = time.time()
    simTime = toc - tic
    
   
    return posVec, forceVec, Jcomp, powerVec #, P, Qconv, Qrad, T, Bfield, simTime
#    return posVec, forceField[:,2], Jcomp, powerVec #, P, Qconv, Qrad, T, Bfield, simTime

