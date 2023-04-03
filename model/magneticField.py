# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 12:14:57 2015

@author: Suzanne
"""
import numpy as np, time

from model import mutualInduc as mi

def selfAddsZero(aVector):
    nanVec             = np.isnan(aVector)
    notNanVec          = np.abs(nanVec-1)
    alsoNotInfVec      = aVector*notNanVec
    vecWithZeroForSelf = np.nan_to_num(alsoNotInfVec)
    return vecWithZeroForSelf
    

# Magnetic flux density calculation
def magFieldCalc(Jcart,Icart,sourceSample,rPRIMEsample,MAGrPRIMEsample,UNITrPRIMEsample,sourceCoil,rPRIMEcoils,MAGrPRIMEcoils,UNITrPRIMEcoils,Sample):
    (rs, cs, ds)  = np.shape(rPRIMEsample)
    (rc, cc, dc)  = np.shape(rPRIMEcoils)
    row, numNodes = np.shape(MAGrPRIMEsample)
    nSlices       = row/numNodes
        
    # convert J to cartesian coordinates    
    
    B = np.zeros((numNodes, 3)) + np.zeros((numNodes, 3))*1j
    tic = time.time() # time the cross product calculation
    for f in np.arange(numNodes):
        integral = np.zeros(3) + np.zeros(3)*1j # initialise integral for each new field point
        
        sampleCrosses = np.cross(Jcart[:,:],UNITrPRIMEsample[:,f,:])
        coilCrosses   = np.cross(Icart[:,:],UNITrPRIMEcoils[:,f,:])
        
        # cross product test
        #sampleCrosses_phi = np.sum(np.arctan2(np.real(sampleCrosses[:,1]),np.real(sampleCrosses[:,0])))
        #print('np.real(sampleCrosses[:,0])',np.real(sampleCrosses[:,0]))
        
        for i in np.arange(3):
            integral[i] = np.sum(selfAddsZero(sampleCrosses[:,i]* sourceSample[:,3] / (MAGrPRIMEsample[:,f])**2.0)) \
            + np.sum(selfAddsZero(coilCrosses[:,i] * sourceCoil[:,3] / (MAGrPRIMEcoils[:,f])**2.0))
            
                  
        # save result in B component vectors
        B[f,:] = integral
    toc = time.time()
    timeVectorized = toc - tic
    #print('timeVectorized:', timeVectorized)
        
            
    B *= (Sample.mu0/(4.0*np.pi))
                    
    return B


# Magnetic flux density calculation - not vectorized
def magFieldCalcOld(Jcart,Icart,sourceSample,rPRIMEsample,MAGrPRIMEsample,UNITrPRIMEsample,sourceCoil,rPRIMEcoils,MAGrPRIMEcoils,UNITrPRIMEcoils,Sample):
    (rs, cs, ds) = np.shape(rPRIMEsample)
    (rc, cc, dc) = np.shape(rPRIMEcoils)
    row, numNodes = np.shape(MAGrPRIMEsample)
    nSlices       = row/numNodes
    
    # convert J to cartesian coordinates    
    
    B = np.zeros((numNodes, 3)) + np.zeros((numNodes, 3))*1j
    tic = time.time() # time the cross product calculation
    for f in np.arange(numNodes):
        integral = np.zeros(3) + np.zeros(3)*1j # initialise integral for each new field point
        
        # integrate over the current sources in the sample
        for s in np.arange(rs):
            if f != s:
                integral += np.cross(Jcart[s,:],UNITrPRIMEsample[s,f,:]) * sourceSample[s,3] / (MAGrPRIMEsample[s,f])**2.0
                    
        # integrate over the current sources in the coil
        for s in np.arange(rc):
            integral += np.cross(Icart[s,:],UNITrPRIMEcoils[s,f,:]) * sourceCoil[s,3] / (MAGrPRIMEcoils[s,f])**2.0
        
        # save result in B component vectors
        B[f,:] = integral
    toc = time.time()
    timeNotVectorized = toc - tic
    print('timeNotVectorized:', timeNotVectorized)
        
            
    B *= (Sample.mu0/(4.0*np.pi))
    
    return B
