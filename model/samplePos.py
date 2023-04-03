# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 16:45:30 2015

@author: Suzanne
"""

import numpy as np, time

import coilsRprime as cr
from solveCurrent import solveCurrent
from Jcart import Jconvert2cart
from magneticField import magFieldCalc

def samplePos(aSamplePos, nodes, layers, sections, slices, coil, sample, sourceSample, rPrimeSample, MAGrPRIMEsample, UNITrPRIMEsample, verts):
    
    nodesPos, sourceSamplePos, rPrimeCoils, magRprimeCoils, unitRprimeCoils, sourceCoil, vertsPos \
     = cr.coilsRPrime(aSamplePos, coil.loops, coil.loopsSph, nodes, sourceSample, layers, sections, slices, verts)
    
    # Solve for induced current desity
    Jcomp, C = solveCurrent(nodesPos, coil, sample)
        
    Jmag  = np.sqrt(np.real(Jcomp)**2.0 + np.imag(Jcomp)**2.0)
#    print('Total current induced:',sum(np.real(Jcomp)*nodesPos[:,3]))
    
    
    # Magnetic flux density
    Jcart, Icart = Jconvert2cart(Jcomp, coil.Ivec, slices, coil, nodesPos)
    
    tic = time.time()
    Bfield = magFieldCalc(Jcart,Icart,sourceSample,rPrimeSample,MAGrPRIMEsample,UNITrPRIMEsample,sourceCoil,rPrimeCoils,magRprimeCoils,unitRprimeCoils,sample)
    toc = time.time()
    
   
    # Lifting force calculation
#    forceField = 0.5*np.real(np.cross(Jcart[:(layers*sections),:], np.conj(Bfield)))
    forceField = np.real(np.cross(Jcart[:(layers*sections),:], np.conj(Bfield)))
#    F_lift_x   = np.sum(forceField[:,0]*(nodesPos[:,3]*2.0*np.pi*nodesPos[:,0])) # loop volume not element volume - or is it?
#    F_lift_y   = np.sum(forceField[:,1]*(nodesPos[:,3]*2.0*np.pi*nodesPos[:,0]))
    F_lift_z   = np.sum(forceField[:,2]*(nodesPos[:,3]*2.0*np.pi*nodesPos[:,0]))
    
    # Checks
#    print('F_lift_x', F_lift_x)
#    print('F_lift_y', F_lift_y)
#    print('F_lift_z', F_lift_z)
#    sampleVolume  = (4/3)*np.pi*sample.R**3.0
#    sampleVolume2 = np.sum(sourceSample[:,3])
#    print('The two methods to calculate the sample volume should be the same:', sampleVolume, sampleVolume2)
#    sampleMass    = sample.rho * sampleVolume
#    sampleWeight  = sampleMass*9.81
#    print('weight', sampleWeight)

    
    # Force cross product check: F_lift_phi should be zero
    #F_lift_phi = (np.arctan2(forceField[:,1],forceField[:,0]))
    #print('F_lift_phi',F_lift_phi)
    
    return nodesPos, Jcomp, Jcart, Jmag, Bfield, forceField, F_lift_z, vertsPos
