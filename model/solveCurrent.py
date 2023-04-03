# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 16:06:10 2015

@author: Suzanne
"""
import numpy as np, pylab as pl, unittest
from model import mutualInduc as mi
from model import discSample as ds
from model import defineDesign as dd # only to define coils and samples for unittests

# Function to solve for the induced current distribution in the sample
def solveCurrent(nodes1, coil, sample):
   
    # A - full quater of coefficient matrix
    numNodes = len(nodes1)
    A0 = np.zeros((numNodes,numNodes))
    x = nodes1[:,0]
    z = nodes1[:,2]
    
    # Compute effective radii for self-inductance calculations
    effR = np.sqrt(nodes1[:,3]/np.pi) # from area of a circle, effR = sqrt(A/pi)

    for i in range(numNodes):
        # adding mutual inductance terms to A
        A0[:,i] = mi.MIellip(x,z,x[i],z[i],sample.mu0)*nodes1[i,3] # ellip
        # adding self-inductance terms to A
        if effR[i] > x[i]:
            effR[i] = x[i]
            
        A0[i,i] = mi.SIellip(x[i],effR[i],sample.mu0)*nodes1[i,3] #ellip
        
    A = A0*coil.omega*sample.sigma
        
    # B - diagonal quater of coefficient matrix
    b = 2.0*np.pi*nodes1[:,0]
    B = np.diag(b)
        
    # Assemble LHS
    C = np.zeros((2*numNodes, 2*numNodes))
    C[:numNodes,:numNodes] = A
    C[numNodes:,numNodes:] = A
    C[:numNodes,numNodes:] = B
    C[numNodes:,:numNodes] = -B
 
    
    # Assemble RHS
    numLoops = len(coil.loops)
    RHS = np.zeros((2*numNodes, 1))
    for i in range(numLoops):
        RHS[:numNodes,0] += - mi.MIellip(x,z,coil.loops[i,0],coil.loops[i,2],sample.mu0)*np.real(coil.Ivec[i]) # ellip
        RHS[numNodes:,0] += mi.MIellip(x,z,coil.loops[i,0],coil.loops[i,2],sample.mu0)*np.imag(coil.Ivec[i]) # ellip

    RHS *= coil.omega*sample.sigma

        
    # Solve linear system
    J = np.linalg.solve(C,RHS)
        
    reJ   = J[:numNodes,0]
    imJ   = J[numNodes:,0]
    Jcomp = reJ + imJ*1j
  
    return Jcomp, C
    

# Testing
class testSolveCurrent(unittest.TestCase):
    def testSpyMatrix(self):
        testNodes, nodesSph, sourceSample, sourceSampleSph, rPrimeSample, MAGrPRIMEsample, UNITrPRIMEsample, Verts, zeniths, azimuths \
        = ds.discSample(0.005, 20, 20, 20)
        
        x_i = np.array([0.01015, 0.01325, 0.01015, 0.01325, 0.01015, 0.01325, 0.01015, 0.01325]) # m (loop radii)
        z_i = np.array([0.0, 0.0, 0.0031, 0.0031, 0.0131, 0.0131, 0.0162, 0.0162]) # m (relative loop heights)
        k_i = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]) # binary (current direction)

        sz = dd.coil( x_i, z_i, k_i, 
                     450E3, # Hz (frequency of the coil current)
                     200.,  # A (current in coil)
                     0.0031 # m (coil tube radius)
                     )
                     
        Al = dd.sample(0.0934,       # (emissivity of liquid droplet, Royer et al, 2013)
            4252890.,     # S/m (electrical conductivity, Royer et al, 2013)
            4E-7*np.pi,   # H/m (magnetic permeability (free space OR sample? air ~ Al ~ Cu), Royer et al, 2013)
            0.005,        # m (sample radius)
            2702.,        # kg/m^3 (density of Al @ 300K, Cengel)
            660.+273.     # K (liquidus temperature of Al, engineering toolbox)
            ) 
        
        
        testJcomp, testC = solveCurrent(testNodes, sz, Al)
        
        pl.figure()
        pl.spy(testC)
        
        
# Run unittests when module is called as a script
if __name__ == '__main__':

    try:
        unittest.main()
    except SystemExit as inst:
        if inst.args[0] is True: # raised by sys.exit(True) when tests failed
            raise
