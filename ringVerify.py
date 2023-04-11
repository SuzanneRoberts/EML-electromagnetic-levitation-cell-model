import sys
sys.path.append('model')

import defineDesign as dd
import pylab as pl   
from emlcSim import emlcSim
import numpy as np


def ringVerify():

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
    pl.plot(ringRadii, np.array(indImagVec),'-r*')
    pl.xlabel('Ring radii')
    pl.ylabel('Induced current [A]')
    pl.legend(['Re(J) x area','mag(J) x area'], loc = 'upper right')
    
   
    # mesh independence
    meshVec = np.linspace(50,500,11, dtype=int)
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
    
    
# Function call to run investigation
ringVerify()
