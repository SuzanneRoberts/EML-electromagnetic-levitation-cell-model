import defineDesign as dd
import pylab as pl   
from emlcSim import emlcSim
import numpy as np

def plotRootsFigs():

    Layers   = 15
    Sections = 15
    Slices   = 15
    nForce   = 50
    minForce = -0.01
    maxForce =  0.015
    
    myCoil       = dd.loop1dir1
    mySample     = dd.Cu
    myAtmosphere = dd.Ar
    
    thePosVec, theForceVec, Jcomp, powerVec = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    sampleWeight = ((4.0/3.0)*np.pi*(dd.Cu.R**3.0))*dd.Cu.rho*9.81 # W = mg = V(rho)g
    
    # figure
    fig = pl.figure()
    ax = fig.add_subplot(111)
    
    pl.plot(thePosVec, theForceVec, 'k-')
    pl.plot(thePosVec, sampleWeight*np.ones(np.size(thePosVec)), 'k--')
    
    pl.xlabel('Sample position along the z-axis (centerline of the coil) [m]')
    pl.ylabel('Force [N]')
    pl.legend(['Levitation force','Sample weight'], loc = 'lower right')
    
#    # for loop2dir1
#    ax.annotate('coil loop 1', xy=(myCoil.z_i[0], -0.1), xytext=(myCoil.z_i[0], -0.08),
#            arrowprops=dict(facecolor='black', shrink=0.08), fontsize = 16)            
#    ax.annotate('coil loop 2', xy=(myCoil.z_i[1], -0.1), xytext=(myCoil.z_i[1], -0.08),
#            arrowprops=dict(facecolor='black', shrink=0.08), fontsize = 16)
#    ax.annotate('Stable levitation position', xy=(-0.0024, 0.046), xytext=(0.0, 0.07),
#            arrowprops=dict(facecolor='black', shrink=0.08), fontsize = 16)
            
    # for loop1dir1
    ax.annotate('coil loop 1', xy=(myCoil.z_i[0], -0.15), xytext=(myCoil.z_i[0], -0.1),
            arrowprops=dict(facecolor='black', shrink=0.08), fontsize = 16)            
#    ax.annotate('Stable levitation position', xy=(0.009, 0.046), xytext=(0.01, 0.07),
#            arrowprops=dict(facecolor='black', shrink=0.08), fontsize = 16)            
    
    pl.show()

