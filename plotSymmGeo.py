import sys
sys.path.append('model')

import defineDesign as dd
import pylab as pl   
from emlcSim import emlcSim
from emlc import meshIndependent


def plotSymmGeo():

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
    
    
# Function call to run investigation
plotSymmGeo()   
