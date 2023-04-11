import sys
sys.path.append('model')

import defineDesign as dd
import pylab as pl   
from emlcSim import emlcSim


def szCenter():

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
        pl.show()
                
#    pl.figure(101)
#    pl.legend([str(meshDensity[0]), str(meshDensity[1]), str(meshDensity[2])])


# Function call to run investigation
szCenter()
