import sys
sys.path.append('model')

import defineDesign as dd
import pylab as pl   
from emlcSim import emlcSim


def fjCase():

    layers   = 25
    sections = 25
    slices   = 25
    nForce   = 5
    minForce = 0.0
    maxForce = 0.01
    
    thePosVec, theForceVec, Jcomp, powerVec = emlcSim(dd.fj, dd.Cu, dd.Ar, layers, sections, slices, nForce, minForce, maxForce)
    
    pl.figure()
    pl.plot(thePosVec, theForceVec)
    pl.show()
    
    
# Function call to run investigation
fjCase()
