import sys
sys.path.append('model')

import defineDesign as dd
from emlcSim import emlcSim


def ring4handCompare():

    layers   = 1
    sections = 1
    slices   = 2
    nForce   = 1
    minForce = 0.0
    maxForce = 0.0
    
    thePosVec, theForceVec, Jcomp, powerVec = emlcSim(dd.ring, dd.CuRing, dd.Ar, layers, sections, slices, nForce, minForce, maxForce)
    
    
# Function call to run investigation
ring4handCompare()
