import sys
sys.path.append('model')

import defineDesign as dd
import pylab as pl   
from emlcSim import emlcSim
import numpy as np
import plotLitData as pld

def plotWithMoF():

    Layers   = 20
    Sections = 20
    Slices   = 20
    nForce   = 25
    minForce = -0.0
    maxForce =  0.017
    
    thePosVec, theForceVec, Jcomp, powerVec = emlcSim(dd.mo, dd.Ni, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
        
    pld.plotWithMoghimiCaseForce(thePosVec, theForceVec)
    
    pl.show()
    
    
# Function call to run investigation
plotWithMoF()
