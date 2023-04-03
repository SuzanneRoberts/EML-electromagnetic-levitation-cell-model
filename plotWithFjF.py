import sys
sys.path.append('model')

from model import defineDesign as dd
import pylab as pl   
from model.emlcSim import emlcSim
from model import plotLitData as pld


def plotWithFjF():

    Layers   = 20
    Sections = 20
    Slices   = 20
    nForce   = 25
    minForce = -0.045
    maxForce =  0.03
    
    baseValue = 1.256637061435917e-006 #1.256637061435917e-006
    percentChangeS = 10
    percentChangeL = 15


    dd.Cu.mu0 = (1 + percentChangeL/100)*baseValue

    thePosVec10, theForceVec10, Jcomp, powerVec = emlcSim(dd.fj, dd.Cu, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
    dd.Cu.R = 0.006
    thePosVec12, theForceVec12, Jcomp, powerVec = emlcSim(dd.fj, dd.Cu, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
    dd.Cu.R = 0.0075
    thePosVec15, theForceVec15, Jcomp, powerVec = emlcSim(dd.fj, dd.Cu, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
    dd.Cu.R = 0.010
    thePosVec20, theForceVec20, Jcomp, powerVec = emlcSim(dd.fj, dd.Cu, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    pld.plotWithFrommJehnCaseForce(thePosVec10, theForceVec10, thePosVec12, theForceVec12, thePosVec15, theForceVec15, thePosVec20, theForceVec20)
     
    
    dd.Cu.mu0 = baseValue
    
    thePosVec10, theForceVec10, Jcomp, powerVec = emlcSim(dd.fj, dd.Cu, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
    dd.Cu.R = 0.006
    thePosVec12, theForceVec12, Jcomp, powerVec = emlcSim(dd.fj, dd.Cu, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
    dd.Cu.R = 0.0075
    thePosVec15, theForceVec15, Jcomp, powerVec = emlcSim(dd.fj, dd.Cu, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
    dd.Cu.R = 0.010
    thePosVec20, theForceVec20, Jcomp, powerVec = emlcSim(dd.fj, dd.Cu, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    pld.plotWithFrommJehnCaseForce(thePosVec10, theForceVec10, thePosVec12, theForceVec12, thePosVec15, theForceVec15, thePosVec20, theForceVec20)
    

    pl.figure()   
    
    pl.show()
    
    
# Function call to run investigation  
plotWithFjF()
