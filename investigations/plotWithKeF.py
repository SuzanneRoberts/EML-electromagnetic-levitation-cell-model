import defineDesign as dd
import pylab as pl   
from emlcSim import emlcSim
import numpy as np
import plotLitData as pld

def plotWithKeF():

    Layers   = 20
    Sections = 20
    Slices   = 20
    nForce   = 25
    minForce = -0.045
    maxForce =  0.030
    
    thePosVec, theForceVec, Jcomp, powerVec = emlcSim(dd.ke, dd.Zn, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)

    
    pl.figure(101)

    znVolume  = (4/3)*np.pi*dd.Zn.R**3.0
    znMass    = dd.Zn.rho * znVolume
    znWeight  = znMass*9.81
   
    pld.plotWithKermanpurCaseForce(thePosVec, (np.array(theForceVec) - znWeight))
    
    pld.plotWithKermanpurCaseForce(thePosVec, (np.array(theForceVec) - znWeight)/35)
    pl.legend(['Fromm & Jehn model (current implementation), \n Force divided by 35','Kermanpur et al. model reported results','Kermanpur et al. reported data points'], loc='lower right', fontsize = 22)
    
    pl.show()
