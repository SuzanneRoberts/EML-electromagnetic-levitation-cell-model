import sys
sys.path.append('model')

from model import defineDesign as dd
import pylab as pl   
from model.emlcSim import emlcSim
from model.emlc import levSim
from model import plotLitData as pld
import numpy as np


def compareRoyerTemp():

    Layers   = 25
    Sections = 25
    Slices   = 25
    
    myCoil = dd.roOpt
    mySample = dd.Al
    myAtmosphere = dd.Ar
    
    # Vary the sample emissivity
    baseValue = 0.1
    percentChangeS = 70
    percentChangeL = 110

    mySample.epsilon = (1 + percentChangeL/100)*baseValue
    
    #Ivec = (np.arange(20)+17)*10
    Ivec = (np.arange(20)+17)*10
    Tvec = []
    zvec = []
    for myI in Ivec:
        myCoil.I = myI
        dd.roOpt.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
        tupOut = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
        # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
        Tvec.append(tupOut[7])
        zvec.append(tupOut[2])
  
    pld.plotWithRoyerOptTemp(Ivec, np.array(Tvec)-273.15) 
    
    pl.show()
      
##    myCoil = dd.ro
##    
##    TvecSeed = []
##    for myI in Ivec:
##        myCoil.I = myI
##        dd.roOpt.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
##        tupOut = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
##        # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
##        TvecSeed.append(tupOut[7])
#    
#    pl.figure()
#    pl.plot(Ivec, np.array(Tvec)-273.15, '-bo')
##    pl.plot(Ivec, np.array(TvecSeed)-273.15, '-ro')
#    pl.show()
#    
#    pl.figure()
#    pl.plot(Ivec, np.array(zvec), '-go')
##    pl.plot(Ivec, np.array(TvecSeed)-273.15, '-ro')
#    pl.show()


# Function call to run investigation
compareRoyerTemp()
