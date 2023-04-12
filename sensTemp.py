import sys
sys.path.append('model')

import defineDesign as dd
import pylab as pl
import numpy as np
from emlc import levSim


def sensTemp():

    Layers   = 10 #30
    Sections = 10 #30
    Slices   = 10 #30
    
    myCoil = dd.roOpt
    mySample = dd.Al
    myAtmosphere = dd.Ar
    
    sensParameter = 'emissivity'
#    sensParameter = 'electricalConductivity'
#    sensParameter = 'magneticPermeability'
#    sensParameter = 'sampleRadius'
#    sensParameter = 'thermalConductivity'
#    sensParameter = 'kinematicViscosity'
#    sensParameter = 'position'
    
    # Vary the parameter
    epsilonBaseValue = 0.1
    sigmaBaseValue   = 1E7
    muBaseValue      = 4E-7*np.pi
    rBaseValue       = 0.0031
    kBaseValue       = 152E-3
    etaBaseValue      = 199E-7
    percentChangeS = 10
    percentChangeL = 50
    positionChangeS = 0.004
    positionChangeL = 0.008

    if sensParameter == 'emissivity':
        mySample.epsilon = (1 - percentChangeL/100)*epsilonBaseValue
    if sensParameter == 'electricalConductivity':
        mySample.sigma = (1 - percentChangeL/100)*sigmaBaseValue
    if sensParameter == 'magneticPermeability':
        mySample.mu0 = (1 - percentChangeL/100)*muBaseValue
    if sensParameter == 'sampleRadius':
        mySample.R = (1 - percentChangeL/100)*rBaseValue
    if sensParameter == 'thermalConductivity':
        myAtmosphere.k = (1 - percentChangeL/100)*kBaseValue
    if sensParameter == 'kinematicViscosity':
        myAtmosphere.eta = (1 - percentChangeL/100)*etaBaseValue
    if sensParameter == 'position':
        mySample.posVar = - positionChangeL
    Ivec = (np.arange(20)+17)*10
    TvecA = []
    zvecA = []
    for myI in Ivec:
        myCoil.I = myI
        dd.roOpt.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
        tupOut = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
        # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
        TvecA.append(tupOut[7])
        zvecA.append(tupOut[2])
    
    if sensParameter == 'emissivity':    
        mySample.epsilon = (1 - percentChangeS/100)*epsilonBaseValue
    if sensParameter == 'electricalConductivity':    
        mySample.sigma = (1 - percentChangeS/100)*sigmaBaseValue
    if sensParameter == 'magneticPermeability':
        mySample.mu0 = (1 - percentChangeS/100)*muBaseValue
    if sensParameter == 'sampleRadius':
        mySample.R = (1 - percentChangeS/100)*rBaseValue
    if sensParameter == 'thermalConductivity':
        myAtmosphere.k = (1 - percentChangeS/100)*kBaseValue
    if sensParameter == 'kinematicViscosity':
        myAtmosphere.eta = (1 - percentChangeS/100)*etaBaseValue
    if sensParameter == 'position':
        mySample.posVar = - positionChangeS
    Ivec = (np.arange(20)+17)*10
    TvecB = []
    zvecB = []
    for myI in Ivec:
        myCoil.I = myI
        dd.roOpt.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
        tupOut = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
        # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
        TvecB.append(tupOut[7])
        zvecB.append(tupOut[2])

    if sensParameter == 'emissivity':        
        mySample.epsilon = epsilonBaseValue
    if sensParameter == 'electricalConductivity':        
        mySample.sigma = sigmaBaseValue
    if sensParameter == 'magneticPermeability':
        mySample.mu0 = muBaseValue
    if sensParameter == 'sampleRadius':
        mySample.R = rBaseValue
    if sensParameter == 'thermalConductivity':
        myAtmosphere.k = kBaseValue
    if sensParameter == 'kinematicViscosity':
        myAtmosphere.eta = etaBaseValue
    if sensParameter == 'position':
        mySample.posVar = 0.0
    Ivec = (np.arange(20)+17)*10
    TvecC = []
    zvecC = []
    for myI in Ivec:
        myCoil.I = myI
        dd.roOpt.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
        tupOut = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
        # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
        TvecC.append(tupOut[7])
        zvecC.append(tupOut[2])

    if sensParameter == 'emissivity':
        mySample.epsilon = (1 + percentChangeS/100)*epsilonBaseValue
    if sensParameter == 'electricalConductivity':
        mySample.sigma = (1 + percentChangeS/100)*sigmaBaseValue
    if sensParameter == 'magneticPermeability':
        mySample.mu0 = (1 + percentChangeS/100)*muBaseValue
    if sensParameter == 'sampleRadius':
        mySample.R = (1 + percentChangeS/100)*rBaseValue
    if sensParameter == 'thermalConductivity':
        myAtmosphere.k = (1 + percentChangeS/100)*kBaseValue
    if sensParameter == 'kinematicViscosity':
        myAtmosphere.eta = (1 + percentChangeS/100)*etaBaseValue
    if sensParameter == 'position':
        mySample.posVar = positionChangeS
    #Ivec = (np.arange(20)+17)*10
    Ivec = (np.arange(20)+17)*10
    TvecD = []
    zvecD = []
    for myI in Ivec:
        myCoil.I = myI
        dd.roOpt.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
        tupOut = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
        # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
        TvecD.append(tupOut[7])
        zvecD.append(tupOut[2])
        
    if sensParameter == 'emissivity':
        mySample.epsilon = (1 + percentChangeL/100)*epsilonBaseValue
    if sensParameter == 'electricalConductivity':
        mySample.sigma = (1 + percentChangeL/100)*sigmaBaseValue
    if sensParameter == 'magneticPermeability':
        mySample.mu0 = (1 + percentChangeL/100)*muBaseValue
    if sensParameter == 'sampleRadius':
        mySample.R = (1 + percentChangeL/100)*rBaseValue
    if sensParameter == 'thermalConductivity':
        myAtmosphere.k = (1 + percentChangeL/100)*kBaseValue
    if sensParameter == 'kinematicViscosity':
        myAtmosphere.eta = (1 + percentChangeL/100)*etaBaseValue
    if sensParameter == 'position':
        mySample.posVar = positionChangeL
    Ivec = (np.arange(20)+17)*10
    TvecE = []
    zvecE = []
    for myI in Ivec:
        myCoil.I = myI
        dd.roOpt.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
        tupOut = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
        # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
        TvecE.append(tupOut[7])
        zvecE.append(tupOut[2])
  
  
    pl.figure()
    
#    pl.plot(Ivec, np.array(TvecA), '--mo')
#    pl.plot(Ivec, np.array(TvecB), '--r^')
#    pl.plot(Ivec, np.array(TvecC), '-g*')
#    pl.plot(Ivec, np.array(TvecD), '--bv')
#    pl.plot(Ivec, np.array(TvecE), '--cs')
    
    pl.plot(Ivec, np.array(TvecA), '--ko')
    pl.plot(Ivec, np.array(TvecB), '--k^')
    pl.plot(Ivec, np.array(TvecC), '-k*')
    pl.plot(Ivec, np.array(TvecD), '--kv')
    pl.plot(Ivec, np.array(TvecE), '--ks')

    pl.xlabel('Coil current  [A]')
    pl.ylabel('Temperature [$^\circ$C]')
    
    if sensParameter == 'emissivity':
        pl.legend(['50% decrease in sample emissivity',
                   '10% decrease in sample emissivity',
                   'Emissivity = 0.1',
                   '10% increase in sample emissivity',
                   '50% increase in sample emissivity'], loc='upper right')    
    if sensParameter == 'electricalConductivity':
        pl.legend([str(percentChangeL) + '% decrease in electrical conductivity',
                   str(percentChangeS) + '% decrease in electrical conductivity',
                   'Sample conductivity = %.1E S/m' %sigmaBaseValue,
                   str(percentChangeS) + '% increase in electrical conductivity',
                   str(percentChangeL) + '% increase in electrical conductivity'], loc='upper right') 
        pl.ylim([800,1600])    
    if sensParameter == 'magneticPermeability':
        pl.legend([str(percentChangeL) + '% decrease in magnetic permeability',
                   str(percentChangeS) + '% decrease in magnetic permeability',
                   'Magnetic permeability = %.1E H/m' %muBaseValue,
                   str(percentChangeS) + '% increase in magnetic permeability',
                   str(percentChangeL) + '% increase in magnetic permeability'], loc='upper right')
        pl.ylim([1000,2200])
    if sensParameter == 'sampleRadius':    
        pl.legend([str(percentChangeL) + '% decrease in sample radius',
                   str(percentChangeS) + '% decrease in sample radius',
                   'Sample radius = ' + str(rBaseValue) + ' m',
                   str(percentChangeS) + '% increase in sample radius',
                   str(percentChangeL) + '% increase in sample radius'], loc='upper right')  
        pl.ylim([600,1600])
    if sensParameter == 'thermalConductivity':
        pl.legend([str(percentChangeL) + '% decrease in thermal conductivity',
                   str(percentChangeS) + '% decrease in thermal conductivity',
                   'Atmosphere thermal conductivity = %.1E W/mK' %sigmaBaseValue,
                   str(percentChangeS) + '% increase in thermal conductivity',
                   str(percentChangeL) + '% increase in thermal conductivity'], loc='upper right')
        pl.ylim([1000,1900])
    if sensParameter == 'kinematicViscosity':
        pl.legend([str(percentChangeL) + '% decrease in kinematic viscosity',
                   str(percentChangeS) + '% decrease in kinematic viscosity',
                   'Atmosphere kinematic viscosity = %.1E m$^2$/s' %etaBaseValue,
                   str(percentChangeS) + '% increase in kinematic viscosity',
                   str(percentChangeL) + '% increase in kinematic viscosity'], loc='upper right')        
    if sensParameter == 'position':
        pl.legend([str(positionChangeL*1000) + 'mm decrease in levitation position',
                   str(positionChangeS*1000) + 'mm decrease in levitation position',
                   'Actual levitation position',
                   str(positionChangeS*1000) + 'mm increase in levitation position',
                   str(positionChangeL*1000) + 'mm increase in levitation position'], loc='upper left')
        pl.ylim([0,5000])
    
    pl.grid(True)               
    pl.tick_params(labelsize=20)
    
    pl.show()
    
    
# Function call to run investigation
sensTemp()
