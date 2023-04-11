import sys
sys.path.append('model')

import defineDesign as dd
import pylab as pl   
from emlcSim import emlcSim
import numpy as np

def sensPropPlot():

    Layers   = 15
    Sections = 15
    Slices   = 15
    nForce   = 25
    minForce = -0.005 # -0.005 # -0.045 # -0.0023   #
    maxForce = 0.005 # 0.005  #  0.030 # 0.017 #
    
    coilDesign = dd.sz #mo #sz
    material   = dd.Fe #Ni #Fe
    atmosphere = dd.Ar  
    
    #varProp = material.sigma - does not work to change the parameter in one place only -> Python pointer error?
    baseValue = 4E-7*np.pi #1E7 #7870 #0.0031 #4E-7*np.pi
    percentChangeS = 10
    percentChangeL = 50

    material.mu0 = (1 - percentChangeL/100)*baseValue
    thePosVecA, theForceVecA, Jcomp, powerVec = emlcSim(coilDesign, material, atmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    material.mu0 = (1 - percentChangeS/100)*baseValue
    thePosVecB, theForceVecB, Jcomp, powerVec = emlcSim(coilDesign, material, atmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    material.mu0 = baseValue
    thePosVecC, theForceVecC, Jcomp, powerVec = emlcSim(coilDesign, material, atmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)

    material.mu0 = (1 + percentChangeS/100)*baseValue
    thePosVecD, theForceVecD, Jcomp, powerVec = emlcSim(coilDesign, material, atmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)

    material.mu0 = (1 + percentChangeL/100)*baseValue
    thePosVecE, theForceVecE, Jcomp, powerVec = emlcSim(coilDesign, material, atmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    pl.figure()   
#    pl.plot(thePosVecA, (np.array(theForceVecA) - dd.Zn.weight), '-r*')
#    pl.plot(thePosVecB, (np.array(theForceVecB) - dd.Zn.weight), '-b^')
#    pl.plot(thePosVecC, (np.array(theForceVecC) - dd.Zn.weight), '-go')
#    pl.plot(thePosVecD, (np.array(theForceVecD) - dd.Zn.weight), '-mv')
#    pl.plot(thePosVecE, (np.array(theForceVecE) - dd.Zn.weight), '-cs')
    
    pl.plot(thePosVecA, np.array(theForceVecA), '--mo')
    pl.plot(thePosVecB, np.array(theForceVecB), '--r^')
    pl.plot(thePosVecC, np.array(theForceVecC), '-g*')
    pl.plot(thePosVecD, np.array(theForceVecD), '--bv')
    pl.plot(thePosVecE, np.array(theForceVecE), '--cs')
    
#    pl.plot(thePosVecA, np.array(theForceVecA), '--ko')
#    pl.plot(thePosVecB, np.array(theForceVecB), '--k^')
#    pl.plot(thePosVecC, np.array(theForceVecC), '-k*')
#    pl.plot(thePosVecD, np.array(theForceVecD), '--kv')
#    pl.plot(thePosVecE, np.array(theForceVecE), '--ks')
    
    
#    # Compare with Fromm & Jehn's experimental data
#    fj10cmDynes = np.loadtxt('FrommJehnExperimentalData/FrommJehnFig4-10mmDiam.txt')
#    
#    fj10 = np.zeros(np.shape(fj10cmDynes))
#    fj10[:,0] = fj10cmDynes[:,0]*1E-2
#    fj10[:,1] = fj10cmDynes[:,1]*1E-5
#    
#    pl.plot(fj10[:,0],fj10[:,1],'kd', markersize = 10)
    
    
#    # Compare with El-Kaddah and Szekely model results
#    sz200 = np.loadtxt('SzekelyModel/szFig4force.txt')
# 
#    pl.plot(sz200[:,0],sz200[:,1],'kd', markersize = 10)
    
    
#    # Compare with Moghimi model
#    moghimiModel = np.loadtxt('MoghimiModel/MoghimiFig8.txt')
#    moghimiZ     = moghimiModel[:,0]
#    moghimiF     = moghimiModel[:,1]  
#    
#    pl.plot(moghimiZ, moghimiF, 'kd', markersize = 10)
    
    
#    # Compare with Kermanpur model case
#    KermanpurModel = np.loadtxt('KermanpurModel/KermanpurFig8.txt')
#    KermanpurZ     = KermanpurModel[:,0]
#    KermanpurF     = KermanpurModel[:,1]  
#    
#    pl.plot(KermanpurZ, KermanpurF, 'kd', markersize = 10)
    
    
    pl.xlabel('Sample position along the centerline of the coil (symmetric coil symmetry plane at zero)  [m]')
    #pl.xlabel('Sample position along the centerline of the coil  [m]')
    pl.ylabel('Lifting force acting on the sample [N]')
#    pl.legend(['25% decrease in coil current',
#               '10% decrease in coil current',
#               'Coil current = 250 A',
#               '10% increase in coil current',
#               '25% increase in coil current'], loc='upper right')    
    
    pl.legend([str(percentChangeL) + '% decrease in magnetic permeability',
               str(percentChangeS) + '% decrease in magnetic permeability',
               'Magnetic permeability = %.1E H/m' %baseValue,
               str(percentChangeS) + '% increase in magnetic permeability',
               str(percentChangeL) + '% increase in magnetic permeability'], loc='upper right')
    
#    pl.legend([str(percentChangeL) + '% decrease in electrical conductivity',
#               str(percentChangeS) + '% decrease in electrical conductivity',
#               'Sample conductivity = %.1E S/m' %baseValue,
#               str(percentChangeS) + '% increase in electrical conductivity',
#               str(percentChangeL) + '% increase in electrical conductivity'], loc='upper right')    
               
#    pl.legend([str(percentChangeL) + '% decrease in sample radius',
#               str(percentChangeS) + '% decrease in sample radius',
#               'Sample radius = ' + str(baseValue) + ' m',
#               str(percentChangeS) + '% increase in sample radius',
#               str(percentChangeL) + '% increase in sample radius'], loc='upper right')                 

#    pl.legend(['-25% A model with Kermanpur et al. coil design',
#               '-10% A model with Kermanpur et al. coil design',
#               '250 A model with Kermanpur et al. coil design',
#               '+10% A model with Kermanpur et al. coil design',
#               '+25% A model with Kermanpur et al. coil design',
#               'Kermanpur et al. model results with 500 A'], loc='lower right')
    
    pl.grid(True)               
    pl.tick_params(labelsize=20)
        
    pl.show()


# function call to run the investigation
sensPropPlot()
