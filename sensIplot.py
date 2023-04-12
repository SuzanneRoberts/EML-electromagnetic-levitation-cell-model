import sys
sys.path.append('model')

import defineDesign as dd
import pylab as pl   
from emlcSim import emlcSim
import numpy as np


def sensIplot():

    Layers   = 15
    Sections = 15
    Slices   = 15
    nForce   = 25
    minForce = -0.005 # -0.005 # -0.045 # -0.0023   #
    maxForce =  0.005 # 0.005  #  0.030 # 0.017 #
    
    coilDesign = dd.sz #sz
    material   = dd.Fe #Fe
    atmosphere = dd.Ar
    
    coilDesign.I = (1 - 50/100)*250
    coilDesign.Ivec = (coilDesign.loops[:,4]*coilDesign.I) + (np.abs(coilDesign.loops[:,4]-1)*(-coilDesign.I))
    thePosVecA, theForceVecA, Jcomp, powerVec = emlcSim(coilDesign, material, atmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    coilDesign.I = (1 - 10/100)*250
    coilDesign.Ivec = (coilDesign.loops[:,4]*coilDesign.I) + (np.abs(coilDesign.loops[:,4]-1)*(-coilDesign.I))
    thePosVecB, theForceVecB, Jcomp, powerVec = emlcSim(coilDesign, material, atmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    coilDesign.I = 250
    coilDesign.Ivec = (coilDesign.loops[:,4]*coilDesign.I) + (np.abs(coilDesign.loops[:,4]-1)*(-coilDesign.I))
    thePosVecC, theForceVecC, Jcomp, powerVec = emlcSim(coilDesign, material, atmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)

    coilDesign.I = (1 + 10/100)*250
    coilDesign.Ivec = (coilDesign.loops[:,4]*coilDesign.I) + (np.abs(coilDesign.loops[:,4]-1)*(-coilDesign.I))
    thePosVecD, theForceVecD, Jcomp, powerVec = emlcSim(coilDesign, material, atmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)

    coilDesign.I = (1 + 50/100)*250
    coilDesign.Ivec = (coilDesign.loops[:,4]*coilDesign.I) + (np.abs(coilDesign.loops[:,4]-1)*(-coilDesign.I))
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
    
#    # Plot Szekely coil loop positions
#    pl.plot(  0.005*np.ones(2), np.array([-0.015, 0.015]), '-k')
#    pl.plot( 0.0081*np.ones(2), np.array([-0.015, 0.015]), '-k')
#    pl.plot( -0.005*np.ones(2), np.array([-0.015, 0.015]), '-k')
#    pl.plot(-0.0081*np.ones(2), np.array([-0.015, 0.015]), '-k')
    
    
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
#    pl.xlabel('Sample position along the z-axis (centerline of the coil) [m]. \n Symmetric coil loops are at z = - 0.008, z = - 0.005, z = 0.005, z = 0.008,  with the symmetry plane at z = 0.')
#    pl.xlabel('Sample position along the centerline of the coil  [m]')
    pl.ylabel('Lifting force acting on the sample [N]')
    pl.legend(['50% decrease in coil current',
               '10% decrease in coil current',
               'Coil current = 250 A',
               '10% increase in coil current',
               '50% increase in coil current'], loc='upper right')    
#    pl.legend(['-25% A model with Kermanpur et al. coil design',
#               '-10% A model with Kermanpur et al. coil design',
#               '250 A model with Kermanpur et al. coil design',
#               '+10% A model with Kermanpur et al. coil design',
#               '+25% A model with Kermanpur et al. coil design',
#               'Kermanpur et al. model results with 500 A'], loc='lower right', fontsize=18)
    
    pl.grid(True)               
    pl.tick_params(labelsize=20)
        
    pl.show()


# Function call to run the investigation
sensIplot()