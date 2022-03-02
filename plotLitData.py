# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 07:52:12 2016

@author: Suzanne
"""

import numpy as np, pylab as pl

def plotWithMoghimiCaseForce(modelPosVec, modelForceVec):
    pl.figure()
    pl.plot(modelPosVec, modelForceVec, '-k', linewidth = 2)
    
    # Compare with Moghimi model
    moghimiModel = np.loadtxt('MoghimiModel/MoghimiFig8.txt')
    moghimiZ     = moghimiModel[:,0]
    moghimiF     = moghimiModel[:,1]  
    
    pl.plot(moghimiZ, moghimiF, 'kd', markersize = 10)
    
    # Compare with Kermanpur model of Moghimi case
    KermanpurModelMo = np.loadtxt('KermanpurModel/KermanpurFig6mo.txt')
    KermanpurMoZ     = KermanpurModelMo[:,0]
    KermanpurMoF     = KermanpurModelMo[:,1]  
    
    pl.plot(KermanpurMoZ, KermanpurMoF, 'b*', markersize = 10)
    
    pl.xlabel('Sample position along the centerline of the coil (center of the bottom loop is zero) [m]', fontsize=24)
    pl.ylabel('Lifting force acting on the sample [N]', fontsize=24)
    pl.legend(['Model (current implementation)','Moghimi et al. model reported results','Kermanpur et al. model reported results'], loc='upper right', fontsize=22)
    
    pl.tick_params(labelsize=20)
    

    
def plotWithKermanpurCaseForce(modelPosVec, modelForceVec):
    pl.figure()
    pl.plot(modelPosVec, modelForceVec, '-k', linewidth = 2)
    
    # Compare with Kermanpur model case
    KermanpurModel = np.loadtxt('KermanpurModel/KermanpurFig8.txt')
    KermanpurZ     = KermanpurModel[:,0]
    KermanpurF     = KermanpurModel[:,1]  
    
    pl.plot(KermanpurZ, KermanpurF, 'b*', markersize = 10)
    
    # Compare with Kermanpur model case data points
    KermanpurModelpnts = np.loadtxt('KermanpurModel/KermanpurFig8pnts.txt')
    KermanpurZpnts     = KermanpurModelpnts[:,0]
    KermanpurFpnts     = KermanpurModelpnts[:,1]  
    
    pl.plot(KermanpurZpnts, KermanpurFpnts, 'ko', markersize = 10)
    
    pl.xlabel('Sample position along the centerline of the coil (center of the second layer of loops is zero) [m]', fontsize=24)
    pl.ylabel('Total force (lifting plus gravitational force) \n acting on the sample [N]', fontsize=24)
    pl.legend(['Model (current implementation)','Kermanpur et al. model reported results','Kermanpur et al. reported data points'], loc='lower right', fontsize=22)
    
    pl.tick_params(labelsize=20)
    
    
def plotWithFrommJehnCaseForce(modelPosVec10, modelForceVec10, modelPosVec12, modelForceVec12, modelPosVec15, modelForceVec15, modelPosVec20, modelForceVec20):
    # load Fromm & Jehn's experimental data from the text file
    fj10cmDynes = np.loadtxt('FrommJehnExperimentalData/FrommJehnFig4-10mmDiam.txt')
    fj12cmDynes = np.loadtxt('FrommJehnExperimentalData/FrommJehnFig4-12mmDiam.txt')
    fj15cmDynes = np.loadtxt('FrommJehnExperimentalData/FrommJehnFig4-15mmDiam.txt')
    fj20cmDynes = np.loadtxt('FrommJehnExperimentalData/FrommJehnFig4-20mmDiam.txt')
    
    
    fj10 = np.zeros(np.shape(fj10cmDynes))
    fj10[:,0] = fj10cmDynes[:,0]*1E-2
    fj10[:,1] = fj10cmDynes[:,1]*1E-5
    
    fj12 = np.zeros(np.shape(fj12cmDynes))
    fj12[:,0] = fj12cmDynes[:,0]*1E-2
    fj12[:,1] = fj12cmDynes[:,1]*1E-5
    
    fj15 = np.zeros(np.shape(fj15cmDynes))
    fj15[:,0] = fj15cmDynes[:,0]*1E-2
    fj15[:,1] = fj15cmDynes[:,1]*1E-5
    
    fj20 = np.zeros(np.shape(fj20cmDynes))
    fj20[:,0] = fj20cmDynes[:,0]*1E-2
    fj20[:,1] = fj20cmDynes[:,1]*1E-5
    
    pl.figure()
    pl.plot(fj20[:,0],fj20[:,1],'kv', markersize = 10)
    pl.plot(modelPosVec20,modelForceVec20,'k', linewidth = 2)
    pl.plot(fj15[:,0],fj15[:,1],'gs', markersize = 10)
    pl.plot(modelPosVec15,modelForceVec15,'g', linewidth = 2)
    pl.plot(fj12[:,0],fj12[:,1],'r*', markersize = 10)
    pl.plot(modelPosVec12,modelForceVec12,'r', linewidth = 2)
    pl.plot(fj10[:,0],fj10[:,1],'bo', markersize = 10)
    pl.plot(modelPosVec10,modelForceVec10,'b', linewidth = 2)
  
    
    #pl.title('Lifting force along the centerline of the coil')
    pl.xlabel('Sample position along the centerline of the coil (center of the single coil loop is zero) [m]', fontsize=24)
    pl.ylabel('Lifting force acting on the sample [N]', fontsize=24)
    pl.legend(['Fromm & Jehn experimental data, d=20mm',
               'Model (current implementation), d=20mm',
               'Fromm & Jehn experimental data, d=15mm',
               'Model (current implementation), d=15mm',
               'Fromm & Jehn experimental data, d=12mm',
               'Model (current implementation), d=12mm',
               'Fromm & Jehn experimental data, d=10mm',
               'Model (current implementation), d=10mm'], loc='upper left', fontsize=18)
               
    pl.tick_params(labelsize=20)
    
    
    
def plotWithFrommJehnCaseForce0(modelPosVec10, modelForceVec10, modelPosVec12, modelForceVec12, modelPosVec15, modelForceVec15, modelPosVec20, modelForceVec20):
    # load Fromm & Jehn's experimental data from the text file
    fj10cmDynes = np.loadtxt('FrommJehnExperimentalData/FrommJehnFig4-10mmDiam.txt')
    fj12cmDynes = np.loadtxt('FrommJehnExperimentalData/FrommJehnFig4-12mmDiam.txt')
    fj15cmDynes = np.loadtxt('FrommJehnExperimentalData/FrommJehnFig4-15mmDiam.txt')
    fj20cmDynes = np.loadtxt('FrommJehnExperimentalData/FrommJehnFig4-20mmDiam.txt')
    
    
    fj10 = np.zeros(np.shape(fj10cmDynes))
    fj10[:,0] = fj10cmDynes[:,0]*1E-2
    fj10[:,1] = fj10cmDynes[:,1]*1E-5
    
    fj12 = np.zeros(np.shape(fj12cmDynes))
    fj12[:,0] = fj12cmDynes[:,0]*1E-2
    fj12[:,1] = fj12cmDynes[:,1]*1E-5
    
    fj15 = np.zeros(np.shape(fj15cmDynes))
    fj15[:,0] = fj15cmDynes[:,0]*1E-2
    fj15[:,1] = fj15cmDynes[:,1]*1E-5
    
    fj20 = np.zeros(np.shape(fj20cmDynes))
    fj20[:,0] = fj20cmDynes[:,0]*1E-2
    fj20[:,1] = fj20cmDynes[:,1]*1E-5
    
#    pl.figure()
#    pl.plot(fj20[:,0],fj20[:,1],'kv', markersize = 10)
#    pl.plot(modelPosVec20,modelForceVec20,'k-', linewidth = 2)
#    pl.plot(fj15[:,0],fj15[:,1],'ks', markersize = 10)
#    pl.plot(modelPosVec15,modelForceVec15,'k--', linewidth = 2)
#    pl.plot(fj12[:,0],fj12[:,1],'k*', markersize = 10)
#    pl.plot(modelPosVec12,modelForceVec12,'k:', linewidth = 2)
#    pl.plot(fj10[:,0],fj10[:,1],'ko', markersize = 10)
#    pl.plot(modelPosVec10,modelForceVec10,'k.', linewidth = 2)

    pl.figure()
    pl.plot(fj20[:,0],fj20[:,1],'bv', markersize = 10)
    pl.plot(modelPosVec20,modelForceVec20,'b-', linewidth = 4)
    pl.plot(fj15[:,0],fj15[:,1],'rs', markersize = 10)
    pl.plot(modelPosVec15,modelForceVec15,'r-', linewidth = 4)
    pl.plot(fj12[:,0],fj12[:,1],'g*', markersize = 10)
    pl.plot(modelPosVec12,modelForceVec12,'g-', linewidth = 4)
    pl.plot(fj10[:,0],fj10[:,1],'ko', markersize = 10)
    pl.plot(modelPosVec10,modelForceVec10,'k-', linewidth = 4)  
    
    #pl.title('Lifting force along the centerline of the coil')
    pl.xlabel('Sample position along the centerline of the coil (center of the single coil loop is zero) [m]', fontsize=24)
    pl.ylabel('Lifting force acting on the sample [N]', fontsize=24)
#    pl.legend(['Fromm & Jehn experiment, d=20mm',
#               'Model, d=20mm',
#               'Fromm & Jehn experiment, d=15mm',
#               'Model, d=15mm',
#               'Fromm & Jehn experiment, d=12mm',
#               'Model, d=12mm',
#               'Fromm & Jehn experiment, d=10mm',
#               'Model, d=10mm'], loc='upper right', fontsize=18)

    pl.legend(['Experiment, d=20mm',
               'Model, d=20mm',
               'Experiment, d=15mm',
               'Model, d=15mm',
               'Experiment, d=12mm',
               'Model, d=12mm',
               'Experiment, d=10mm',
               'Model, d=10mm'], loc='upper right', fontsize=18)
               
    pl.tick_params(labelsize=20)
    
    
    
def plotWithElKaddahSzekelyCaseForce(modelPosVec200, modelForceVec200, modelPosVec250, modelForceVec250, modelPosVec300, modelForceVec300):
    # load El-Kaddah and Szekely model results from the text file
    sz200 = np.loadtxt('SzekelyModel/szFig4force.txt')
    sz250 = np.loadtxt('SzekelyModel/szFig5force.txt')
    sz300 = np.loadtxt('SzekelyModel/szFig6force.txt')
    

    pl.figure()
    pl.plot(sz300[:,0],sz300[:,1],'kv', markersize = 10)
    pl.plot(modelPosVec300,modelForceVec300,'k', linewidth = 2)
    pl.plot(sz250[:,0],sz250[:,1],'r*', markersize = 10)
    pl.plot(modelPosVec250,modelForceVec250,'r', linewidth = 2)
    pl.plot(sz200[:,0],sz200[:,1],'bo', markersize = 10)
    pl.plot(modelPosVec200,modelForceVec200,'b', linewidth = 2)
  
    
    #pl.title('Lifting force along the centerline of the coil')
    pl.xlabel('Sample position along the centerline of the coil (symmetry plane of the symmetrical coil is zero) [m]', fontsize=24)
    pl.ylabel('Lifting force acting on the sample [N]', fontsize=24)
    pl.legend(['El-Kaddah & Szekely model results, I=300A',
               'Model (current implementation), I=300A',
               'El-Kaddah & Szekely model results, I=250A',
               'Model (current implementation), I=250A',
               'El-Kaddah & Szekely model results, I=200A',
               'Model (current implementation), I=200A'], loc='upper right', fontsize=18)
               
    pl.tick_params(labelsize=20)
        
    
    
def plotWithElKaddahSzekelyCasePower(modelPosVec200, modelPowerVec200, modelPosVec250, modelPowerVec250, modelPosVec300, modelPowerVec300):
    # load El-Kaddah and Szekely model results from the text file
    sz200 = np.loadtxt('SzekelyModel/szFig4power.txt')
    sz250 = np.loadtxt('SzekelyModel/szFig5power.txt')
    sz300 = np.loadtxt('SzekelyModel/szFig6power.txt')
    

    pl.figure()
    pl.plot(sz300[:,0],sz300[:,1],'kv', markersize = 10)
    pl.plot(modelPosVec300,modelPowerVec300,'k', linewidth = 2)
    pl.plot(sz250[:,0],sz250[:,1],'r*', markersize = 10)
    pl.plot(modelPosVec250,modelPowerVec250,'r', linewidth = 2)
    pl.plot(sz200[:,0],sz200[:,1],'bo', markersize = 10)
    pl.plot(modelPosVec200,modelPowerVec200,'b', linewidth = 2)
  
    
    #pl.title('Lifting force along the centerline of the coil')
    pl.xlabel('Sample position along the centerline of the coil (symmetry plane of the symmetrical coil is zero) [m]', fontsize=24)
    pl.ylabel('Power absorbed by the sample [W]', fontsize=24)
    pl.legend(['El-Kaddah & Szekely model results, J=300A',
               'Model (current implementation), J=300A',
               'El-Kaddah & Szekely model results, J=250A',
               'Model (current implementation), J=250A',
               'El-Kaddah & Szekely model results, J=200A',
               'Model (current implementation), J=200A'], loc='upper left', fontsize=18)
               
    pl.tick_params(labelsize=20)
    
    
    
def plotWithRoyerOptTemp(iVec, tVec):
    pl.figure()
    pl.plot(iVec, tVec, '-k', linewidth = 2)
    
    # Compare with Moghimi model
    royerOptExp  = np.loadtxt('RoyerTempExpData/royer_opt_exp.txt')
    royerOptExpI = royerOptExp[:,0]
    royerOptExpT = royerOptExp[:,1]  
    
    pl.plot(royerOptExpI, royerOptExpT, 'kd', markersize = 10)
    
    
    pl.xlabel('Coil current [A]', fontsize=24)
    pl.ylabel('Temperature [$^\circ$C]', fontsize=24)
    pl.legend(['Model','Royer et al. experiment'], loc='upper right', fontsize=22)
    
    pl.tick_params(labelsize=20)