import sys
sys.path.append('model')

import defineDesign as dd
import pylab as pl   
from emlcSim import emlcSim
import numpy as np


def meshIndependence():

    # mesh independence study of the force curve
    meshDensity = [10, 15, 20, 25, 30, 35]
    #meshDensity = [8, 9, 10, 11, 12, 13]
    
    forceIndVecFj = []
    forceIndVecSz = []
    forceIndVecMo = []
    # forceIndVecKe = []
    
    powerAbsVecFj = []
    powerAbsVecSz = []
    powerAbsVecMo = []
    # powerAbsVecKe = []
    
    # initialize old vectors
    theForceVecFj = []
    theForceVecSz = []
    theForceVecMo = []
    
    powerVecFj = []
    powerVecSz = []
    powerVecMo = []
    
    # initialize old vectors
    oldForceVecFj = []
    oldForceVecSz = []
    oldForceVecMo = []
    
    oldpowerVecFj = []
    oldpowerVecSz = []
    oldpowerVecMo = []
    
    for idx, n in enumerate(meshDensity):
                
        oldForceVecFj = theForceVecFj
        oldForceVecSz = theForceVecSz
        oldForceVecMo = theForceVecMo
        
        oldpowerVecFj = powerVecFj
        oldpowerVecSz = powerVecSz
        oldpowerVecMo = powerVecMo        
        
        Layers   = n
        Sections = n
        Slices   = n
        nForce   = 17
        minForce = -0.04
        maxForce =  0.03
        
        thePosVecFj, theForceVecFj, JcompFj, powerVecFj = emlcSim(dd.fj, dd.Cu, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
        thePosVecSz, theForceVecSz, JcompSz, powerVecSz = emlcSim(dd.sz, dd.Fe, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
        thePosVecMo, theForceVecMo, JcompMo, powerVecMo = emlcSim(dd.mo, dd.Ni, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
        # thePosVecKe, theForceVecKe, JcompKe, powerVecKe = emlcSim(dd.ke, dd.Zn, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
        
        forceIndVecFj.append(np.max(theForceVecFj))
        forceIndVecSz.append(np.max(theForceVecSz))
        forceIndVecMo.append(np.max(theForceVecMo))
        # forceIndVecKe.append(np.max(theForceVecKe))
        
        powerAbsVecFj.append(np.max(powerVecFj))
        powerAbsVecSz.append(np.max(powerVecSz))
        powerAbsVecMo.append(np.max(powerVecMo))
        # powerAbsVecKe.append(powerVecKe[8])
        
        pl.figure(100)
        pl.plot(thePosVecFj, theForceVecFj, 'b', linewidth = 1)
        pl.plot(thePosVecSz, theForceVecSz, 'r', linewidth = 1)
        pl.plot(thePosVecMo, theForceVecMo, 'g', linewidth = 1)
        # pl.plot(thePosVecKe, theForceVecKe)
                
        pl.figure(101)
        pl.plot(thePosVecFj, powerVecFj, 'b', linewidth = 1)
        pl.plot(thePosVecSz, powerVecSz, 'r', linewidth = 1)
        pl.plot(thePosVecMo, powerVecMo, 'g', linewidth = 1)
        # pl.plot(thePosVecKe, powerVecKe)
        
        
        if idx > 0:
            pl.figure(102)
            pl.plot(thePosVecFj, np.abs(theForceVecFj - oldForceVecFj), 'b', linewidth = 1)
            pl.plot(thePosVecSz, np.abs(theForceVecSz - oldForceVecSz), 'r', linewidth = 1)
            pl.plot(thePosVecMo, np.abs(theForceVecMo - oldForceVecMo), 'g', linewidth = 1)
            # pl.plot(thePosVecKe, theForceVecKe)
                        
            pl.figure(103)
            pl.plot(thePosVecFj, np.abs(powerVecFj - oldpowerVecFj), 'b', linewidth = 1)
            pl.plot(thePosVecSz, np.abs(powerVecSz - oldpowerVecSz), 'r', linewidth = 1)
            pl.plot(thePosVecMo, np.abs(powerVecMo - oldpowerVecMo), 'g', linewidth = 1)
            # pl.plot(thePosVecKe, powerVecKe)
            
    # convert lists to np.arrays
    meshDensity   = np.array(meshDensity)
    forceIndVecFj = np.array(forceIndVecFj)
    forceIndVecSz = np.array(forceIndVecSz)
    forceIndVecMo = np.array(forceIndVecMo)
    # forceIndVecKe = np.array(forceIndVecKe)

    pl.figure(100)
    pl.xlabel('Sample position along the centerline of the coil [m]', fontsize=24)
    pl.ylabel('Lifting force acting on the sample [N]', fontsize=24)
    pl.legend(['Fromm & Jehn','El-Kaddah & Szekely','Moghimi et al.'], loc='upper right', fontsize = 22)
    # pl.legend([str(meshDensity[0]), str(meshDensity[1]), str(meshDensity[2]), str(meshDensity[3]), str(meshDensity[4]), str(meshDensity[5])])

    
    pl.figure(101)
    pl.xlabel('Sample position along the centerline of the coil [m]', fontsize=24)
    pl.ylabel('Power absorbed by the sample [W]', fontsize=24)
    pl.legend(['Fromm & Jehn','El-Kaddah & Szekely','Moghimi et al.'], loc='upper right', fontsize = 22)
    # pl.legend([str(meshDensity[0]), str(meshDensity[1]), str(meshDensity[2]), str(meshDensity[3]), str(meshDensity[4]), str(meshDensity[5])])
    pl.show()
    
    
    # lifting force mesh independence plots
    pl.figure()
    pl.plot(meshDensity**3.0, forceIndVecFj, '-b*', markersize = 10, linewidth = 2)
    pl.plot(meshDensity**3.0, forceIndVecSz, '-ro', markersize = 10, linewidth = 2)
    pl.plot(meshDensity**3.0, forceIndVecMo, '-g^', markersize = 10, linewidth = 2)
    # pl.plot(meshDensity**3.0, forceIndVecKe, '-ks', markersize = 10)
    pl.xlabel('Number of cells', fontsize = 22)
    pl.ylabel('Maximum lifting force [N]', fontsize = 22)
    pl.legend(['Fromm & Jehn','El-Kaddah & Szekely','Moghimi et al.'], loc='right', fontsize = 22)
    
    pl.figure()
    pl.plot(meshDensity**3.0, np.abs(forceIndVecFj - forceIndVecFj[-1])*100.0/forceIndVecFj[-1], '-b*', markersize = 10, linewidth = 2)
    pl.plot(meshDensity**3.0, np.abs(forceIndVecSz - forceIndVecSz[-1])*100.0/forceIndVecSz[-1], '-ro', markersize = 10, linewidth = 2)
    pl.plot(meshDensity**3.0, np.abs(forceIndVecMo - forceIndVecMo[-1])*100.0/forceIndVecMo[-1], '-g^', markersize = 10, linewidth = 2)
    # pl.plot(meshDensity**3.0, np.abs(forceIndVecKe - forceIndVecKe[-1])*100.0/forceIndVecKe[-1], '-ks', markersize = 10)
    pl.xlabel('Number of cells', fontsize = 22)
    pl.ylabel('Absolute difference in maximum lifting force \n relative to case with the most cells [%]', fontsize = 22)
    pl.legend(['Fromm & Jehn','El-Kaddah & Szekely','Moghimi et al.'], loc='upper right', fontsize = 22)
    
    # power absorbed mesh independence plots
    pl.figure()
    pl.plot(meshDensity**3.0, powerAbsVecFj, '-b*', markersize = 10, linewidth = 2)
    pl.plot(meshDensity**3.0, powerAbsVecSz, '-ro', markersize = 10, linewidth = 2)
    pl.plot(meshDensity**3.0, powerAbsVecMo, '-g^', markersize = 10, linewidth = 2)
    # pl.plot(meshDensity**3.0, powerAbsVecKe, '-ks', markersize = 10)
    pl.xlabel('Number of cells', fontsize = 22)
    pl.ylabel('Maximum power absorbed [W]', fontsize = 22)
    pl.legend(['Fromm & Jehn','El-Kaddah & Szekely','Moghimi et al.'], loc='right', fontsize = 22)
    
    pl.figure()
    pl.plot(meshDensity**3.0, np.abs(powerAbsVecFj - powerAbsVecFj[-1])*100.0/powerAbsVecFj[-1], '-b*', markersize = 10, linewidth = 2)
    pl.plot(meshDensity**3.0, np.abs(powerAbsVecSz - powerAbsVecSz[-1])*100.0/powerAbsVecSz[-1], '-ro', markersize = 10, linewidth = 2)
    pl.plot(meshDensity**3.0, np.abs(powerAbsVecMo - powerAbsVecMo[-1])*100.0/powerAbsVecMo[-1], '-g^', markersize = 10, linewidth = 2)
    # pl.plot(meshDensity**3.0, np.abs(powerAbsVecKe - powerAbsVecKe[-1])*100.0/powerAbsVecKe[-1], '-ks', markersize = 10)
    pl.xlabel('Number of cells', fontsize = 22)
    pl.ylabel('Absolute difference in maximum power absorbed \n relative to case with the most cells [%]', fontsize = 22)
    pl.legend(['Fromm & Jehn','El-Kaddah & Szekely','Moghimi et al.'], loc='upper right', fontsize = 22)
    
#    # induced current mesh independence plots
#    pl.figure()
#    pl.plot(meshDensity**3.0, np.sum(np.real(JcompFj)), '-b*')
#    pl.plot(meshDensity**3.0, np.sum(np.real(JcompSz)), '-ro')
#    pl.plot(meshDensity**3.0, np.sum(np.real(JcompMo)), '-g^')
#    pl.plot(meshDensity**3.0, np.sum(np.real(JcompKe)), '-ks')
#    pl.xlabel('Number of cells', fontsize = 22)
#    pl.ylabel('Total real current induced [A]', fontsize = 22)
#    pl.legend(['Fromm & Jehn','El-Kaddah & Szekely','Moghimi','Kermanpur'], loc='upper right', fontsize = 22)
   
    pl.show()   
    
    
# Function call to run the investigation
meshIndependence()
