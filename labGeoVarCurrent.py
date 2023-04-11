import sys
sys.path.append('model')

import defineDesign as dd
import pylab as pl   
from emlcSim import emlcSim
from emlc import levSim
import numpy as np
import plotting as pp
import time


def labGeoVarCurrent():

    ticAll = time.time()    
    
    Layers   = 25
    Sections = 25
    Slices   = 25
    
    myCoil = dd.lb
    mySample = dd.Fe
    myAtmosphere = dd.Ar
    
    nForce   = 20
    minForce = np.min(myCoil.z_i) - 2*myCoil.R
    maxForce = np.max(myCoil.z_i) + 2*myCoil.R
    
    currentVec = np.array([100, 150, 200, 250, 300])
    
    mySample.R = 0.005
    
    tic = time.time()
    myCoil.I = currentVec[0]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut1 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
    thePosVec1, theForceVec1, Jcomp1, powerVec1 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    toc = time.time()
    
    myCoil.I = currentVec[1]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut2 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec2, theForceVec2, Jcomp2, powerVec2 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    myCoil.I = currentVec[2]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut3 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec3, theForceVec3, Jcomp3, powerVec3 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    myCoil.I = currentVec[3]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut4 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec4, theForceVec4, Jcomp4, powerVec4 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    myCoil.I = currentVec[4]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut5 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec5, theForceVec5, Jcomp5, powerVec5 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    

    mySample.R = 0.007
    
    myCoil.I = currentVec[0]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut6 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
    thePosVec6, theForceVec6, Jcomp6, powerVec6 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    toc = time.time()
    
    myCoil.I = currentVec[1]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut7 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec7, theForceVec7, Jcomp7, powerVec7 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
   
    myCoil.I = currentVec[2]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut8 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec8, theForceVec8, Jcomp8, powerVec8 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    myCoil.I = currentVec[3]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut9 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec9, theForceVec9, Jcomp9, powerVec9 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
   
    myCoil.I = currentVec[4]
    myCoil.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut10 = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    thePosVec10, theForceVec10, Jcomp10, powerVec10 = emlcSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices, nForce, minForce, maxForce)
    
    
#    pl.figure()
#    pl.plot(tupOut1[0], tupOut1[1], '-bo')
#    pl.plot(tupOut1[2], tupOut1[3], 'rx')
#    pl.plot(tupOut1[0], tupOut1[4]*np.ones(len(tupOut1[0])), 'k')
    
    pp.plotDissertation((np.array(thePosVec1)*1000, theForceVec1, 
#                         np.array(tupOut1[0])*1000, tupOut1[1],
#                         tupOut1[2]*1000, tupOut1[3],
                         np.array(thePosVec2)*1000, theForceVec2,
#                         np.array(tupOut2[0])*1000, tupOut2[1],
#                         tupOut2[2]*1000, tupOut2[3], 
                         np.array(thePosVec3)*1000, theForceVec3,
#                         np.array(tupOut3[0])*1000, tupOut3[1],
#                         tupOut3[2]*1000, tupOut3[3],
                         np.array(thePosVec4)*1000, theForceVec4,
#                         np.array(tupOut4[0])*1000, tupOut4[1],
#                         tupOut4[2]*1000, tupOut4[3], 
                         np.array(thePosVec5)*1000, theForceVec5,
#                         np.array(tupOut5[0])*1000, tupOut5[1],
#                         tupOut5[2]*1000, tupOut5[3],
                         np.array(thePosVec5)*1000, tupOut5[4]*np.ones(len(thePosVec5)) ), 
                        'forcePlotVarCurrent5mmFe', 
                        xString='Sample position along the centerline of the coil [mm]', 
                        yString='Force [N]', 
                        legendList=['Levitation force, '+str(currentVec[0])+'A, '+'Not levitated', 
                                    'Levitation force, '+str(currentVec[1])+'A, '+'Sample @ '+str(round(tupOut2[2][0]*1000,1))+' mm', 
                                    'Levitation force, '+str(currentVec[2])+'A, '+'Sample @ '+str(round(tupOut3[2][0]*1000,1))+' mm',
                                    'Levitation force, '+str(currentVec[3])+'A, '+'Sample @ '+str(round(tupOut4[2][0]*1000,1))+' mm',
                                    'Levitation force, '+str(currentVec[4])+'A, '+'Sample @ '+str(round(tupOut5[2][0]*1000,1))+' mm',
                                    'Sample weight (radius = 5 mm)'], 
                        locString='lower right')
                        
    pp.plotDissertation((np.array(thePosVec6)*1000, theForceVec6, 
#                         np.array(tupOut6[0])*1000, tupOut6[1],
#                         tupOut6[2]*1000, tupOut6[3],
                         np.array(thePosVec7)*1000, theForceVec7,
#                         np.array(tupOut7[0])*1000, tupOut7[1],
#                         tupOut7[2]*1000, tupOut7[3], 
                         np.array(thePosVec8)*1000, theForceVec8,  
#                         np.array(tupOut8[0])*1000, tupOut8[1],
#                         tupOut8[2]*1000, tupOut8[3],
                         np.array(thePosVec9)*1000, theForceVec9,
#                         np.array(tupOut9[0])*1000, tupOut9[1],
#                         tupOut9[2]*1000, tupOut9[3], 
                         np.array(thePosVec10)*1000, theForceVec10, 
#                         np.array(tupOut10[0])*1000, tupOut10[1],
#                         tupOut10[2]*1000, tupOut10[3],
                         np.array(thePosVec10)*1000, tupOut10[4]*np.ones(len(thePosVec10)) ), 
                        'forcePlotVarCurrent7mmFe', 
                        xString='Sample position along the centerline of the coil [mm]', 
                        yString='Force [N]', 
                        legendList=['Levitation force, '+str(currentVec[0])+'A, '+'Not levitated', 
                                    'Levitation force, '+str(currentVec[1])+'A, '+'Sample @ '+str(round(tupOut7[2][0]*1000,1))+' mm', 
                                    'Levitation force, '+str(currentVec[2])+'A, '+'Sample @ '+str(round(tupOut8[2][0]*1000,1))+' mm',
                                    'Levitation force, '+str(currentVec[3])+'A, '+'Sample @ '+str(round(tupOut9[2][0]*1000,1))+' mm',
                                    'Levitation force, '+str(currentVec[4])+'A, '+'Sample @ '+str(round(tupOut10[2][0]*1000,1))+' mm',
                                    'Sample weight (radius = 7 mm)',], 
                        locString='lower right') 
                        
    pp.plotDissertation((currentVec[1:], np.array([ tupOut2[7], tupOut3[7], tupOut4[7], tupOut5[7] ])-273.15,
                         currentVec[1:], np.array([ tupOut7[7], tupOut8[7], tupOut9[7], tupOut10[7] ])-273.15 ), 
                        'tempPlotVarCurrentFe', 
                        xString='Coil Current [A]', 
                        yString='Temperature [$^{\circ}$C]', 
                        legendList=['5mm sample radius', 
                                    '7mm sample radius'], 
                        locString='upper left')

    pl.show()
    
    tocAll = time.time()
    
    print('Time to run one simmulation', toc-tic)
    print('Time to run the full labGeo', tocAll-ticAll)
    
    #meshIndependent(myCoil, mySample, myAtmosphere)
    
    
# Function call to run the investigation
labGeoVarCurrent()
