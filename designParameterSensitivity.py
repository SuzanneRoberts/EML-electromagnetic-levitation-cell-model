import sys
sys.path.append('model')

import defineDesign as dd
import pylab as pl
import numpy as np 
from emlc import levSim


def designParameterSensitivity():

### Comparing different coil designs, all with a
### copper sample and
### argon atmosphere; and 
### using the same 15x15x15 discretisation in all cases.


# Instantiating a sample and an atmosphere
    mySample     = dd.Cu
    myAtmosphere = dd.Ar

# Discretizing the sample
    Layers   = 15
    Sections = 15
    Slices   = 15
    
    
# Instantiating and plotting the different coil designs

    # Szekely paper coil design
    coilSz      = dd.sz # dd.loop2dir1
    tupOutSz    = levSim(coilSz, mySample, myAtmosphere, Layers, Sections, Slices)
    pl.plot(tupOutSz[0], tupOutSz[1], '-bo')
    pl.plot(tupOutSz[2], tupOutSz[3], 'rx')
    pl.plot(tupOutSz[0], tupOutSz[4]*np.ones(len(tupOutSz[0])), 'k')
    
#    myCoil.I = 130
#    dd.ro.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    
#    tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
    
    print('powerVec', tupOutSz[5])
    print('tempVec', tupOutSz[6])
    print('T', tupOutSz[7])
    print('simTime', tupOutSz[8])
    

    myCoil2 = dd.sz_doubleLoops
    tupOut2 = levSim(myCoil2, mySample, myAtmosphere, Layers, Sections, Slices)
    pl.plot(tupOut2[0], tupOut2[1], '-go')
    pl.plot(tupOut2[2], tupOut2[3], 'cx')
    pl.plot(tupOut2[0], tupOut2[4]*np.ones(len(tupOut2[0])), 'm')
    
    myCoil3 = dd.sz_doubleBottomLoops
    tupOut3 = levSim(myCoil3, mySample, myAtmosphere, Layers, Sections, Slices)
    pl.plot(tupOut3[0], tupOut3[1], '--go')
    pl.plot(tupOut3[2], tupOut3[3], 'cx')
    pl.plot(tupOut3[0], tupOut3[4]*np.ones(len(tupOut3[0])), 'm')
    
    myCoil4 = dd.sz_8b2t
    print(myCoil4.I)
    tupOut4 = levSim(myCoil4, mySample, myAtmosphere, Layers, Sections, Slices)
    print(tupOut4[7])
    pl.plot(tupOut4[0], tupOut4[1], '-bs')
    pl.plot(tupOut4[2], tupOut4[3], 'cx')
    pl.plot(tupOut4[0], tupOut4[4]*np.ones(len(tupOut4[0])), 'm')
    
    myCoil4.I = 200.0
    print(myCoil4.I)
    tupOut5 = levSim(myCoil4, mySample, myAtmosphere, Layers, Sections, Slices)
    pl.plot(tupOut5[0], tupOut5[1], '--mo')
    pl.plot(tupOut5[2], tupOut5[3], 'cx')
    pl.plot(tupOut5[0], tupOut5[4]*np.ones(len(tupOut5[0])), 'm')
    
    myCoil4.I = 400.0
    print(myCoil4.I)
    myCoil4.f = 450.0
    print(myCoil4.f)
    tupOut6 = levSim(myCoil4, mySample, myAtmosphere, Layers, Sections, Slices)
    print(tupOut6[7])
    pl.plot(tupOut6[0], tupOut6[1], '--cx')
    pl.plot(tupOut6[2], tupOut6[3], 'cx')
    pl.plot(tupOut6[0], tupOut6[4]*np.ones(len(tupOut5[0])), 'm')
    
    myCoil7 = dd.sr_1b
    tupOut7 = levSim(myCoil7, mySample, myAtmosphere, Layers, Sections, Slices)
    pl.plot(tupOut7[0], tupOut7[1], '--bx')
    pl.plot(tupOut7[2], tupOut7[3], 'rx')
    pl.plot(tupOut7[0], tupOut7[4]*np.ones(len(tupOut7[0])), 'k')

# tube radius
# size of the gap between the different current direction loops
# adding loops to the side (as opposed to vertically)
# paper #1 new title: EMLC design based on mathematical modelling

# paper #2: 
# max. stable levitation region (i.e. aSamplePos - Fmax) and / or
# max. range of sample weights
# (the "and" case will be maximizing the intergral under the Fz - weight curve)
   

    pl.show()
    
    
# Function call to run investigation
designParameterSensitivity()
