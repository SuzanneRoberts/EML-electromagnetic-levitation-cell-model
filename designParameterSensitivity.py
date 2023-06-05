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
    coil_sz      = dd.sz # dd.loop2dir1
    tupOut_sz    = levSim(coil_sz, mySample, myAtmosphere, Layers, Sections, Slices)
    pl.plot(tupOut_sz[0], tupOut_sz[1], '-bo')
    pl.plot(tupOut_sz[2], tupOut_sz[3], 'rx')
    pl.plot(tupOut_sz[0], tupOut_sz[4]*np.ones(len(tupOut_sz[0])), 'k')
    
#    myCoil.I = 130
#    dd.ro.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    
#    tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
    
    print('powerVec', tupOut_sz[5])
    print('tempVec', tupOut_sz[6])
    print('T', tupOut_sz[7])
    print('simTime', tupOut_sz[8])
    
    # Double the number of loops for the Szekely paper coil design
    coil_sz2x    = dd.sz_doubleLoops
    tupOut_sz2x  = levSim(coil_sz2x, mySample, myAtmosphere, Layers, Sections, Slices)
    pl.plot(tupOut_sz2x[0], tupOut_sz2x[1], '-go')
    pl.plot(tupOut_sz2x[2], tupOut_sz2x[3], 'cx')
    pl.plot(tupOut_sz2x[0], tupOut_sz2x[4]*np.ones(len(tupOut_sz2x[0])), 'm')
    
    coil_sz8b4t   = dd.sz_doubleBottomLoops
    tupOut_sz8b4t = levSim(coil_sz8b4t, mySample, myAtmosphere, Layers, Sections, Slices)
    pl.plot(tupOut_sz8b4t[0], tupOut_sz8b4t[1], '--go')
    pl.plot(tupOut_sz8b4t[2], tupOut_sz8b4t[3], 'cx')
    pl.plot(tupOut_sz8b4t[0], tupOut_sz8b4t[4]*np.ones(len(tupOut_sz8b4t[0])), 'm')
    
    coil_sz8b2t   = dd.sz_8b2t
    print(coil_sz8b2t.x_i)
    print(coil_sz8b2t.z_i)
    print(coil_sz8b2t.I)
    tupOut_sz8b2t = levSim(coil_sz8b2t, mySample, myAtmosphere, Layers, Sections, Slices)
    print(tupOut_sz8b2t[7])
    pl.plot(tupOut_sz8b2t[0], tupOut_sz8b2t[1], '-bs')
    pl.plot(tupOut_sz8b2t[2], tupOut_sz8b2t[3], 'cx')
    pl.plot(tupOut_sz8b2t[0], tupOut_sz8b2t[4]*np.ones(len(tupOut_sz8b2t[0])), 'm')
    
    coil_sz8b2t.I = 200.0
    dd.sz_8b2t.Ivec = (coil_sz8b2t.loops[:,4]*coil_sz8b2t.I) + (np.abs(coil_sz8b2t.loops[:,4]-1)*(-coil_sz8b2t.I))
    print(coil_sz8b2t.I)
    tupOut_sz8b2t_I = levSim(coil_sz8b2t, mySample, myAtmosphere, Layers, Sections, Slices)
    pl.plot(tupOut_sz8b2t_I[0], tupOut_sz8b2t_I[1], '--mo')
    pl.plot(tupOut_sz8b2t_I[2], tupOut_sz8b2t_I[3], 'cx')
    pl.plot(tupOut_sz8b2t_I[0], tupOut_sz8b2t_I[4]*np.ones(len(tupOut_sz8b2t_I[0])), 'm')
    
    coil_sz8b2t.I = 400.0
    dd.sz_8b2t.Ivec = (coil_sz8b2t.loops[:,4]*coil_sz8b2t.I) + (np.abs(coil_sz8b2t.loops[:,4]-1)*(-coil_sz8b2t.I))
    print(coil_sz8b2t.I)
    coil_sz8b2t.f = 450.0
    print(coil_sz8b2t.f)
    print(coil_sz8b2t.x_i)
    print(coil_sz8b2t.z_i)
    tupOut_sz8b2t_fI = levSim(coil_sz8b2t, mySample, myAtmosphere, Layers, Sections, Slices)
    print(tupOut_sz8b2t_fI[7])
    pl.plot(tupOut_sz8b2t_fI[0], tupOut_sz8b2t_fI[1], '--cx')
    pl.plot(tupOut_sz8b2t_fI[2], tupOut_sz8b2t_fI[3], 'cx')
    pl.plot(tupOut_sz8b2t_fI[0], tupOut_sz8b2t_fI[4]*np.ones(len(tupOut_sz8b2t_fI[0])), 'm')
    
    coil_sr1b = dd.sr_1b
    tupOut_sr1b = levSim(coil_sr1b, mySample, myAtmosphere, Layers, Sections, Slices)
    pl.plot(tupOut_sr1b[0], tupOut_sr1b[1], '--bx')
    pl.plot(tupOut_sr1b[2], tupOut_sr1b[3], 'rx')
    pl.plot(tupOut_sr1b[0], tupOut_sr1b[4]*np.ones(len(tupOut_sr1b[0])), 'k')

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
