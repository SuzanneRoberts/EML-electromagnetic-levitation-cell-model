import sys
sys.path.append('model')

import defineDesign as dd
import pylab as pl
import numpy as np 
from emlc import levSim


def levPos():

    Layers   = 15
    Sections = 15
    Slices   = 15
    
    myCoil       = dd.sz2 # dd.loop2dir1
    mySample     = dd.Cu
    myAtmosphere = dd.Ar

#    myCoil = dd.ro
#    mySample = dd.Al
#    myAtmosphere = dd.Ar
#    
#    myCoil.I = 130
#    dd.ro.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    tupOut = levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    # tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime)
    
    print('powerVec', tupOut[5])
    print('tempVec', tupOut[6])
    print('T', tupOut[7])
    print('simTime', tupOut[8])
    
    pl.figure()
    pl.plot(tupOut[0], tupOut[1], '-bo')
    pl.plot(tupOut[2], tupOut[3], 'rx')
    pl.plot(tupOut[0], tupOut[4]*np.ones(len(tupOut[0])), 'k')
    pl.show()
    
    
# Function call to run investigation
levPos()
