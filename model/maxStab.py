# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 22:31:05 2016

@author: suzanne
"""

import scipy.optimize as op, pylab as pl, numpy as np
import emlc
import defineDesign as dd

def levSimMinGradObj(x, mySample, myAtmosphere, Layers, Sections, Slices):
    # x = [current, freq, r1, r2, r3, r4, r5, r6, r7, r8, z1, z2, z3, z4, z5, z6, z7, z8]
    
    myCoil = dd.sz
    myCoil.x_i = np.array(x[2:10])
    myCoil.z_i = np.array(x[10:18])
  
    myCoil.I   = x[0]
    dd.sz.Ivec = (myCoil.loops[:,4]*myCoil.I) + (np.abs(myCoil.loops[:,4]-1)*(-myCoil.I))
    
    myCoil.f   = x[1]
    
    tupOut = emlc.levSim(myCoil, mySample, myAtmosphere, Layers, Sections, Slices)
    #tupOut = (posVec, forceVec, aSamplePos, F_lift_z, sampleWeight, powerVec, tempVec, T, simTime, aSamplePosGrad)
    return (tupOut[9], tupOut[3], tupOut[4]) # aSamplePosGrad


#    myCoil.k_i = np.ones(np.size(myCoil.z_i))
#    myCoil.k_i[-x[2]:] = 0

Layers   = 15
Sections = 15
Slices   = 15

mySample     = dd.Al
myAtmosphere = dd.Ar

func = lambda X: levSimMinGradObj(X, mySample, myAtmosphere, Layers, Sections, Slices)[0]

#x0 = np.zeros(18)
#x0[0]    = dd.sz.I
#x0[1]    = dd.sz.f
#x0[2:10] = dd.sz.x_i
#x0[10:18] = dd.sz.z_i

x0 = [100.0, 10E3, 0.007, 0.014, 0.021, 0.028, 0.035, 0.042, 0.049, 0.056, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05]

bnds = [(0.0, 1000.0),(10E3, 90E3),(0.007,0.07),(0.007,0.07),(0.007,0.07),(0.007,0.07),(0.007,0.07),(0.007,0.07),(0.007,0.07),(0.007,0.07),(-0.05,0.05),(-0.05,0.05),(-0.05,0.05),(-0.05,0.05),(-0.05,0.05),(-0.05,0.05),(-0.05,0.05),(-0.05,0.05)]

cons = ({'type': 'ineq', 'fun': lambda X: levSimMinGradObj(X, mySample, myAtmosphere, Layers, Sections, Slices)[1] - levSimMinGradObj(X, mySample, myAtmosphere, Layers, Sections, Slices)[2] + 1E-6})
cons = ({'type': 'ineq', 'fun': lambda X: - levSimMinGradObj(X, mySample, myAtmosphere, Layers, Sections, Slices)[1] + levSimMinGradObj(X, mySample, myAtmosphere, Layers, Sections, Slices)[2] + 1E-6})

#func = lambda x: (x[0]**2.0 - 4)**2.0
#x0   = [3.0]
#print(len(x0))
#bnds = [(0.0,10.0)]
#print(len(bnds))

fmin = op.minimize(fun=func, x0=x0, method='COBYLA', bounds=bnds, constraints=cons, options={'ftol': 1e-12})
#fmin = op.minimize(fun=func, x0=x0, method='SLSQP', bounds=bnds, constraints=cons)

print(fmin)
print('END')