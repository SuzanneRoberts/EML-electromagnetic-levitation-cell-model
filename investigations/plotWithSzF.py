import defineDesign as dd
import pylab as pl   
from emlcSim import emlcSim
import numpy as np
import plotLitData as pld

def plotWithSzF():

    Layers   = 20
    Sections = 20
    Slices   = 20
    nForce   = 25
    minForce = -0.005
    maxForce =  0.005
    
    dd.sz.I = 300
    dd.sz.Ivec = (dd.sz.loops[:,4]*dd.sz.I) + (np.abs(dd.sz.loops[:,4]-1)*(-dd.sz.I))
    thePosVec300, theForceVec300, Jcomp, powerVec300 = emlcSim(dd.sz, dd.Fe, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
    dd.sz.I = 250
    dd.sz.Ivec = (dd.sz.loops[:,4]*dd.sz.I) + (np.abs(dd.sz.loops[:,4]-1)*(-dd.sz.I))
    thePosVec250, theForceVec250, Jcomp, powerVec250 = emlcSim(dd.sz, dd.Fe, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
    dd.sz.I = 200
    dd.sz.Ivec = (dd.sz.loops[:,4]*dd.sz.I) + (np.abs(dd.sz.loops[:,4]-1)*(-dd.sz.I))
    thePosVec200, theForceVec200, Jcomp, powerVec200 = emlcSim(dd.sz, dd.Fe, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)

    pl.figure()
    pld.plotWithElKaddahSzekelyCaseForce(thePosVec200, theForceVec200, thePosVec250, theForceVec250, thePosVec300, theForceVec300)
    
    pl.figure()
    pl.plot(thePosVec200, powerVec200)
    pl.plot(thePosVec250, powerVec250)
    pl.plot(thePosVec300, powerVec300)
    pl.xlabel('Sample position along the centerline of the coil (the symmetry plane of the symmetrical coil is at zero) [m]')
    pl.ylabel('Power absorbed by the sample [W]')
    pl.legend(['I = 200','I = 250','I = 300'], loc = 'upper center')
    
    pl.show()
