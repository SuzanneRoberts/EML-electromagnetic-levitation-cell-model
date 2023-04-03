# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 09:05:57 2016

@author: Suzanne
"""

import numpy as np, pylab as pl
import defineDesign as dd
import plotLitData as pld

def forceCurve(s, c, f):
    # change opposite current direction to -1 (from 0)
    # c.k_i += (c.k_i - 1) <- the problem with this is that each time the function is called the
    # the c.k_i attribute is changed (and the new values is saved in c.k_i)
    # Therefore (I think) never update an attribute in a function. If you want to do it, write a method for it.

    Idir  = c.k_i + (c.k_i - 1)
    
    phi   = s.R*np.sqrt(np.pi * c.f * s.mu0 * s.sigma)
    
    gterm = (np.sinh(2*phi) - np.sin(2*phi)) / ((np.sinh(phi))**2.0 + (np.sin(phi))**2.0)
    g     = 1 - (3/(4*phi))*gterm
    
    zmin =  np.min(c.z_i) #  -0.045 #  -0.01 #
    zmax =  np.max(c.z_i) + 4*c.R #  0.03 #  0.01 #
    zvec = np.linspace(zmin,zmax,(np.ceil(zmax-zmin)/c.R))
    Fvec = []
    Hvec = []
    for z in zvec:
        hterm1 = np.sum( (Idir*c.x_i**2.0) / (c.x_i**2.0 + (z - c.z_i)**2.0)**(3.0/2.0) )
        hterm2 = np.sum( (Idir*c.x_i**2.0*(z - c.z_i)) / (c.x_i**2.0 + (z - c.z_i)**2.0)**(5.0/2.0) )
        h      = hterm1 * hterm2
        
        F = (3.0/2.0) * np.pi * s.mu0 * (c.I**2.0) * (s.R**3.0) * g * h
        
        Fvec.append(F)
        
        H = 0.5 * c.I * hterm1
        Hvec.append(H)
        
    return zvec, Fvec, Hvec


plotWithFj    = 0
plotWithSz    = 0
plotWithKe    = 0
plotWithMo    = 0
plotRootsFigs = 1


if plotWithFj == True:
    # Compute model force curves
    zvecfj10, Fvecfj10, Hvecfj10  = forceCurve(dd.Cu,dd.fj,dd.Ar)
    dd.Cu.R = 0.006
    zvecfj12, Fvecfj12, Hvecfj12  = forceCurve(dd.Cu,dd.fj,dd.Ar)
    dd.Cu.R = 0.0075
    zvecfj15, Fvecfj15, Hvecfj15  = forceCurve(dd.Cu,dd.fj,dd.Ar)
    dd.Cu.R = 0.010
    zvecfj20, Fvecfj20, Hvecfj20  = forceCurve(dd.Cu,dd.fj,dd.Ar)
    
    pld.plotWithFrommJehnCaseForce(zvecfj10, np.array(Fvecfj10), zvecfj12, np.array(Fvecfj12), zvecfj15, np.array(Fvecfj15), zvecfj20, np.array(Fvecfj20))
    
    #-----------
    pl.figure()
    #-----------
    pl.plot(zvecfj10,np.array(Hvecfj10))
    pl.title('Magnetic field strength along the centerline of the coil')


if plotWithSz == True:
    # Compute model force curves
    dd.sz.I = 200
    dd.sz.Ivec = (dd.sz.loops[:,4]*dd.sz.I) + (np.abs(dd.sz.loops[:,4]-1)*(-dd.sz.I))
    zvecSz200, FvecSz200, HvecSz200  = forceCurve(dd.Fe,dd.sz,dd.Ar)
    
    dd.sz.I = 250
    dd.sz.Ivec = (dd.sz.loops[:,4]*dd.sz.I) + (np.abs(dd.sz.loops[:,4]-1)*(-dd.sz.I))
    zvecSz250, FvecSz250, HvecSz250  = forceCurve(dd.Fe,dd.sz,dd.Ar)
    
    dd.sz.I = 300
    dd.sz.Ivec = (dd.sz.loops[:,4]*dd.sz.I) + (np.abs(dd.sz.loops[:,4]-1)*(-dd.sz.I))
    zvecSz300, FvecSz300, HvecSz300  = forceCurve(dd.Fe,dd.sz,dd.Ar)
    
    pld.plotWithElKaddahSzekelyCaseForce(zvecSz200, FvecSz200, zvecSz250, FvecSz250, zvecSz300, FvecSz300)
    
    #-----------
    pl.figure()
    #-----------
    pl.plot(zvecSz200,np.array(HvecSz200))
    pl.title('Magnetic field strength along the centerline of the coil')


if plotWithMo == True:
    # Compute model force curve
    zvecMo, FvecMo, HvecMo  = forceCurve(dd.Ni,dd.mo,dd.Ar)
        
    pld.plotWithMoghimiCaseForce(zvecMo, np.array(FvecMo))
    
    pld.plotWithMoghimiCaseForce(zvecMo, np.array(FvecMo)/275)
    pl.legend(['Fromm & Jehn model (current implementation), \n Force divided by 275','Moghimi et al. model reported results','Kermanpur et al. model reported results'], loc='upper right', fontsize = 22)
        

if plotWithKe == True:
    # Compute model force curve
    zvecKe, FvecKe, HvecKe  = forceCurve(dd.Zn,dd.ke,dd.Ar)

    pld.plotWithKermanpurCaseForce(zvecKe, (np.array(FvecKe) - dd.Zn.weight))
    
    pld.plotWithKermanpurCaseForce(zvecKe, (np.array(FvecKe) - dd.Zn.weight)/30000)
    pl.legend(['Fromm & Jehn model (current implementation), \n Force divided by 30000','Kermanpur et al. model reported results','Kermanpur et al. reported data points'], loc='lower right', fontsize = 22)


if plotRootsFigs == True:
    # Compute model force curve for a single loop coil
    zvec1loop1dir, Fvec1loop1dir, Hvec1loop1dir  = forceCurve(dd.Cu,dd.loop1dir1,dd.Ar)
    zvec2loop1dir, Fvec2loop1dir, Hvec2loop1dir  = forceCurve(dd.Cu,dd.loop2dir1,dd.Ar)
    zvec2loop2dir, Fvec2loop2dir, Hvec2loop2dir  = forceCurve(dd.Cu,dd.loop2dir2,dd.Ar)
    
    #-----------
    pl.figure()
    #-----------
    pl.plot(zvec2loop1dir, np.array(Fvec2loop1dir))
    pl.title('Single loop coil force vector')

pl.show()