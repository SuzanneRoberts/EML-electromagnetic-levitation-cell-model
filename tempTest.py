# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 09:26:52 2016

@author: Suzanne
"""

import numpy as np, pylab as pl
import defineDesign as dd
import plotLitData as pld
import royerPos as rp

dd.sz.I = 200
dd.sz.Ivec = (dd.sz.loops[:,4]*dd.sz.I) + (np.abs(dd.sz.loops[:,4]-1)*(-dd.sz.I))
zvecSz200, FvecSz200, HvecSz200  = rp.forceCurve(dd.Fe,dd.sz,dd.Ar)

dd.sz.I = 250
dd.sz.Ivec = (dd.sz.loops[:,4]*dd.sz.I) + (np.abs(dd.sz.loops[:,4]-1)*(-dd.sz.I))
zvecSz250, FvecSz250, HvecSz250  = rp.forceCurve(dd.Fe,dd.sz,dd.Ar)

dd.sz.I = 300
dd.sz.Ivec = (dd.sz.loops[:,4]*dd.sz.I) + (np.abs(dd.sz.loops[:,4]-1)*(-dd.sz.I))
zvecSz300, FvecSz300, HvecSz300  = rp.forceCurve(dd.Fe,dd.sz,dd.Ar)

pld.plotWithElKaddahSzekelyCaseForce(zvecSz200, FvecSz200, zvecSz250, FvecSz250, zvecSz300, FvecSz300)

pl.show()