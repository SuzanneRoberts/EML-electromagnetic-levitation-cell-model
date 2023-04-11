import sys
sys.path.append('model')

import defineDesign as dd
import pylab as pl   
from emlcSim import emlcSim
import numpy as np

def varyMatProps():

    Layers   = 20
    Sections = 20
    Slices   = 20
    nForce   = 1
    minForce =  0.010 # -0.045
    maxForce =  0.010 # 0.030

    scaleFactor = np.array([0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 4.0, 6.0, 8.0, 10])
    force10s = []
    
    for n in scaleFactor:
        
#        #dd.Fe.mu0 = n*4E-7*np.pi
#        dd.Cu.sigma = n*5.96E7      
#        thePosVec10, theForceVec10, Jcomp, powerVec = emlcSim(dd.fj, dd.Cu, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
#        force10s.append(theForceVec10[0])
        
#        dd.Cu.sigma = n*5.96E7
#        dd.Cu.R = 0.006
#        thePosVec12, theForceVec12, Jcomp, powerVec = emlcSim(dd.fj, dd.Cu, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
#        
#        dd.Cu.sigma = n*5.96E7
#        dd.Cu.R = 0.0075
#        thePosVec15, theForceVec15, Jcomp, powerVec = emlcSim(dd.fj, dd.Cu, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
#        
        dd.Cu.sigma = n*5.96E7
        dd.Cu.R = 0.010
        thePosVec20, theForceVec20, Jcomp, powerVec = emlcSim(dd.fj, dd.Cu, dd.Ar, Layers, Sections, Slices, nForce, minForce, maxForce)
        force10s.append(theForceVec20[0])
        
#        pl.figure(101)
#        #pl.plot(thePosVec, theForceVec)
#       
#        #pld.plotWithMoghimiCaseForce(thePosVec, theForceVec)
#        #pld.plotWithKermanpurCaseForce(thePosVec, (theForceVec - dd.Zn.weight)) #/20
#        #pl.legend(['Model (current implementation) divided by 20','Kermanpur et al. model reported results','Kermanpur et al. reported data points'], loc='lower right', fontsize=22)
#        pld.plotWithFrommJehnCaseForce(thePosVec10, theForceVec10, thePosVec12, theForceVec12, thePosVec15, theForceVec15, thePosVec20, theForceVec20)
#        pl.legend(['Fromm & Jehn experimental data, d=20mm',
#               'Model (current implementation), d=20mm',
#               'Fromm & Jehn experimental data, d=15mm',
#               'Model (current implementation), d=15mm',
#               'Fromm & Jehn experimental data, d=12mm',
#               'Model (current implementation), d=12mm',
#               'Fromm & Jehn experimental data, d=10mm',
#               'Model (current implementation), d=10mm'], loc='upper right', fontsize=18)
    pl.figure(102)
    pl.plot(scaleFactor, np.array(force10s), '-bo')
        
    pl.show()
    
    
# Function call to run investigation
varyMatProps()
