# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 14:34:00 2015

@author: Suzanne
"""

import numpy as np, pylab as pl

def Jconvert2cart(Jvec, Ivec, Slices, Coil, NodesPos):
    numNodes = len(Jvec)
    numLoops = len(Ivec)
    Jcart = np.zeros((numNodes*Slices, 3)) + np.zeros((numNodes*Slices, 3))*1j
    Icart = np.zeros((numLoops*Slices, 3)) + np.zeros((numLoops*Slices, 3))*1j
    for n in np.arange(Slices):
        Jcart[n*numNodes:(n*(numNodes)+numNodes), 0] = -np.sin(n*(2.0*np.pi/Slices))*Jvec
        Jcart[n*numNodes:(n*(numNodes)+numNodes), 1] =  np.cos(n*(2.0*np.pi/Slices))*Jvec
        #Jcart[n*numNodes:(n*(numNodes)+numNodes), 0] = -np.sin(n*(2.0*np.pi/Slices))*np.real(Jvec) - np.sin(n*(2.0*np.pi/Slices))*np.imag(Jvec)*1j
        #Jcart[n*numNodes:(n*(numNodes)+numNodes), 1] =  np.cos(n*(2.0*np.pi/Slices))*np.real(Jvec) + np.cos(n*(2.0*np.pi/Slices))*np.imag(Jvec)*1j
        Icart[n*numLoops:(n*(numLoops)+numLoops), 0] = -np.sin(n*(2.0*np.pi/Slices))*Ivec/Coil.loops[:,3]
        Icart[n*numLoops:(n*(numLoops)+numLoops), 1] =  np.cos(n*(2.0*np.pi/Slices))*Ivec/Coil.loops[:,3]
        

#    fig = pl.figure()
#    from mpl_toolkits.mplot3d import Axes3D
#    ax = fig.add_subplot(111, projection='3d')
#    ax.scatter(NodesPos[:,0], NodesPos[:,1], NodesPos[:,2], c=np.abs(np.real(Jvec)) )
#    #pl.legend(['nodesSph','nodes'])
#    pl.axis('equal')
#    pl.show()
#    
#    # first slice
#    n = 0
#    slice1x = Jcart[n*numNodes:(n*(numNodes)+numNodes), 0]
#    slice1y = Jcart[n*numNodes:(n*(numNodes)+numNodes), 1]
#        
#    fig = pl.figure()
#    from mpl_toolkits.mplot3d import Axes3D
#    ax = fig.add_subplot(111, projection='3d')
#    ax.scatter(NodesPos[:,0], NodesPos[:,1], NodesPos[:,2], c=(np.sqrt(np.real(slice1x)**2.0 + np.real(slice1y)**2.0)) )
#    #pl.legend(['nodesSph','nodes'])
#    pl.axis('equal')
#    pl.show()
#    
#    pl.figure()
#    pl.plot(np.abs(np.real(Jvec)))
#    pl.plot(np.sqrt(np.real(slice1x)**2.0 + np.real(slice1y)**2.0))
#    pl.legend(['Jvec','Jcart'])
    
    return Jcart, Icart
    