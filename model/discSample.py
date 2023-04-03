# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 11:45:09 2015

@author: Suzanne
"""

import numpy as np, pylab as pl, unittest

def discSampleNoGrading(R, layers, sections, slices):
    # REQUIRED MODULES: math, numpy as np, pylab as pl

    # Populate nodes
    numBlocks = layers*sections

    layerR       = R/layers
    sectionTheta = np.pi/sections
    slicePhi     = 2*np.pi/slices

    nodesSph   = np.zeros((layers*sections,5))
    nodes      = np.zeros((layers*sections,4))

    count = 0
    for i in np.arange(layers):
        nodesSph[count:count+sections,0] = (i + 0.5)*layerR
        nodesSph[count:count+sections,3] = 0.5*(np.pi*(((i+1)*layerR)**2.0-(i*layerR)**2.0))/sections
        nodesSph[count:count+sections,4] = (((4.0/3.0)*np.pi*(((i+1)*layerR)**3.0 - (i*layerR)**3.0))/sections)/slices  # element volume
        count += sections
        
    zeniths = (np.arange(layers) + 0.5)*layerR
        
    for i in np.arange(sections):
        nodesSph[sections*np.arange(layers)+i,1] = i*sectionTheta + 0.5*sectionTheta
        
    azimuths = sectionTheta*np.arange(sections) + 0.5*sectionTheta
        
    nodes[:,0] = nodesSph[:,0]*np.sin(nodesSph[:,1])
    nodes[:,2] = nodesSph[:,0]*np.cos(nodesSph[:,1])
    nodes[:,3] = nodesSph[:,3]
    
    # Populate sourceSample
    sourceSampleSph = np.zeros((layers*sections*slices,4))
    sourceSample    = np.zeros((layers*sections*slices,4))
    for i in np.arange(slices):
        sourceSampleSph[i*numBlocks:(i*numBlocks)+numBlocks, :] = nodesSph[:,0:4]
        sourceSampleSph[i*numBlocks:(i*numBlocks)+numBlocks, 2] = i*slicePhi
        sourceSampleSph[i*numBlocks:(i*numBlocks)+numBlocks, 3] = nodesSph[:,4]  # element volume
        
    sourceSample[:,0] = sourceSampleSph[:,0]*np.sin(sourceSampleSph[:,1])*np.cos(sourceSampleSph[:,2])
    sourceSample[:,1] = sourceSampleSph[:,0]*np.sin(sourceSampleSph[:,1])*np.sin(sourceSampleSph[:,2])
    sourceSample[:,2] = sourceSampleSph[:,0]*np.cos(sourceSampleSph[:,1])
    sourceSample[:,3] = sourceSampleSph[:,3] # element volume
    
    # compute separation matrix for sample
    rPrimeSample = np.zeros((layers*sections*slices, layers*sections, 3))
    for block in np.arange(numBlocks):
        rPrimeSample[:, block, 0] = nodes[block,0] - sourceSample[:,0]
        rPrimeSample[:, block, 1] = nodes[block,1] - sourceSample[:,1]
        rPrimeSample[:, block, 2] = nodes[block,2] - sourceSample[:,2]
            
    magRprimeSample  = np.sqrt(rPrimeSample[:,:,0]**2.0 + rPrimeSample[:,:,1]**2.0 + rPrimeSample[:,:,2]**2.0) 
    
    unitRprimeSample = np.zeros(np.shape(rPrimeSample))
    for i in np.arange(3):
        unitRprimeSample[:,:,i] = rPrimeSample[:,:,i]/magRprimeSample
    
    return nodes, nodesSph, sourceSample, sourceSampleSph, rPrimeSample, magRprimeSample, unitRprimeSample, zeniths, azimuths


def discSample(R, layers, sections, slices):
    # REQUIRED MODULES: math, numpy as np, pylab as pl

    # Populate nodes
    numBlocks = layers*sections

    sectionTheta = np.pi/sections
    slicePhi     = 2*np.pi/slices

    nodesSph   = np.zeros((layers*sections,5))
    nodes      = np.zeros((layers*sections,4))

    # Add radial mesh grading
    linLayers      = np.linspace(0.0, R, layers+1)
    tempGrading    = linLayers**(1.0/3.0) # Gradingfunction of linLayers: linLayers**(1.0/4.0)
    tempMax        = np.max(tempGrading)
    finalMin       = np.min(linLayers)
    finalMax       = np.max(linLayers)
    gradedLayers   = finalMin + (finalMax - finalMin)*tempGrading/tempMax
    layerThickness = gradedLayers[1:len(gradedLayers)] - gradedLayers[0:(len(gradedLayers)-1)]
    
    # Matrix of vertices for matplotlib.PolyCollection plot
    vertsSph = np.zeros((layers*sections, 4, 2))
    verts    = np.zeros((layers*sections, 4, 2))
    
    count = 0
    for i in np.arange(layers):
        nodesSph[count:count+sections,0] = gradedLayers[i] + 0.5*layerThickness[i]
        nodesSph[count:count+sections,3] = 0.5*(np.pi*(gradedLayers[i+1]**2.0 - gradedLayers[i]**2.0))/sections
        nodesSph[count:count+sections,4] = (((4.0/3.0)*np.pi*(gradedLayers[i+1]**3.0 - gradedLayers[i]**3.0))/sections)/slices
        count += sections

        # Populate radii of vertsSph
        vertsSph[count:count+sections,0,0] = gradedLayers[i]
        vertsSph[count:count+sections,1,0] = gradedLayers[i]
        vertsSph[count:count+sections,2,0] = gradedLayers[i+1]
        vertsSph[count:count+sections,3,0] = gradedLayers[i+1]
        
    zeniths = gradedLayers[0:layers] + 0.5*layerThickness

        
    for i in np.arange(sections):
        nodesSph[sections*np.arange(layers)+i,1] = i*sectionTheta + 0.5*sectionTheta
        
        # Populate theta angles of vertsSph
        vertsSph[sections*np.arange(layers)+i,0,1] = i*sectionTheta
        vertsSph[sections*np.arange(layers)+i,3,1] = i*sectionTheta
        vertsSph[sections*np.arange(layers)+i,1,1] = (i+1)*sectionTheta
        vertsSph[sections*np.arange(layers)+i,2,1] = (i+1)*sectionTheta
        
    azimuths = sectionTheta*np.arange(sections) + 0.5*sectionTheta

                
    nodes[:,0] = nodesSph[:,0]*np.sin(nodesSph[:,1])
    nodes[:,2] = nodesSph[:,0]*np.cos(nodesSph[:,1])
    nodes[:,3] = nodesSph[:,3]
    
    # Convert vertsSph to verts
    verts[:,:,0] = vertsSph[:,:,0]*np.sin(vertsSph[:,:,1])
    verts[:,:,1] = vertsSph[:,:,0]*np.cos(vertsSph[:,:,1])
  
    
    # Populate sourceSample
    sourceSampleSph = np.zeros((layers*sections*slices,4))
    sourceSample    = np.zeros((layers*sections*slices,4))
    for i in np.arange(slices):
        sourceSampleSph[i*numBlocks:(i*numBlocks)+numBlocks, :] = nodesSph[:,0:4]
        sourceSampleSph[i*numBlocks:(i*numBlocks)+numBlocks, 2] = i*slicePhi
        sourceSampleSph[i*numBlocks:(i*numBlocks)+numBlocks, 3] = nodesSph[:,4]
        
    sourceSample[:,0] = sourceSampleSph[:,0]*np.sin(sourceSampleSph[:,1])*np.cos(sourceSampleSph[:,2])
    sourceSample[:,1] = sourceSampleSph[:,0]*np.sin(sourceSampleSph[:,1])*np.sin(sourceSampleSph[:,2])
    sourceSample[:,2] = sourceSampleSph[:,0]*np.cos(sourceSampleSph[:,1])
    sourceSample[:,3] = sourceSampleSph[:,3]
    
    # compute separation matrix for sample
    rPrimeSample = np.zeros((layers*sections*slices, layers*sections, 3))
    for block in np.arange(numBlocks):
        rPrimeSample[:, block, 0] = nodes[block,0] - sourceSample[:,0]
        rPrimeSample[:, block, 1] = nodes[block,1] - sourceSample[:,1]
        rPrimeSample[:, block, 2] = nodes[block,2] - sourceSample[:,2]
            
    magRprimeSample  = np.sqrt(rPrimeSample[:,:,0]**2.0 + rPrimeSample[:,:,1]**2.0 + rPrimeSample[:,:,2]**2.0) 
    
    unitRprimeSample = np.zeros(np.shape(rPrimeSample))
    for i in np.arange(3):
        unitRprimeSample[:,:,i] = rPrimeSample[:,:,i]/magRprimeSample
    
    return nodes, nodesSph, sourceSample, sourceSampleSph, rPrimeSample, magRprimeSample, unitRprimeSample, verts, zeniths, azimuths


def discRing(R, r, slices):
    # REQUIRED MODULES: math, numpy as np, pylab as pl

    # Populate nodes
    numBlocks = 1.0
    
    slicePhi     = 2.0*np.pi/slices
    
    # Set up nodes and nodesSph
    nodesSph   = np.zeros((numBlocks,5))
    nodes      = np.zeros((numBlocks,4))

    nodesSph[0,0] = R
    nodesSph[0,1] = np.pi/2.0
    nodesSph[0,3] = np.pi*r**2.0 # element area
    nodesSph[0,4] = slicePhi*R * nodesSph[0,3] # element volume (arc length * area, where arc length = angle in randians * radius)
    print('nodesSph', nodesSph)
        
    zeniths = nodesSph[:,0]
            
    azimuths = nodesSph[:,1]
    
    verts    = np.zeros((numBlocks, 4, 2))
    verts[0, 0, 0] = R + r
    verts[0, 1, 0] = R
    verts[0, 2, 0] = R - r
    verts[0, 3, 0] = R
    
    verts[0, 0, 1] = 0.0
    verts[0, 1, 1] = r
    verts[0, 2, 1] = 0.0
    verts[0, 3, 1] = - r
        
    nodes[:,0] = nodesSph[:,0]*np.sin(nodesSph[:,1])
    nodes[:,2] = nodesSph[:,0]*np.cos(nodesSph[:,1])
    nodes[:,3] = nodesSph[:,3]
    print('nodes', nodes)
    
    # Populate sourceSample
    sourceSampleSph = np.zeros((numBlocks*slices,4))
    sourceSample    = np.zeros((numBlocks*slices,4))
    for i in np.arange(slices):
        sourceSampleSph[i*numBlocks:(i*numBlocks)+numBlocks, :] = nodesSph[:,0:4]
        sourceSampleSph[i*numBlocks:(i*numBlocks)+numBlocks, 2] = i*slicePhi
        sourceSampleSph[i*numBlocks:(i*numBlocks)+numBlocks, 3] = nodesSph[:,4]  # element volume
        
    sourceSample[:,0] = sourceSampleSph[:,0]*np.sin(sourceSampleSph[:,1])*np.cos(sourceSampleSph[:,2])
    sourceSample[:,1] = sourceSampleSph[:,0]*np.sin(sourceSampleSph[:,1])*np.sin(sourceSampleSph[:,2])
    sourceSample[:,2] = sourceSampleSph[:,0]*np.cos(sourceSampleSph[:,1])
    sourceSample[:,3] = sourceSampleSph[:,3] # element volume
    
    # compute separation matrix for sample
    rPrimeSample = np.zeros((numBlocks*slices, numBlocks, 3))
    for block in np.arange(numBlocks):
        rPrimeSample[:, block, 0] = nodes[block,0] - sourceSample[:,0]
        rPrimeSample[:, block, 1] = nodes[block,1] - sourceSample[:,1]
        rPrimeSample[:, block, 2] = nodes[block,2] - sourceSample[:,2]
            
    magRprimeSample  = np.sqrt(rPrimeSample[:,:,0]**2.0 + rPrimeSample[:,:,1]**2.0 + rPrimeSample[:,:,2]**2.0) 
    
    unitRprimeSample = np.zeros(np.shape(rPrimeSample))
    for i in np.arange(3):
        unitRprimeSample[:,:,i] = rPrimeSample[:,:,i]/magRprimeSample
    
    return nodes, nodesSph, sourceSample, sourceSampleSph, rPrimeSample, magRprimeSample, unitRprimeSample, verts, zeniths, azimuths



class testDiscSample(unittest.TestCase):
    def testAreas(self):
        # TEST - areas: do the areas add up to the correct value?
        R = 0.005
        nodes, nodesSph, sourceSample, sourceSampleSph, rPrimeSample, magRprimeSample, unitRprimeSample, verts, zeniths, azimuths = discSample(R, 5, 8, 3)
        areaTotal = 0.5*(np.pi*R**2.0)
        self.assertTrue(abs(areaTotal - np.sum(nodes[:,3])) < 1E-14)
    def testGeo(self):
        # TEST - coordinates: plot the geometry (only the nodes)
        nodes, nodesSph, sourceSample, sourceSampleSph, rPrimeSample, magRprimeSample, unitRprimeSample, verts, zeniths, azimuths = discSample(0.005, 4, 4, 8)
        pl.figure()
        pl.plot(nodes[:,0], nodes[:,2], 'ko')
        pl.axis('equal')
        pl.show()
        #pl.figure()
        #pl.plot(nodesSph[:,0]*np.sin(nodesSph[:,1]),nodesSph[:,0]*np.cos(nodesSph[:,1]))
        #pl.plot(nodes[:,0],nodes[:,2],'*')
        fig = pl.figure()
        from mpl_toolkits.mplot3d import Axes3D
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(nodes[:,0], nodes[:,1], nodes[:,2], c=np.arange(len(nodes)))
        #pl.legend(['nodesSph','nodes'])
        pl.axis('equal')
        ax.elev = 0
        ax.azim = -90
        pl.show()
    def testSourceSample(self):
        # TEST - sample source coordinates: plot the sample source coordinates (with area indicated by the colour)
        nodes, nodesSph, sourceSample, sourceSampleSph, rPrimeSample, magRprimeSample, unitRprimeSample, verts, zeniths, azimuths = discSample(0.005, 4, 4, 8)        
        fig = pl.figure()
        from mpl_toolkits.mplot3d import Axes3D
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(sourceSample[:,0], sourceSample[:,1], sourceSample[:,2], c=np.arange(len(sourceSample)))
        ax.elev = 0
        ax.azim = 90
        pl.show()
    def testVolumes(self):
        # TEST - volumes: do the volumes add up to the correct value?
        R = 0.005
        nodes, nodesSph, sourceSample, sourceSampleSph, rPrimeSample, magRprimeSample, unitRprimeSample, verts, zeniths, azimuths = discSample(R, 10, 10, 10)
        volTotal = (4.0/3.0)*(np.pi*R**3.0)
        self.assertTrue(abs(volTotal - np.sum(sourceSample[:,3])) < 1E-14)


# Run unittests when module is called as a script
if __name__ == '__main__':
    
    try:
        unittest.main()
    except SystemExit as inst:
        if inst.args[0] is True: # raised by sys.exit(True) when tests failed
            raise