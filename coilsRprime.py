# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 15:08:33 2015

@author: Suzanne
"""
import numpy as np, unittest, pylab as pl, copy

# Function to add the coils as source points to the separation vector matrix
def coilsRPrime(samplePos1, loops, loopsSph, nodes, sourceSample, layers, sections, slices, verts):
    
    # Update nodes for current sample position
    nodesPos              = copy.deepcopy(nodes)
    sourceSamplePos       = copy.deepcopy(sourceSample)
    vertsPos              = copy.deepcopy(verts)
    nodesPos[:,2]        += samplePos1   
    sourceSamplePos[:,2] += samplePos1
    vertsPos[:,:,1]         += samplePos1
    
    slicePhi = 2*np.pi/slices
    
    (numLoops,c1) = np.shape(loops)
    numBlocks     = layers*sections
        
    sourceCoilSph = np.zeros((slices*numLoops,4))
    sourceCoil    = np.zeros((slices*numLoops,4))
    for i in np.arange(slices):
        sourceCoilSph[i*numLoops:(i*numLoops)+numLoops, :] = loopsSph[:,0:4]
        sourceCoilSph[i*numLoops:(i*numLoops)+numLoops, 2] = i*slicePhi
        sourceCoilSph[i*numLoops:(i*numLoops)+numLoops, 3] = loops[:,3]*slicePhi*loops[:,0]

    sourceCoil[:,0] = sourceCoilSph[:,0]*np.sin(sourceCoilSph[:,1])*np.cos(sourceCoilSph[:,2])
    sourceCoil[:,1] = sourceCoilSph[:,0]*np.sin(sourceCoilSph[:,1])*np.sin(sourceCoilSph[:,2])
    sourceCoil[:,2] = sourceCoilSph[:,0]*np.cos(sourceCoilSph[:,1])
    sourceCoil[:,3] = sourceCoilSph[:,3]
    
    rPrimeCoils   = np.zeros((slices*numLoops,numBlocks,4)) 
    
    # compute separation matrix for coil
    rPrimeCoils = np.zeros((numLoops*slices, layers*sections, 3))
    for block in np.arange(numBlocks):
        rPrimeCoils[:, block, 0] = nodesPos[block,0] - sourceCoil[:,0] ##nodesPos was nodes!
        rPrimeCoils[:, block, 1] = nodesPos[block,1] - sourceCoil[:,1] ##nodesPos was nodes!
        rPrimeCoils[:, block, 2] = nodesPos[block,2] - sourceCoil[:,2] ##nodesPos was nodes!
    
    magRprimeCoils = np.sqrt(rPrimeCoils[:,:,0]**2.0 + rPrimeCoils[:,:,1]**2.0 + rPrimeCoils[:,:,2]**2.0) 
    
    unitRprimeCoils = np.zeros(np.shape(rPrimeCoils))
    for i in np.arange(3):
        unitRprimeCoils[:,:,i] = rPrimeCoils[:,:,i]/magRprimeCoils    
    
    
    plotting = 0
    if plotting == 1:
        # Plot nodes in the correct positions
        fig = pl.figure()
        from mpl_toolkits.mplot3d import Axes3D
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(loops[:,0], loops[:,1], loops[:,2], c=np.arange(len(loops)))
        ax.scatter(nodesPos[:,0], nodesPos[:,1], nodesPos[:,2], c=np.arange(len(nodesPos)))
        ax.elev = 0
        ax.azim = 270
        pl.show()
        
        # Plot nodes in the correct positions - with all the slices
        fig = pl.figure()
        from mpl_toolkits.mplot3d import Axes3D
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(sourceCoil[:,0], sourceCoil[:,1], sourceCoil[:,2], c=sourceCoil[:,3])
        ax.scatter(sourceSamplePos[:,0], sourceSamplePos[:,1], sourceSamplePos[:,2], c=sourceSamplePos[:,3])
        ax.elev = 90
        ax.azim = 0
        pl.show()
    
    return nodesPos, sourceSamplePos, rPrimeCoils, magRprimeCoils, unitRprimeCoils, sourceCoil, vertsPos

    
class testDiscSample(unittest.TestCase):
    def inputs(self):
        layers1   = 4 #6 
        sections1 = 4 #7
        slices1   = 6 #10
        
        import defineDesign as dd
        import discSample as ds
    
        # Instantiate coil, sample and fluid
        x_i = np.array([0.01015, 0.01325, 0.01015, 0.01325, 0.01015, 0.01325, 0.01015, 0.01325]) # m (loop radii)
        z_i = np.array([0.0, 0.0, 0.0031, 0.0031, 0.0131, 0.0131, 0.0162, 0.0162]) # m (relative loop heights)
        k_i = np.array([1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]) # binary (current direction)
        
        sz = dd.coil( x_i, z_i, k_i, 
                  178E3, # Hz (frequency of the coil current)
                  250.,  # A (current in coil)
                  0.003  # m (coil tube radius)
                  )        
        
        nodes1, nodesSph1, sourceSample, sourceSampleSph, rPrimeSample, magRprimeSample, unitRprimeSample, verts, zeniths, azimuths = ds.discSample(0.006, layers1, sections1, slices1)
            
        nodesPos1, sourceSamplePos1, rPrimeCoils, magRprimeCoils, unitRprimeCoils, sourceCoil, vertsPos = coilsRPrime(0.01, sz.loops, sz.loopsSph, nodes1, sourceSample, layers1, sections1, slices1, verts)
        loops1 = sz.loops
         
        return loops1, sourceCoil, nodesPos1, sourceSamplePos1
        
    def testGeo(self):
        # TEST - loops: plot the geometry
        loops1, sourceCoil, nodesPos1, sourceSamplePos1 = self.inputs()
        fig = pl.figure()
        from mpl_toolkits.mplot3d import Axes3D
        ax = fig.add_subplot(111, projection='3d')
        #ax.scatter(loops1[:,0], loops1[:,1], loops1[:,2], c=np.arange(len(loops1)))
        ax.scatter(nodesPos1[:,0], nodesPos1[:,1], nodesPos1[:,2], c=np.arange(len(nodesPos1)))
        pl.axis('equal')
        pl.axis([-0.015, 0.015, -0.015, 0.015])
        pl.show()
    def testSourceSample(self):
        # TEST - sample source coordinates: plot the sample source coordinates (with area indicated by the colour)
        loops1, sourceCoil, nodesPos1, sourceSamplePos1 = self.inputs()
        fig = pl.figure()
        from mpl_toolkits.mplot3d import Axes3D
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(sourceCoil[:,0], sourceCoil[:,1], sourceCoil[:,2], c=sourceCoil[:,3])
        ax.scatter(sourceSamplePos1[:,0], sourceSamplePos1[:,1], sourceSamplePos1[:,2], c=sourceSamplePos1[:,3])
        pl.show()

        
# Run unittests when module is called as a script
if __name__ == '__main__':
    
    try:
        unittest.main()
    except SystemExit as inst:
        if inst.args[0] is True: # raised by sys.exit(True) when tests failed
            raise