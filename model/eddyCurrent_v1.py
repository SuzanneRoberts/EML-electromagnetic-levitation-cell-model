import math, numpy as np, pylab as pl, unittest, mutualInduc as mi
# For: Hack found online to remove the perspective view from the plot
from mpl_toolkits.mplot3d import proj3d
# For plotting
import matplotlib.pyplot as plt
from matplotlib import cm


def discSample(R, layers, sections, slices):
    # REQUIRED MODULES: math, numpy as np, pylab as pl
    layerR       = R/layers
    sectionTheta = math.pi/sections

    nodes      = np.zeros((layers*sections,7))
    nodeNbrs   = np.arange(layers*sections)
    nodes[:,0] = nodeNbrs

    count = 0
    for i in np.arange(layers):
        nodes[count:count+sections,4] = i*layerR + 0.5*layerR
        nodes[count:count+sections,3] = 0.5*(math.pi*(((i+1)*layerR)**2.0-(i*layerR)**2.0))/sections
        count += sections

    for i in np.arange(sections):
        nodes[sections*np.arange(layers)+i,5] = i*sectionTheta + 0.5*sectionTheta - math.pi/2
        
    nodes[:,1] = nodes[:,4]*np.cos(nodes[:,5])
    nodes[:,2] = nodes[:,4]*np.sin(nodes[:,5])
    
    # compute separation matrix for sample
    numBlocks = layers*sections
    rPrimeSample = np.zeros((layers*sections*slices, layers*sections, 3))
    for phi in np.arange(slices):
        for block in nodeNbrs:
            rPrimeSample[phi*numBlocks:((phi*numBlocks)+numBlocks), block, 0] = nodes[block,4] - nodes[:,4]
            rPrimeSample[phi*numBlocks:((phi*numBlocks)+numBlocks), block, 1] = nodes[block,5] - nodes[:,5]
            rPrimeSample[phi*numBlocks:((phi*numBlocks)+numBlocks), block, 2] = 0.0 - phi*(2*np.pi/slices)
            
    magRprimeSample  = np.sqrt(rPrimeSample[:,:,0]**2.0 + rPrimeSample[:,:,1]**2.0 + rPrimeSample[:,:,2]**2.0) 
    
    unitRprimeSample = np.zeros(np.shape(rPrimeSample))
    for i in np.arange(3):
        unitRprimeSample[:,:,i] = rPrimeSample[:,:,i]/magRprimeSample
    
    return nodes, rPrimeSample, magRprimeSample, unitRprimeSample


class testDiscSample(unittest.TestCase):
    def testAreas(self):
        # TEST - areas: do the areas add up to the correct value?
        R = 0.005
        nodes, rPrimeSample, MAGrPRIMEsample, UNITrPRIMEsample = discSample(R, 5, 8, 2)
        areaTotal = 0.5*(math.pi*R**2.0)
        self.assertTrue(abs(areaTotal - np.sum(nodes[:,3])) < 1E-12)
    def testGeo(self):
        # TEST - coordinates: plot the geometry (only the nodes)
        nodes, rPrimeSample, MAGrPRIMEsample, UNITrPRIMEsample = discSample(0.006, 6, 7, 2)
        pl.figure()
        pl.plot(nodes[:,4]*np.cos(nodes[:,5]),nodes[:,4]*np.sin(nodes[:,5]))
        pl.plot(nodes[:,1],nodes[:,2],'*')
        pl.axis('equal')
        pl.show()

# Classes (structures) to define the input: coil, sample and fluid

class coil(object):
    # coil geometry class
    def __init__(self, coilLoops, freq, current, tubeRadius):
        self.loops = coilLoops
        self.omega = freq
        self.I     = current
        self.R     = tubeRadius
        
class sample(object):
    # sample class
    def __init__(self, emissivity, electricalConductivity, magneticPermeability, radius, density, liquidusTemp):
        self.epsilon = emissivity
        self.sigma   = electricalConductivity
        self.mu0     = magneticPermeability
        self.R       = radius
        self.rho     = density
        self.Tm      = liquidusTemp
        
class fluid(object):
    def __init__(self, fluidTemp, envTemp, thermalCond, velocity, density, kinematicViscosity, specificHeat):
        self.Tf  = fluidTemp
        self.T0  = envTemp
        self.k   = thermalCond
        self.v   = velocity
        self.rho = density
        self.eta = kinematicViscosity
        self.Cp  = specificHeat


# Function to solve for the induced current distribution in the sample
def solveCurrent(gridpoints, coil, sample):
    # A - full quater of coefficient matrix
    numNodes = len(gridpoints)
    A = np.zeros((numNodes,numNodes))
    for i in range(numNodes):
        # adding mutual inductance terms to A
        A[:,i] = mi.MIellip(x,y,x[i],y[i])*gridpoints[i,3]
        # adding self-inductance terms to A
        A[i,i] = mi.MIformula(x[i],y[i],x[i],y[i])*gridpoints[i,3]
    
    A *= coil.omega*sample.sigma
    
    # B - diagonal quater of coefficient matrix
    b = 2*np.pi*gridpoints[:,1] # El-Kaddah & Szekely use a LH system, I use RH, therefore cos instead of their sin
    B = np.diag(b)
    
    # Assemble LHS
    C = np.zeros((2*numNodes, 2*numNodes))
    C[:numNodes,:numNodes] = A
    C[numNodes:,numNodes:] = A
    C[:numNodes,numNodes:] = B
    C[numNodes:,:numNodes] = -B
    pl.figure()
    pl.spy(C)
    
    
    # Assemble RHS
    numLoops = len(coil.loops)
    I = np.zeros((2*numNodes, 1))
    for i in range(numLoops):
        I[:numNodes,0] += - mi.MIellip(x,y,coil.loops[i,1],coil.loops[i,2])*coil.loops[i,3]*coil.I
        I[numNodes:,0] += mi.MIellip(x,y,coil.loops[i,1],coil.loops[i,2])*np.abs(coil.loops[i,3]-1)*coil.I
        
    I *= coil.omega*sample.sigma
    
    # Solve linear system
    J = np.linalg.solve(C,I)
    
    reJ   = J[:numNodes,0]
    imJ   = J[numNodes:,0]
    Jcomp = reJ + imJ*1j

    return Jcomp
        

# Function to add the coils as source points to the separation vector matrix
def coilsRPrime(loops, nodes, layers, sections, slices):
    
    (numLoops,c1) = np.shape(loops)
    numBlocks = layers*sections
    rPrimeCoils = np.zeros((slices*numLoops,numBlocks,3))
    
    for phi in np.arange(slices):
        for block in np.arange(numBlocks):
            rPrimeCoils[(phi*numLoops):((phi*numLoops)+numLoops), block, 0] = nodes[block,4] - np.sqrt(loops[:,1]**2.0 + loops[:,2]**2.0)
            rPrimeCoils[(phi*numLoops):((phi*numLoops)+numLoops), block, 1] = nodes[block,5] - np.tan(loops[:,2]/loops[:,1])
            rPrimeCoils[(phi*numLoops):((phi*numLoops)+numLoops), block, 2] = 0.0 - phi*(2*np.pi/slices)
    
    magRprimeCoils = np.sqrt(rPrimeCoils[:,:,0]**2.0 + rPrimeCoils[:,:,1]**2.0 + rPrimeCoils[:,:,2]**2.0) 
    
    unitRprimeCoils = np.zeros(np.shape(rPrimeCoils))
    for i in np.arange(3):
        unitRprimeCoils[:,:,i] = rPrimeCoils[:,:,i]/magRprimeCoils    
    
    #rPrimeWithCoils = np.vstack((rPrimeNoCoils, rPrimeCoils))    
    
    return rPrimeCoils, magRprimeCoils, unitRprimeCoils

        
# Hack found online to remove the perspective view from the plot
def orthogonal_proj(zfront, zback):
    a = (zfront+zback)/(zfront-zback)
    b = -2*(zfront*zback)/(zfront-zback)
    return np.array([[1,0,0,0],
                        [0,1,0,0],
                        [0,0,a,b],
                        [0,0,-0.0001,zback]])
                            
                            
# Function to plot a field over the sample
def plotSampleField(xvec, yvec, fieldvec, titleString, colorbarString):
    
    # Hack found online to remove the perspective view from the plot
    proj3d.persp_transformation = orthogonal_proj
    
    # Plotting the field
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    fullX = np.concatenate((xvec,-xvec))
    fullY = np.concatenate((yvec,yvec))
    fullZ = np.concatenate((fieldvec,fieldvec))
    surf=ax.plot_trisurf(fullX*1000, fullY*1000, fullZ, cmap=cm.jet, linewidth=0.0)
    
    ax.view_init(90,-90)
    ax.set_zticks([])
    cb = fig.colorbar(surf, shrink=0.5, aspect=5)
    cb.set_label(colorbarString)
    ax.set_title(titleString)
    ax.set_xlabel('mm')
    ax.set_ylabel('mm')
    

    
#==============================================================================
# Main program
#==============================================================================

# Instantiate coil, sample and fluid
r_i = np.array([0.01015, 0.01325, 0.01015, 0.01325, 0.01015, 0.01325, 0.01015, 0.01325]) # m (loop radii)
z_i = np.array([0.0, 0.0, 0.0031, 0.0031, 0.0131, 0.0131, 0.0162, 0.0162]) # m (relative loop heights)
k_i = np.array([1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]) # binary (current direction)
l = np.zeros((len(r_i),4))
l[:,0] = np.arange(len(r_i))
l[:,1] = r_i
l[:,2] = z_i
l[:,3] = k_i

sz = coil(l, 
          178E3, # Hz (frequency of the coil current)
          250.,  # A (current in coil)
          0.003  # m (coil tube radius)
          )        
Al = sample(0.0934,       # (emissivity of liquid droplet, Royer et al, 2013)
            4252890.,     # S/m (electrical conductivity, Royer et al, 2013)
            4E-7*np.pi,   # H/m (magnetic permeability (free space OR sample? air ~ Al ~ Cu), Royer et al, 2013)
            0.005,        # m (sample radius)
            2702.,        # kg/m^3 (density of Al @ 300K, Cengel)
            660.+273.     # K (liquidus temperature of Al, engineering toolbox)
            )            
Cu = sample(0.15,       # (emissivity of liquid droplet, Cengel - copper commercial sheet)
            5.96E7,     # S/m (electrical conductivity, Wikipedia - Electrical conductivity)
            4E-7*np.pi, # H/m (magnetic permeability (free space OR sample? air ~ Al ~ Cu), internet)
            0.005,      # m (sample radius)
            8933.,      # kg/m^3 (density of Cu @ 300K, Cengel - pure copper)
            1083.+273.  # K (liquidus temperature of Cu, internet: http://www.azom.com/article.aspx?ArticleID=2856#_Melting_Point_of)
            )            
Ar = fluid(298.,      # K (fluid temperature, Royer et al, 2013)
           298.,      # K (temperature of surroundings, Royer et al, 2013)
           0.028687,  # W/(mK) (thermal conductivity of fluid, Royer et al, 2013)
           0.0184,    # m/s (fluid velocity, Royer et al, 2013) 
           0.1797,    # kg/m^3 (fluid density, Royer et al, 2013)
           1.9089E-5, # m^2/s (kinematic viscosity, Royer et al, 2013)
           520.3      # J/kgK (specific heat of argon @ 298K, Cengel)
           )

# Discretize sample
NODES, rPRIMEsample, MAGrPRIMEsample, UNITrPRIMEsample = discSample(0.005, 40, 40, 6)
x = NODES[:,1]
y = NODES[:,2]
numNodes, col = np.shape(NODES)

pl.figure()
pl.plot(x,y,'o')
pl.axis('equal')

# Discretize coils
rPRIMEcoils, MAGrPRIMEcoil, UNITrPRIMEcoil = coilsRPrime(sz.loops, NODES, 40, 40, 6)

# Solve for induced current desity
Jcomp = solveCurrent(NODES, sz, Al)

Jmag  = np.sqrt(np.real(Jcomp)**2.0 + np.imag(Jcomp)**2.0)
print('Total current induced:',sum(Jmag*NODES[:,3]))

# Power calculation
specP = 0.5*np.real(Jcomp*Jcomp.conj())/Al.sigma
P  = specP*(NODES[:,3]*2.0*np.pi*NODES[:,1])
Ptotal = np.sum(P)
print('Total power absorbed',Ptotal)


# Magnetic flux density calculation
(rs, cs, ds) = np.shape(rPRIMEsample)
(rc, cc, dc) = np.shape(rPRIMEcoils)

B = np.zeros((numNodes, 3))

for f in np.arange(numNodes):
    integral = np.zeros(3) # initialise integral for each new field point
    # integrate over the current sources in the sample
    for s in np.arange(rs):
        integral += 0.0
    # integrate over the current sources in the coil
    for s in np.arange(rc):
        integral += 0.0
    
    # save result in B component vectors
    B[f,:] = integral
    
        
B *= (mi.mu0()/(4.0*np.pi))

# Lifting force calculation

# Find sample levitation position

# Stiring force calculation

# Sample temperature calculation
tol    = 1E-8
update = 2.0*tol
T      = 1000 #K - initialize temperature
Boltz = 5.67E-8
area = 4.0*np.pi*Al.R**2.0
while abs(update) > tol:
    Qrad = Boltz*Al.epsilon*area*(T**4.0 - Ar.T0**4.0)
    Nu = 2.0 + 0.6 * (2.0*Al.R*Ar.v*Ar.rho/Ar.eta)**0.5 * (Ar.Cp*Ar.eta/Ar.k)**(1/3)
    h = Ar.k * Nu / (2.0*Al.R)
    Qconv = h*area*(T - Ar.Tf)
    
    R = Ptotal - Qrad - Qconv 
    dRdT = - (4.0*Boltz*Al.epsilon*area*T**3.0) - (h*area) 
    
    update = - R/dRdT
    T += update


# Plots
plotSampleField(x, y, Jmag*1E-6, 'Current density, J', 'A/mm$^2$')
plotSampleField(x, y, specP*1E-9, 'Specific power absorbed, P', 'W/mm$^3$')
plotSampleField(x, y, np.arctan(np.imag(Jcomp)/np.real(Jcomp)), 'Phase, theta', 'Radians')


# Print tests
print('Qconv',Qconv)
print('Qrad',Qrad)
print('Sample temperature',T)

NODES1, rPRIMEsample1, MAGrPRIMEsample1, UNITrPRIMEsample1 = discSample(0.005, 3, 3, 3)
rPrimeCoils1, MAGrPRIMEcoil1, UNITrPRIMEcoil1 = coilsRPrime(sz.loops, NODES1, 3, 3, 3)
#print('UNITrPRIMEsample1', UNITrPRIMEsample1)
#print('UNITrPRIMEcoil1', UNITrPRIMEcoil1)


# Run unittests when module is called as a script
if __name__ == '__main__':
    
    try:
        unittest.main()
    except SystemExit as inst:
        if inst.args[0] is True: # raised by sys.exit(True) when tests failed
            raise
