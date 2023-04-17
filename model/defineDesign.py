# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 13:27:33 2015

@author: Suzanne
"""
import numpy as np

# Classes (structures) to define the input: coil, sample and fluid
class coil(object):
    # coil geometry class
    def __init__(self, x_i, z_i, k_i, freq, current, tubeRadius):
        self.f     = freq
        self.I     = current
        self.R     = tubeRadius
        self.x_i   = x_i
        self.z_i   = z_i
        self.k_i   = k_i
        
        self.omega = freq*(2.0*np.pi)
        
        l = np.zeros((len(x_i),5))
        l[:,0] = x_i
        l[:,2] = z_i
        l[:,3] = (np.pi*self.R**2.0) # loop cross-sectional area
        l[:,4] = k_i
        self.loops = l
        
        self.loopsSph        = np.zeros(np.shape(l))
        self.loopsSph[:,3:5] = l[:,3:5]
        self.loopsSph[:,0]   = np.sqrt(l[:,0]**2.0 + l[:,2]**2.0)
        self.loopsSph[:,1]   = np.arctan2(l[:,0],l[:,2])
        
        self.Ivec = (self.loops[:,4]*self.I) + (np.abs(self.loops[:,4]-1)*(-self.I))
                        
        
class sample(object):
    # sample class
    def __init__(self, emissivity, electricalConductivity, magneticPermeability, radius, density, liquidusTemp):
        self.epsilon = emissivity
        self.sigma   = electricalConductivity
        self.mu0     = magneticPermeability
        self.R       = radius
        self.rho     = density
        self.Tm      = liquidusTemp
        self.posVar  = 0.0

        
class fluid(object):
    def __init__(self, fluidTemp, envTemp, thermalCond, velocity, density, kinematicViscosity, specificHeat):
        self.Tf  = fluidTemp
        self.T0  = envTemp
        self.k   = thermalCond
        self.v   = velocity
        self.rho = density
        self.eta = kinematicViscosity
        self.Cp  = specificHeat
        
        
# Instantiate coil, sample and fluid
x_i = np.array([0.01015, 0.01325, 0.01015, 0.01325, 0.01015, 0.01325, 0.01015, 0.01325]) # m (loop radii)
#z_i = np.array([0.0, 0.0, 0.0031, 0.0031, 0.0131, 0.0131, 0.0162, 0.0162])#-0.0081 # m (relative loop heights)
z_i = np.array([-0.0081, -0.0081, -0.005, -0.005, 0.005, 0.005, 0.0081, 0.0081]) # m (relative loop heights)
k_i = np.array([1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]) # binary (current direction)

ring = coil( np.array([10.0]), np.array([0.6]), np.array([1.0]),
          27E3,     # Hz (frequency of the coil current)
          1.,       # A (current in coil)
          0.226     # m (coil tube radius)
          )
sz = coil( x_i, z_i, k_i, 
          450E3,    # Hz (frequency of the coil current)
          400.,#200.,     # A (current in coil)
          0.0031/2  # m (coil tube radius)
          )
sz2 = coil( np.array([0.01015, 0.01325, 0.01015, 0.01325, 0.01015, 0.01325, 0.01015, 0.01325, 0.01015, 0.01325, 0.01015, 0.01325, 0.01015, 0.01325, 0.01015, 0.01325]),
	   np.array([-0.0143, -0.0143, -0.0112, -0.0112, -0.0081, -0.0081, -0.005, -0.005, 0.005, 0.005, 0.0081, 0.0081, 0.0112, 0.0112, 0.0143, 0.0143]),
	   np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), 
          450E3,    # Hz (frequency of the coil current)
          400.,#200.,     # A (current in coil)
          0.0031/2  # m (coil tube radius)
          )
fj = coil( np.array([0.021]), np.array([0.0]), np.array([1.0]),
          270E3,    # Hz (frequency of the coil current)
          500.,     # A (current in coil)
          0.003     # m (coil tube radius)
          )
mo = coil( np.array([0.0116, 0.0162, 0.0208, 0.0116, 0.0162, 0.0208, 0.0254, 0.030, 0.0116, 0.0162, 0.0208]), 
          #np.array([0.0, 0.0, 0.0, 0.0046, 0.0046, 0.0046, 0.0046, 0.0046, 0.0146, 0.0146, 0.0146]), 
          #np.array([0.0, 0.0, 0.0, 0.005, 0.005, 0.005, 0.005, 0.005, 0.015, 0.015, 0.015]),
          np.array([0.0, 0.0, 0.0, 0.0046, 0.0046, 0.0046, 0.0046, 0.0046, 0.0146, 0.0146, 0.0146]) - 0.0023,
          #np.array([0.0, 0.0, 0.0, 0.005, 0.005, 0.005, 0.005, 0.005, 0.015, 0.015, 0.015]) - 0.0023,
          np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]),
          450E3,    # Hz (frequency of the coil current)
          np.sqrt(0.87*15000/0.195),     # = 258.69494955077289 # A (current in coil)
          0.0046     # m (coil tube radius)
          )          
ke = coil( np.array([0.020, 0.027, 0.034, 0.041, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020]), 
          np.array([0.0, 0.0, 0.0, 0.0, 0.007, 0.014, 0.021, 0.028, 0.035, 0.055, 0.062]) - 0.035, 
          np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0]),
          450E3,    # Hz (frequency of the coil current)
          500.,     # A (current in coil)
          0.003     # m (coil tube radius)
          )
ro = coil( np.array([0.012, 0.012, 0.012, 0.016, 0.012, 0.016]), 
          np.array([0.004, -0.004, 0.0, 0.004, 0.012, 0.016]), 
          np.array([1.0, 1.0, 1.0, 1.0, 0.0, 0.0]),
          178E3,    # Hz (frequency of the coil current)
          170.,     # A (current in coil)
          0.0032    # m (coil tube radius)
          )
roOpt = coil( np.array([0.012, 0.016, 0.012, 0.016, 0.020, 0.012, 0.016]), 
          np.array([0.0, 0.0, 0.004, 0.004, 0.004, 0.012, 0.012]), 
          np.array([1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0]),
          178E3,    # Hz (frequency of the coil current)
          170.,     # A (current in coil)
          0.0032    # m (coil tube radius)
          )
lb = coil( np.array([0.0081, 0.0083, 0.0085, 0.0088, 0.0091, 0.0109, 0.0141, 0.0164, 0.0172, 0.0134, 0.0135, 0.0135]), 
          np.array([0.0006, 0.0089, 0.0048, 0.0124, 0.0163, 0.0202, 0.0243, 0.0289, 0.0334, 0.0398, 0.0492, 0.0443]), 
          np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]),
          179E3,    # Hz (frequency of the coil current)
          350.,     # A (current in coil)
          0.003     # m (coil tube radius)
          )          
fsangle = 10*np.pi/180 # radians
fs = coil( np.array([0.008, 0.012, 0.008+(0.004*np.tan(fsangle)), 0.008+(0.008*np.tan(fsangle)), 0.010, 0.010]), 
          np.array([0.0, 0.0, 0.004, 0.008, 0.020, 0.024]), 
          np.array([1.0, 1.0, 1.0, 1.0, 0.0, 0.0]),
          200E3,    # Hz (frequency of the coil current)
          210.,     # A (current in coil)
          0.003     # m (coil tube radius)
          )          
loop1dir1 = coil( np.array([0.01]), np.array([0.0]), np.array([0.0]), 
          450E3,    # Hz (frequency of the coil current)
          900.,     # A (current in coil)
          0.0031    # m (coil tube radius)
          )
loop1dir2 = coil( np.array([0.01]), np.array([0.0]), np.array([1.0]), 
          450E3,    # Hz (frequency of the coil current)
          900.,     # A (current in coil)
          0.0031    # m (coil tube radius)
          )
loop2dir1 = coil( np.array([0.01, 0.01]), np.array([-0.010, 0.010]), np.array([0.0, 1.0]), 
          450E3,    # Hz (frequency of the coil current)
          900.,     # A (current in coil)
          0.0031    # m (coil tube radius)
          )
loop2dir2 = coil( np.array([0.01, 0.01]), np.array([-0.010, 0.010]), np.array([0.0, 1.0]), 
          450E3,    # Hz (frequency of the coil current)
          900.,     # A (current in coil)
          0.0031    # m (coil tube radius)
          )
Al = sample(0.1, #0.0934,#        # 0.1 # (emissivity of liquid droplet, Royer et al, 2013)
            4.573E6, #4252890.,#      # 4.573E6 # S/m (electrical conductivity, Royer et al, 2013)
            4E-7*np.pi,   # H/m (magnetic permeability (free space OR sample? air ~ Al ~ Cu), Royer et al, 2013)
            0.005,        # m (sample radius)
            2702.,        # kg/m^3 (density of Al @ 300K, Cengel)
            660.+273.     # K (liquidus temperature of Al, engineering toolbox)
            )            
Cu = sample(0.11,       # (emissivity of liquid droplet, Cengel - copper commercial sheet)
            5.96E7,     # S/m (electrical conductivity, Wikipedia - Electrical conductivity)
            1.256637061435917e-006, #4E-7*np.pi, # H/m (magnetic permeability (free space OR sample? air ~ Al ~ Cu), internet)
            0.005,      # m (sample radius)
            8933.,      # kg/m^3 (density of Cu @ 300K, Cengel - pure copper)
            1083.+273.  # K (liquidus temperature of Cu, internet: http://www.azom.com/article.aspx?ArticleID=2856#_Melting_Point_of)
            )     
CuRing = sample(0.11,   # (emissivity of liquid droplet, Cengel - copper commercial sheet)
            5.96E7,     # S/m (electrical conductivity, Wikipedia - Electrical conductivity)
            1.256637061435917e-006, #4E-7*np.pi, # H/m (magnetic permeability (free space OR sample? air ~ Al ~ Cu), internet)
            10.0,       # m (sample radius)
            8933.,      # kg/m^3 (density of Cu @ 300K, Cengel - pure copper)
            1083.+273.  # K (liquidus temperature of Cu, internet: http://www.azom.com/article.aspx?ArticleID=2856#_Melting_Point_of)
            )
Fe = sample(0.7,       # (emissivity of liquid droplet, Cengel: oxidized iron @ 500-900K, 0.64-0.78)
            1E7,#/22,  # S/m (electrical conductivity, Wikipedia @ 20 deg C)
            4E-7*np.pi,    # H/m (magnetic permeability (free space OR sample? air ~ Al ~ Cu), Wikipedia)
            0.0031,    # m (sample radius)
            7870.,     # kg/m^3 (density of Fe @ 300K, Cengel)
            1810.      # K (liquidus temperature of Fe, Cengel)
            )
Ni = sample(0.47,     # (emissivity of liquid droplet, Cengel: oxidized nickel @ 450-1000K, 0.37-0.57)
            14.3E6,   # S/m (electrical conductivity, Wikipedia @ 20 deg C)
            100*4E-7*np.pi, #1.26E-4,   # 1.26E-4 - 7.54E-4 # H/m (magnetic permeability (free space OR sample? air ~ Al ~ Cu), Wikipedia)
            (3.0*0.001/(4.0*np.pi*8900.))**(1.0/3.0),  # m (sample radius)
            8900.,    # kg/m^3 (density of Ni @ 300K, Cengel)
            1728.     # K (liquidus temperature of Pure Ni, Cengel)
            )
Zn = sample(0.25,     # (emissivity of liquid droplet, Cengel: oxidized zinc @ 300K)
            2.75E6,   # S/m (electrical conductivity, Kermanpur et al. @ 900 K)
            200*4E-7*np.pi, #1.26E-4,   # 1.26E-4 - 7.54E-4 # H/m (magnetic permeability (free space OR sample?), Kermanpur et al.)
            (3.0*0.003/(4.0*np.pi*6290.))**(1.0/3.0),  # m (sample radius)
            6290.,    # kg/m^3 (density of Zn @ 1000K, Kermanpur et al.)
            692.     # K (liquidus temperature of Zn, Kermanpur et al.)
            )            
Ar = fluid(298.,      # K (fluid temperature, Royer et al, 2013)
           298.,      # K (temperature of surroundings, Royer et al, 2013)
           152E-3, # #0.028687,#  # 152E-3, # W/(mK) (thermal conductivity of fluid, Royer et al, 2013)
           0.0184,    # m/s (fluid velocity, Royer et al, 2013) 
           0.1625,  # 0.1797, #0.1625, # kg/m^3 (fluid density, Royer et al, 2013)
           199E-7, #  #1.9089E-5, # 199E-7, # m^2/s (kinematic viscosity, Royer et al, 2013)
           520.3      # J/kgK (specific heat of argon @ 298K, Cengel)
           )
