import numpy as np, pylab as pl, scipy.special as sp, unittest, copy, time


# Set up matplotlib setting to format figures by setting mpl.rcParams, 
# see http://matplotlib.org/users/customizing.html
# and http://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
pl.rcParams['lines.linewidth'] = 2
pl.rcParams['lines.markersize'] = 10
pl.rcParams['legend.fontsize'] = 16
pl.rcParams['axes.labelsize'] = 18
pl.rcParams['axes.titlesize'] = 18
pl.rcParams['xtick.labelsize'] = 16
pl.rcParams['ytick.labelsize'] = 16



def R1greaterR2(R1, z1, R2, z2):
    # Check that loop 1 is the loop with the larger radius
    if R1 < R2:
        tempR = R1
        R1    = R2
        R2    = tempR
        tempz = z1
        z1    = z2
        z2    = tempz
    return R1, z1, R2, z2


def coPlanarTestCase(R1, z1, R2, z2, mu0):
    # REQUIRED MODULES: numpy as np
    # Function to compute the mutual inductance
    # when provided with the radii and axial positions
    # of two concentric, coplanar loops.
    # Uses an analytical expression derived for this case as available in
    # MIT Inductance and Magnetic energy notes, eq. 11.1.11.
    # Assumptions:  1. concentric loops
    #               2. coplanar loops
    #               3. R1 >> R2
    
    # Loop 2 has to be the loop with the smaller radius
    (R1, z1, R2, z2) = R1greaterR2(R1, z1, R2, z2)       

    # co-planar test case
    return mu0 * np.pi * R2**2.0 / (2.0*R1) # Nm/A^2 (MIT notes, eq. 11.1.11)


def MIformula(R1, z1, R2, z2, mu0):
    # REQUIRED MODULES: numpy as np
    # Function to compute the mutual inductance of two coaxial loops
    # when provided with the radii and axial positions
    # using the formula used by Fromm & Jehn reduces to co-planar test case
    
    # Loop 2 has to be the loop with the smaller radius
    if (np.size(R1) != np.size(R2)): # & ((np.size(R1)|np.size(R2)) == 1):
        if np.size(R1) == 1:
            R1value = R1
            R1 = R1value*np.ones(np.size(R2))
            z1value = z1
            z1 = z1value*np.ones(np.size(z2))
        if np.size(R2) == 1:
            R2value = R2
            R2 = R2value*np.ones(np.size(R1))
            z2value = z2
            z2 = z2value*np.ones(np.size(z1))
    
    R1new = np.zeros(np.size(R1))
    R2new = np.zeros(np.size(R2))
    z1new = np.zeros(np.size(z1))
    z2new = np.zeros(np.size(z2))
    
    if np.size(R1) > 1:
        for n in np.arange(len(R2)):
            (R1new[n], z1new[n], R2new[n], z2new[n]) = R1greaterR2(R1[n], z1[n], R2[n], z2[n])

    else:
        (R1new, z1new, R2new, z2new) = R1greaterR2(R1, z1, R2, z2)
    
    # ring-ring formula (according to Fromm & Jehn, reduces to co-planar test case)
    return 0.5*mu0*np.pi*R2new**2.0 * R1new**2.0/(R1new**2.0 + (z2new-z1new)**2.0)**(3.0/2.0) # Nm/A^2 (Fromm & Jehn Eqs. 5&7)


def MIint(R1, z1, R2, z2, mu0, n):
    # REQUIRED MODULES: numpy as np
    # Function to compute the mutual inductance [M]
    # when provided with the radii and axial positions
    # of two concentric loops
    # using numerical integration of the Neumann formula
    

    # numerical integration
    dtheta   = 2*np.pi/n
    integral = 0.0 # initialization
    x1prev = R1*np.cos(-dtheta)
    x2prev = R2*np.cos(-dtheta)
    y1prev = R1*np.sin(-dtheta)
    y2prev = R1*np.sin(-dtheta)
    for theta1 in dtheta*np.arange(n):
        for theta2 in dtheta*np.arange(n):
            # use the distance formula to find rPrime
            x1 = R1*np.cos(theta1)
            x2 = R2*np.cos(theta2)
            y1 = R1*np.sin(theta1)
            y2 = R2*np.sin(theta2)
            dl1 = np.array([[x1-x1prev],[y1-y1prev]])
            dl2 = np.array([[x2-x2prev],[y2-y2prev]])
            ldotl = np.sum(dl1*dl2);
            #print(ldotl)
            rPrime = np.sqrt((x1-x2)**2.0 + (y1-y2)**2.0 + (z1-z2)**2.0)
            func   = ldotl/rPrime 
            integral += func
            x2prev = copy.copy(x2)
            y2prev = copy.copy(y2)
                
        x1prev = copy.copy(x1)
        y1prev = copy.copy(y1)
            
    mutualInducInt = integral * mu0 / (4.0*np.pi) # Eq.12 El-Kaddah & Szekely
          
    return mutualInducInt


def MIellip(R1, z1, R2, z2, mu0):
    # REQUIRED MODULES: numpy as np, scipy.special as sp
    # Function to compute the mutual inductance of two coaxial loops
    # when provided with the radii and axial positions
    # using Maxwell's formulae in elliptical integrals for the Neumann formula
    # Reference: Ross & Cohen, Tarapore & Evans
    
    # Maxwell's formulae in elliptical integrals
    d  = abs(z2-z1)
    k = (2.0*np.sqrt(R1*R2))/np.sqrt((R1+R2)**2.0+d**2.0) # modulus of the complete elliptic integrals F1 and E1
    m = k**2.0 # parameter of the elliptic integral
    F = sp.ellipk(m)
    E = sp.ellipe(m)
    mutualInducEllip = mu0*np.sqrt(R1*R2)*((((2.0/k)-k)*F)-((2.0/k)*E)) # Ross & Cohen Eq.1 * mu0/(4*pi)

    # check that Ross & Cohen Eq.1 and Eq.2 gives the same answer
    #r1 = np.sqrt((R1+R2)**2.0 + d**2.0)
    #r2 = np.sqrt((R1-R2)**2.0 + d**2.0)
    #k1 = (r1-r2)/(r1+r2) # modulus of the complete elliptic integrals F1 and E1 (see Wikipedia)
    #m1 = k1**2.0 # parameter of the elliptic integral (see Wikipedia)
    #F1 = sp.ellipk(m1)
    #E1 = sp.ellipe(m1)
    #mutualInducEllip1 = 2.0*mu0()*np.sqrt(R1*R2)*(F1-E1)/np.sqrt(k1) # Ross & Cohen Eq.2 * mu0/(4*pi)
    #print(mutualInducEllip1,mutualInducEllip)
    
    return mutualInducEllip
    
    
def SIellip(R, a, mu0):
    # REQUIRED MODULES: numpy as np, scipy.special as sp
    # Function to compute the self-inductance of circular loop
    # when provided with the loop radius and wire radius
    # using elliptical integrals
    # Reference: Ross & Cohen, Tarapore & Evans

    # R   - loop radius
    # a   - wire radius
    
    k = np.sqrt(4.0*R*(R - a)/((2.0*R - a)**2.0)) # modulus of the complete elliptic integrals F1 and E1
    m = k**2.0 # parameter of the elliptic integral
    F = sp.ellipk(m)
    E = sp.ellipe(m)
    selfInducEllip = mu0*(2*R - a)*((1 - (k**2.0/2.0))*F - E) # Ross & Cohen Eq.1 * mu0/(4*pi)
    
    return selfInducEllip
    
    
def SIwien(R, rho, mu0):
    # REQUIRED MODULES: numpy as np
    # Function to compute the self-inductance of a circle
    # when provided with the radius and cross-sectional radius of the circular loop
    # using Wien formula
    # Reference: Ross & Cohen, Tarapore & Evans

    # R   - loop radius
    # rho - wire radius

    #term1 = (1 + (1.0/8.0)*(rho**2.0/R**2.0)) * np.log10(8*R/rho)
    #L = mu0()*R*(term1 - 0.0083*(rho**2.0/R**2.0) - 1.75)
    
    L = mu0*R*(np.log(8*R/rho) - 2.0 )
    
    return L



class testMutualInduc(unittest.TestCase):
    # test that Fromm & Jehn's formula reduces to coPlanarTestCase for the planar case
    # now only use Fromm & Jehn's formula as a test
    def testFrommJehnPlanar(self):
        self.assertEqual(coPlanarTestCase(0.05,0.01,0.02,0.01,4.0E-7*np.pi), MIformula(0.05,0.01,0.02,0.01,4.0E-7*np.pi))

#    # test numerical integration
#    def testNumInt(self):
#        # compare numerical integral with Fromm and Jehn formula
#        self.assertTrue(abs(MIint(0.05,0.01,0.02,0.01,4.0E-7*np.pi,10) - MIformula(0.05,0.01,0.02,0.01,4.0E-7*np.pi)) < 1E-8)
#                
#        # Convergence
#        M_convVec  = np.zeros((10,1))
#        x = 50.0*(np.arange(10.0)+1)
#        x = np.reshape(x,(10,1))
#        for a in np.arange(10):
#            M_ell = MIint(0.02, 0.1, 0.04, 0.1, 4.0E-7*np.pi, x[((a,0))])
#            M_convVec[((a,0))] = M_ell
#
#        pl.figure()
#        pl.plot(x, M_convVec)
#        pl.plot(x, MIellip(0.02, 0.1, 0.04, 0.1, 4.0E-7*np.pi)*np.ones(len(x)), 'r')
#        pl.xlabel('Number of integration intervals, n')
#        pl.ylabel('Mutual inductance from numerical integration')
#                
#        pl.figure()
#        #absErr = abs(M_convVec - M_convVec[9])
#        absErr = abs(M_convVec - MIellip(0.02, 0.1, 0.04, 0.1, 4.0E-7*np.pi))
#        pl.semilogy(x, absErr,'o')
#        pl.xlabel('Number of integration intervals, n')
#        pl.ylabel('Absolute error (Numerical integration)')
#        pl.show() # display all figures up to here
#
#        # Try: Check that the results stays the same when theta1 and theta2 are swapped
        

    # test elliptical integrals
    def testMIellip(self):
        # compare numerical integral with Fromm and Jehn formula
        self.assertTrue(abs(MIellip(0.05,0.03,0.02,0.01, 4.0E-7*np.pi) - MIformula(0.05,0.03,0.02,0.01, 4.0E-7*np.pi)) < 1E-9)
        
        # plot values, abs error and % error for different ratios of radii
        M_intVec  = np.zeros((20,1))
        M_forVec  = np.zeros((20,1))
        smallR    = 0.02 # radius of the smaller loop
        largeRvec = smallR*(1.1*np.arange(20.0) + 2.0)
        largeRvec = np.reshape(largeRvec,(20,1))
        for a in np.arange(20):
            largeR    = largeRvec[(a)]
            M_ell = MIellip(smallR, 0.3, largeR, 0.1, 4.0E-7*np.pi)
            M_for = MIformula(smallR, 0.3, largeR, 0.1, 4.0E-7*np.pi)
            M_intVec[((a,0))] = M_ell
            M_forVec[((a,0))] = M_for

        pl.figure()
        pl.plot(largeRvec/smallR,M_intVec,'o')
        pl.plot(largeRvec/smallR,M_forVec,'s')
        pl.xlabel('Ratio of radii, R1/R2')
        pl.ylabel('Mutual inductance, M [H]')
        pl.legend(('Elliptical integrals','Fromm & Jehn formula'),loc='right')

        pl.figure()
        absErr = abs(M_intVec - M_forVec)
        pl.semilogy(absErr,'o')
        pl.xlabel('Ratio of radii, R1/R2')
        pl.ylabel('Absolute error (Elliptical integrals)')

        pl.figure()
        pcentErr = abs(M_intVec - M_forVec)/M_intVec*100
        pl.plot(pcentErr,'o')
        pl.xlabel('Ratio of radii, R1/R2')
        pl.ylabel('Percentage error (Elliptical integrals)')
        
        pl.show() # display all figures up to here

    def testEllipVecArgs(self):
        # test elliptical integral function with vectors as arguments
        scalar1 = MIellip(0.01, 0.01, 0.03, 0.04, 4.0E-7*np.pi)
        scalar2 = MIellip(0.02, 0.01, 0.03, 0.04, 4.0E-7*np.pi)
        scalar3 = MIellip(0.03, 0.01, 0.03, 0.04, 4.0E-7*np.pi)
        vector  = np.array([0.01,0.02,0.03])
        oneVector = MIellip(vector, 0.01*np.ones(3), 0.03, 0.04, 4.0E-7*np.pi)
        twoVector = MIellip(vector, 0.01*np.ones(3), 0.03*np.ones(3), 0.04*np.ones(3), 4.0E-7*np.pi)
        self.assertEqual(scalar1, oneVector[0])
        self.assertEqual(scalar2, oneVector[1])
        self.assertEqual(scalar3, oneVector[2])
        self.assertEqual(scalar1, twoVector[0])
        self.assertEqual(scalar2, twoVector[1])
        self.assertEqual(scalar3, twoVector[2])
        
    def testFormulaVecArgs(self):
        # test elliptical integral function with vectors as arguments
        scalar1 = MIformula(0.01, 0.01, 0.03, 0.04, 4.0E-7*np.pi)
        scalar2 = MIformula(0.02, 0.01, 0.03, 0.04, 4.0E-7*np.pi)
        scalar3 = MIformula(0.03, 0.01, 0.03, 0.04, 4.0E-7*np.pi)
        vector  = np.array([0.01,0.02,0.03])
        oneVector = MIformula(vector, 0.01*np.ones(3), 0.03, 0.04, 4.0E-7*np.pi)
        print(vector)
        twoVector = MIformula(vector, 0.01*np.ones(3), 0.03*np.ones(3), 0.04*np.ones(3), 4.0E-7*np.pi)
        self.assertEqual(scalar1, oneVector[0])
        self.assertEqual(scalar2, oneVector[1])
        self.assertEqual(scalar3, oneVector[2])
        self.assertEqual(scalar1, twoVector[0])
        self.assertEqual(scalar2, twoVector[1])
        self.assertEqual(scalar3, twoVector[2])
        
    def testFormulaSmallSampleAssumption(self):
        small2largeFormulaVec = []
        small2largeEllipVec = []
        small2largeIntVec = []
        
        secondLoopR = np.arange(10.0) + 1   
        firstLoopR  = 10.0
        
        for Rvar in secondLoopR:
            tic1 = time.time()
            scalar1 = MIformula(firstLoopR, 0, Rvar, 0.6, 4.0E-7*np.pi)
            toc1 = time.time()
            time1 = toc1 - tic1
            print('time1 - formula',time1)
            small2largeFormulaVec.append(scalar1)
            
            tic2 = time.time()
            scalar2 = MIellip(firstLoopR, 0, Rvar, 0.6, 4.0E-7*np.pi)
            toc2 = time.time()
            time2 = toc2 - tic2
            print('time2 - ellip ',time2)
            small2largeEllipVec.append(scalar2)
            
            tic3 = time.time()
            scalar3 = MIint(firstLoopR, 0, Rvar, 0.6, 4.0E-7*np.pi, 200)
            toc3 = time.time()
            time3 = toc3 - tic3
            print('time3 - int ',time3)
            small2largeIntVec.append(scalar3)

        
        pl.figure()
        pl.plot(secondLoopR/firstLoopR, small2largeFormulaVec, '-ko', markerfacecolor = 'w', markeredgewidth = 2)
        pl.plot(secondLoopR/firstLoopR, small2largeIntVec, '-ks', markerfacecolor = 'w', markeredgewidth = 2)
        pl.plot(secondLoopR/firstLoopR, small2largeEllipVec, '-kx', markerfacecolor = 'w', markeredgewidth = 2)
        #pl.xlabel('Ratio of loop radii, R_{loop 2}/R_{loop 1}')
        pl.xlabel('Ratio of loop radii: (radius of loop 2)/(radius of loop 1)')
        pl.ylabel('Mutual Inductance')
        #pl.title('Effect of the assumption used in MIformula that the sample is small relative to the coil')
        pl.legend(['Formula','Trapezoidal rule (200 increments)','Elliptical integrals'], loc = 'upper left')
        pl.show()
        
    def testSIellipAndSIwien(self):
        small2largeFormulaVec = []
        small2largeSIellip = []
        small2largeSIwien = []
        
        loopR = 20.0
        a = (np.arange(20.0) + 1.0)
        aRratio = a/loopR
        
        for aR in aRratio:
            scalar1 = MIformula(loopR, 0, loopR, 0, 4.0E-7*np.pi)
            small2largeFormulaVec.append(scalar1)
            
            scalar2 = SIellip(loopR, aR, 4.0E-7*np.pi)
            small2largeSIellip.append(scalar2)
            
            scalar3 = SIwien(loopR, aR, 4.0E-7*np.pi)
            small2largeSIwien.append(scalar3)
                 
        
        pl.figure()
        pl.plot(aRratio, small2largeFormulaVec, '-bo')
        pl.plot(aRratio, small2largeSIellip, '-gs')
        pl.plot(aRratio, small2largeSIwien, '-r*')
        pl.xlabel('Wire radius to loop radius ratio, a/R')
        pl.ylabel('Self-inductance')
        pl.legend(['MIformula','SIellip','SIwien'], loc = 'upper right')
        pl.show()
        


# Run unittests when module is called as a script
if __name__ == '__main__':

    try:
        unittest.main()
    except SystemExit as inst:
        if inst.args[0] is True: # raised by sys.exit(True) when tests failed
            raise
