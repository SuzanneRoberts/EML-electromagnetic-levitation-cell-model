def mutualInduc(R1, z1, R2, z2, n):
    # REQUIRED MODULES: math, numpy as np, scipy.special as sp
    # Function to compute the mutual inductance [M]
    # when provided with the radii and axial positions
    # of two concentric loops
    # using the Neumann formula
    
    # Loop 2 has to be the loop with the smaller radius
    if R1 < R2:
        tempR = R1
        R1    = R2
        R2    = tempR
        tempz = z1
        z1    = z2
        z2    = tempz        

    # Constant
    mu0 = 4.0E-7*math.pi # H/m or N/A^2 (magnetic permeability of free space)

    # co-planar test case
    testCase = mu0 * math.pi * R2**2.0 / (2.0*R1) # Nm/A^2 (MIT notes, eq. 11.1.11)
    #print('test case',testCase)

    # ring-ring formula (according to Fromm & Jehn, reduces to co-planar test case)
    mutualInducFor = 0.5*mu0*math.pi*R2**2.0 * R1**2.0/(R1**2.0 + (z2-z1)**2.0)**(3.0/2.0) # Nm/A^2 (Fromm & Jehn Eqs. 5&7)

    # Maxwell's formulae in elliptical integrals
    d  = abs(z2-z1)
    r1 = math.sqrt((R1+R2)**2.0 + d**2.0)
    r2 = math.sqrt((R1-R2)**2.0 + d**2.0)
    k1 = (r1-r2)/(r1+r2) # modulus of the complete elliptic integrals F1 and E1
    m1 = k1**2.0 # parameter of the elliptic integral
    F1 = sp.ellipk(m1)
    E1 = sp.ellipe(m1)
    mutualInducEllip1 = 2.0*mu0*math.sqrt(R1*R2)*(F1-E1)/math.sqrt(k1)

    k = (2.0*math.sqrt(R1*R2))/math.sqrt((R1+R2)**2.0+d**2.0) # modulus of the complete elliptic integrals F1 and E1
    m = k**2.0 # parameter of the elliptic integral
    F = sp.ellipk(m)
    E = sp.ellipe(m)
    mutualInducEllip = mu0*math.sqrt(R1*R2)*((((2.0/k)-k)*F)-((2.0/k)*E))

    print(mutualInducEllip1,mutualInducEllip)
    
    # numerical integration
    dtheta   = 2.0*math.pi/n
    ds1      = dtheta*R1 # from: theta_rad = s/r
    # check
    #print(R1,n*dl1/(2*math.pi))
    ds2      = dtheta*R2
    # check
    #print(R2,n*dl2/(2*math.pi))
    integral = 0.0 # initialization
    for theta1 in dtheta*np.arange(n):
        for theta2 in dtheta*np.arange(n):
            # use the distance formula to find rPrime
            x1 = R1*math.cos(theta1)
            x2 = R2*math.cos(theta2)
            y1 = R1*math.sin(theta1)
            y2 = R2*math.sin(theta2)
            rPrime = math.sqrt((x1-x2)**2.0 + (y2-y1)**2.0 + (z1-z2)**2.0)
            func   = 1/rPrime
            integral += ds1*ds2*func
            
    mutualInducInt = integral * mu0 / (4.0*math.pi) # Eq.12 El-Kaddah & Szekely

    absErr = abs(mutualInducEllip - mutualInducFor)
              
    return mutualInducEllip, absErr, mutualInducInt, mutualInducFor

#==============#
# Test program #
#==============#
import math
import numpy as np
import pylab as pl
import scipy.special as sp
print('integration',mutualInduc(0.25, 0.01, 0.5, 0.01, 5))

# 1.  Test integration
# 1.1 Convergence
M_convVec  = np.zeros((10,1))
x = 2**np.arange(10.0)
x = np.reshape(x,(10,1))
for a in np.arange(10.0):
    (M_ell, absErr, M_int, M_for) = mutualInduc(0.02, 0.1, 0.04, 0.1, x[((a,0))])
    M_convVec[((a,0))] = M_ell

pl.figure(1)
pl.semilogy(x, M_convVec)
pl.xlabel('n')
pl.ylabel('M')

# 1.2 Check that sum(dl) = circumference of circle (coded)
# 1.3 Check that the results stays the same when theta1 and theta2 are swapped
# print(mutualInduc(0.02, 0.1, 0.4, 0.1, 40))
# 1.4 Compare to a Matlab result / Matlab symbolic toolbox result
# 1.5 Continue with the rest of the model and compare with experimental data

# 2. Compare numerical integration with the formulas
#    and test the effect of increasing the ratio of radii
M_intVec  = np.zeros((20,1))
M_forVec  = np.zeros((20,1))
smallR    = 0.02
largeRvec = smallR*(1.1*np.arange(20.0) + 2.0)
largeRvec = np.reshape(largeRvec,(20,1))
for a in np.arange(20.0):
    largeR    = largeRvec[(a)]
    (M_ell, absErr, M_int, M_for) = mutualInduc(smallR, 0.1, largeR, 0.1, 40)
    M_intVec[((a,0))] = M_ell
    M_forVec[((a,0))] = M_for

pl.figure(2)
pl.plot(largeRvec/smallR,M_intVec)
pl.plot(largeRvec/smallR,M_forVec) #could *25 to compare shape
pl.xlabel('Ratio of radii, R1/R2')
pl.ylabel('Mutual inductance, M [H]')
pl.legend(('Elliptical integrals','Fromm & Jehn formula'),loc='right')

pl.figure(3)
absErr = abs(M_intVec - M_forVec)
pl.semilogy(absErr)
pl.xlabel('Ratio of radii, R1/R2')
pl.ylabel('Absolute error')

pl.show() # display all figures up to here
