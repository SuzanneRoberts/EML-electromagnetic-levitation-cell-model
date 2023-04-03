import math
import numpy as np
import pylab as pl
import copy

def mutualInduc(R1, z1, R2, z2, n):
    # REQUIRED MODULES: math, numpy as np
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
    print('test case',testCase)

    # ring-ring formula (according to Fromm & Jehn, reduces to co-planar test case)
    mutualInducFor = 0.5*mu0*math.pi*R2**2.0 * R1**2.0/(R1**2.0 + (z2-z1)**2.0)**(3.0/2.0) # Nm/A^2 (Fromm & Jehn Eqs. 5&7)
    
    # numerical integration
    dtheta   = 2*math.pi/n
    ds1      = dtheta*R1 # from: theta_rad = s/r
    # check
    #print(R1,n*dl1/(2*math.pi))
    ds2      = dtheta*R2
    # check
    #print(R2,n*dl2/(2*math.pi))
    integral = 0.0 # initialization
    x1prev = R1*math.cos(-dtheta)
    x2prev = R2*math.cos(-dtheta)
    y1prev = R1*math.sin(-dtheta)
    y2prev = R1*math.sin(-dtheta)
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
            print(ldotl)
            rPrime = np.sqrt((x1-x2)**2.0 + (y1-y2)**2.0 + (z1-z2)**2.0)
            func   = ldotl/rPrime 
            integral += func
            x2prev = copy.copy(x2)
            y2prev = copy.copy(y2)
        x1prev = copy.copy(x1)
        y1prev = copy.copy(y1)
            
    mutualInducInt = integral * mu0 / (4.0*np.pi) # Eq.12 El-Kaddah & Szekely
          
    return mutualInducInt, mutualInducFor


print(mutualInduc(0.01, 0.02, 0.1, 0.02, 1000))

#==============================================================================
# #==============#
# # Test program #
# #==============#
# import math
# import numpy as np
# import pylab as pl
# print('integration',mutualInduc(0.25, 0.01, 0.5, 0.01, 5))
# 
# # 1.  Test integration
# # 1.1 Convergence
# M_convVec  = np.zeros((10,1))
# x = 2**np.arange(10.0)
# x = np.reshape(x,(10,1))
# for a in np.arange(10.0):
#     (M_int, M_for) = mutualInduc(0.02, 0.1, 0.04, 0.1, x[((a,0))])
#     M_convVec[((a,0))] = M_int
# 
# pl.semilogy(x, M_convVec)
# pl.xlabel('n')
# pl.ylabel('M')
# pl.show()
# 
# # 1.2 Check that sum(dl) = circumference of circle (coded)
# # 1.3 Check that the results stays the same when theta1 and theta2 are swapped
# # print(mutualInduc(0.02, 0.1, 0.4, 0.1, 40))
# # 1.4 Compare to a Matlab result / Matlab symbolic toolbox result
# # 1.5 Continue with the rest of the model and compare with experimental data
# 
# # 2. Compare numerical integration with the formulas
# #    and test the effect of increasing the ratio of radii
# M_intVec  = np.zeros((20,1))
# M_forVec  = np.zeros((20,1))
# smallR    = 0.02
# largeRvec = smallR*(1.5*np.arange(20.0) + 2.0)
# largeRvec = np.reshape(largeRvec,(20,1))
# for a in np.arange(20.0):
#     largeR    = largeRvec[(a)]
#     (M_int, M_for) = mutualInduc(smallR, 0.1, largeR, 0.1, 40)
#     M_intVec[((a,0))] = M_int
#     M_forVec[((a,0))] = M_for
# 
# pl.plot(largeRvec/smallR,M_intVec)
# pl.plot(largeRvec/smallR,M_forVec) #could *25 to compare shape
# pl.xlabel('Ratio of radii, R1/R2')
# pl.ylabel('Mutual inductance, M [H]')
# pl.legend(('Numerical integration','Fromm & Jehn formula'),loc='right')
# pl.show()
# 
# absErr = abs(M_intVec - M_forVec)
# pl.plot(absErr)
# pl.xlabel('Ratio of radii, R1/R2')
# pl.ylabel('Absolute error')
# pl.show()
# 
#==============================================================================
