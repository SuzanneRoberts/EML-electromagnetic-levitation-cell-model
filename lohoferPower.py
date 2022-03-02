# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 08:12:08 2015

@author: Suzanne
"""

import scipy.special as sp, numpy as np, pylab as pl, unittest, scipy.optimize as op, scipy
import defineDesign as dd

def scalarSphBessel(n,z):
    [fnValVec, dFdzVec] = sp.sph_jn(n,z)
    return fnValVec[n], dFdzVec[n]

def scalarSphBessel2(z,n):
    [fnValVec, dFdzVec] = sp.sph_jn(n,z)
    return (fnValVec[n])**2.0
    
# anonymous functions in Python 
sphBessel = lambda n, z : scalarSphBessel(n,z)
print(sphBessel(1,3.1))

def sphBesselZeros(n,x0):
    # Newton-Raphson loop
    # To find the next root, the initial guess must be the previous root plus pi
    tol = 10E-14
    update = 2*tol
    while np.abs(update) > tol:
        [fnVal, dfdz] = scalarSphBessel(n,x0)
        update        = - fnVal/dfdz
        x0           += update
        
    return x0
    
# This does not work, because sp.sph_jn only accepts scalars as arguments and op.minimize is vectorized    
#def sphBesselZerosNew(n,xmin,xmax):
#    func = lambda x: (scalarSphBessel2(n,x))**2.0
#    print('w',func(2))
#    lb = xmin
#    ub = xmin + 2.0*np.pi
#    rootList = []
#    while ub < xmax:
#        x0   = [lb + np.pi]
#        bnds = [(lb,ub)]
#        
#        fmin = op.minimize(func, x0=x0, method='SLSQP', bounds=bnds, options={'ftol': 1e-06})
#        
#        rootList.append(fmin)
#        step = ub - lb
#        lb = ub
#        ub += step
#        
#    return rootList
#        
#rootList1 = sphBesselZerosNew(1.0,0.0,10.0)
#print('rootList1',rootList1)

def skinDepth(omega, sigma, mu0):
    return np.sqrt(2/(omega*sigma*mu0))
    
def qn(Rs, dn):
    return Rs/dn

def Hofqn(qn, l):
    tol      = 10E-3
    term     = 2*tol
    H        = 0.0 # initialization
    lastRoot = 0.0 # initialization
    Hlist    = []
    
    # plot Bessel function to check roots
    vecLen = 100
    func   = np.zeros(vecLen)
    for i in np.arange(vecLen):
        [funVali, dfdzi] = scalarSphBessel(l,i)
        func[i]          = funVali
            
    pl.figure(10*l)
    pl.plot(np.arange(vecLen),func)
    pl.title(['Spherical Bessel function, n =',l])
    # end of Bessel function plot check
    
    while np.abs(term) > tol: 
        x0         = lastRoot + np.pi
        BesselRoot = sphBesselZeros(l,x0)
        # Bessel root check plot
        pl.figure(10*l)
        pl.plot(BesselRoot,0.0,'x')
        # End of Bessel root check plot
        term       = 8.0*qn**4.0/((4.0*qn**4.0) + (BesselRoot**4.0))
        H         += term
        Hlist.append(H)
        lastRoot   = x0
        
        return H, Hlist
        
scipy.savetxt('HlistFile.txt',Hlist)
#scipy.loadtxt('HlistFile.txt')


class testLohoferPower(unittest.TestCase):
    def testSphBesselZeros(self):
        # Compare the zeros my function finds to those obtained from the
        # Keisan Online Calculator: http://keisan.casio.com/exec/system/1180573465
        self.assertTrue(sphBesselZeros(1,1572.366+np.pi) - 1,575.50808105948534854 < 10E-14)
        self.assertTrue(sphBesselZeros(50, 92.683+np.pi) - 96.39974404597633670485 < 10E-14)
        
    def testScalarSphBessel(self):
        # Plot the Spherical Bessel function with n=1
        vecLen = 100
        func   = np.zeros(vecLen)
        for i in np.arange(vecLen):
            [funVali, dfdzi] = scalarSphBessel(1,i)
            func[i]          = funVali
            
        pl.figure()
        pl.plot(np.arange(vecLen),func)
        pl.title('Spherical Bessel function, n=1')
        pl.show()
        
    def testHconvergence(self):
        sample = dd.Cu
        coil   = dd.fj
        dn1 = skinDepth(coil.omega,sample.sigma,sample.mu0)
        qn1 = qn(sample.R, dn1)
        [H0,Hlist0]  = Hofqn(qn1, 0)
#        [H1,Hlist1]  = Hofqn(qn1, 1) #!
        [H2,Hlist2]  = Hofqn(qn1, 2)
#        [H3,Hlist3]  = Hofqn(qn1, 3) #!
        [H4,Hlist4]  = Hofqn(qn1, 4)
#        [H5,Hlist5]  = Hofqn(qn1, 5) #!
#        [H6,Hlist6]  = Hofqn(qn1, 6) #!
#        [H7,Hlist7]  = Hofqn(qn1, 7) #!
#        [H8,Hlist8]  = Hofqn(qn1, 8) #!
#        [H9,Hlist9]  = Hofqn(qn1, 9) #!
#        [H10,Hlist10]  = Hofqn(qn1, 10) #!
#        [H11,Hlist11]  = Hofqn(qn1, 11) #!
#        [H12,Hlist12]  = Hofqn(qn1, 12) #!
        
        pl.figure()
        pl.plot(Hlist0)
#        pl.plot(Hlist1)
        pl.plot(Hlist2)
#        pl.plot(Hlist3)
        pl.plot(Hlist4)
#        pl.plot(Hlist5)
#        pl.plot(Hlist6)
#        pl.plot(Hlist7)
#        pl.plot(Hlist8)
#        pl.plot(Hlist9)
#        pl.plot(Hlist10)
#        pl.plot(Hlist11)
#        pl.plot(Hlist12)
        pl.title('Adding terms to H(qn)')
        pl.xlabel('Number of terms')
        pl.ylabel('Value of H')
#        pl.legend(['l = 0','l = 1','l = 2','l = 3','l = 4','l = 5','l = 6','l = 7','l = 8','l = 9','l = 10','l = 11','l = 12'],loc=4)
#        
        pl.figure()
#        pl.plot([H0,H1,H2,H3,H4,H5,H6,H7,H8,H9,H10,H11,H12],'-b*')
        pl.plot([H0,H2,H4],'-b*')
        pl.title('H-values for increasing orders of the Bessel function')
        pl.xlabel('Order (l) of the Bessel function')
        pl.ylabel('Value of H')
        
#        pl.figure()
#        pl.plot([H0,H0+H1,H0+H1+H2,H0+H1+H2+H3,H0+H1+H2+H3+H4,H0+H1+H2+H3+H4+H5,H0+H1+H2+H3+H4+H5+H6,H0+H1+H2+H3+H4+H5+H6+H7,H0+H1+H2+H3+H4+H5+H6+H7+H8,H0+H1+H2+H3+H4+H5+H6+H7+H8+H9,H0+H1+H2+H3+H4+H5+H6+H7+H8+H9+H10,H0+H1+H2+H3+H4+H5+H6+H7+H8+H9+H10+H11,H0+H1+H2+H3+H4+H5+H6+H7+H8+H9+H10+H11+H12],'-b*')
#        pl.title('Accumulative value of H with increase of the Bessel function order')
#        pl.xlabel('Order (l) of the Bessel function')
#        pl.ylabel('Accumulative value of H')


# Run unittests when module is called as a script
if __name__ == '__main__':
    
    try:
        unittest.main()
    except SystemExit as inst:
        if inst.args[0] is True: # raised by sys.exit(True) when tests failed
            raise
