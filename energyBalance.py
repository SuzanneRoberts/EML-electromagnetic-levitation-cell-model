# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 11:31:27 2015

@author: Suzanne
"""
import numpy as np

# Sample temperature calculation
def sampleTemp(powerIn, sa, fl):
    tol    = 1E-8
    update = 2.0*tol
    T      = 1000 #K - initialize temperature
    Boltz = 5.67E-8
    area = 4.0*np.pi*sa.R**2.0
    while abs(update) > tol:
        Qrad = Boltz*sa.epsilon*area*(T**4.0 - fl.T0**4.0)
        Nu = 2.0 + 0.6 * (2.0*sa.R*fl.v*fl.rho/fl.eta)**0.5 * (fl.Cp*fl.eta/fl.k)**(1/3)
        h = fl.k * Nu / (2.0*sa.R)
        Qconv = h*area*(T - fl.Tf)
        
        R = powerIn - Qrad - Qconv 
        dRdT = - (4.0*Boltz*sa.epsilon*area*T**3.0) - (h*area) 
        
        update = - R/dRdT
        T += update
    
    return T, Qrad, Qconv
