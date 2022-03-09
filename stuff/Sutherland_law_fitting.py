#!/usr/bin/python

import numpy as np
from scipy.optimize import curve_fit


### Introduce values of T [K], mu [Pa s] and kappa [W/(m K)] from NIST: https://webbook.nist.gov/chemistry/fluid/
T     = [300.0, 350.0, 400.0, 450.0]
mu    = [1.19e-5, 1.38e-5, 1.57e-5, 1.76e-5]
kappa = [0.0135, 0.0175, 0.0215, 0.0255]


### Obtain reference values
T_0     = T[0]		# Reference temperature [K]
mu_0    = mu[0]		# Reference dynamic viscosity [Pa s]
kappa_0 = kappa[0]	# Reference thermal conductivity [W/(m K)]


### Define Sutherland's Law for dynamic viscosity
def mu_sutherland( T, S_mu ):
    return mu_0*( ( T/T_0 )**1.5 )*( ( T_0 + S_mu )/( T + S_mu ) )
### Define Sutherland's Law for thermal conductivity
def kappa_sutherland( T, S_kappa ):
    return kappa_0*( ( T/T_0 )**1.5 )*( ( T_0 + S_kappa )/( T + S_kappa ) )


### Fit data to Sutherland's Law for dynamic viscosity to obtain S_mu
popt, pcov = curve_fit( mu_sutherland, T, mu )
S_mu = popt[0]
### Fit data to Sutherland's Law for thermal conductivity to obtain S_kappa
popt, pcov = curve_fit( kappa_sutherland, T, kappa )
S_kappa = popt[0]


### Print value of fitted parameters
print( 'mu_0 = ' + str( mu_0 ) + ', ' + 'kappa_0 = ' + str( kappa_0 ) + ', T_0 = ' + str( T_0 ) + ', ' + 'S_mu = ' + str( S_mu ) + ', ' + 'S_kappa = ' + str( S_kappa ) )
