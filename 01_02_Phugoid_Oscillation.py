# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 14:12:54 2015

@author: Jim
"""

import numpy
from matplotlib import pyplot
#%matplotlib inline

# initial conditions
z0 = 100.  # altitude
b0 = 10.  # upward velocity resulting from gust
zt = 100.
g = 9.81


def numerical_solution(t):
    """Uses Euler's method to return an array of numerical solution
    values corresponding to the given t array.

    Parameters
    ----------
    t : array of float
        numerical times.

    Returns
    -------
    z : array of float
        numerical solution of changing elevation values corresponding
        to the given t array.
    """

    N = len(t)
    dt = (t[N-1] - t[0])/(N-1)
    u = numpy.array([z0, b0])

    # initialize an array to hold the changing elevation values
    z = numpy.zeros(N)
    z[0] = z0

    # time-loop using Euler's method
    for n in range(1, N):
        u = u + dt*numpy.array([u[1], g*(1-u[0]/zt)])
        # print u
        z[n] = u[0]

    return z


def exact_solution(t):
    """Returns an array of the exact solution values corresponding
    to the given t array.

    Parameters
    ----------
    t : array of float
        numerical times.

    Returns
    -------
    z_exact : array of float
        exact solution values corresponding to the given t array.
    """
    z_exact = b0*(zt/g)**.5*numpy.sin((g/zt)**.5*t)+\
                (z0-zt)*numpy.cos((g/zt)**.5*t)+zt

    return z_exact


def get_error(z, z_exact, dt):
    """Returns the error relative to analytical solution
    using L-1 norm.

    Parameters
    ----------
    z : array of float
        numerical solution.
    z_exact : array of float
        analytical solution.
    dt : float
        time increment.

    Returns
    -------
    err : float
        L_{1} norm of the error with respect to the exact solution.
    """

    return dt * numpy.sum(numpy.abs(z-z_exact))

T = 10.0
dt = 1.0
N = int(T/dt)+1
t = numpy.linspace(0.0, T, N)

# time-loop using Euler's method
z = numerical_solution(t)
print z

# calculate exact solution for same time frame
z_exact = exact_solution(t)
# print z_exact

# time-increment array
dt_values = numpy.array([0.1, 0.05, 0.01, 0.005, 0.001, 0.0005])
error_values = numpy.empty_like(dt_values)

for i, dt in enumerate(dt_values):
    N = int(T/dt)+1    # number of time-steps
    # discretize the time using numpy.linspace() ###
    t = numpy.linspace(0.0, T, N)

    # time-loop using Euler's method
    z = numerical_solution(t)
    # calculate exact solution for same time frame
    z_exact = exact_solution(t)

    # call the function get_error() ###
    error_values[i] = get_error(z, z_exact, dt)

pyplot.figure(figsize=(15,8))   #set plot size
pyplot.ylim(40,160)             #y-axis plot limits
pyplot.tick_params(axis='both', labelsize=14) #increase font size for ticks
pyplot.xlabel('t', fontsize=14) #x label
pyplot.ylabel('z', fontsize=14) #y label
pyplot.plot(t,z)
pyplot.plot(t, z_exact)
pyplot.legend(['Numerical Solution','Analytical Solution']);

pyplot.figure(figsize=(10, 6))
pyplot.tick_params(axis='both', labelsize=14) #increase tick font size
pyplot.grid(True)                         #turn on grid lines
pyplot.xlabel('$\Delta t$', fontsize=16)  #x label
pyplot.ylabel('Error', fontsize=16)       #y label
pyplot.loglog(dt_values, error_values, 'ko-')  #log-log plot
pyplot.axis('equal')                      #make axes scale equally;