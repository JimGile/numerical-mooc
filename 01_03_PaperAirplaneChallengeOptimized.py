# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 21:40:48 2015

@author: Jim
"""

from math import sin, cos
from scipy.optimize import minimize
import numpy
from matplotlib import pyplot
#%matplotlib inline
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

# model parameters:
g = 9.8      # gravity in m s^{-2}
v_t = 4.9   # trim velocity in m s^{-1}
C_D = 1/5.  # drag coefficient --- or D/L if C_L=1
C_L = 1.0    # for convenience, use C_L = 1


def f(u):
    """Returns the right-hand side of the phugoid system of equations.

    Parameters
    ----------
    u : array of float
        array containing the solution at time n.

    Returns
    -------
    dudt : array of float
        array containing the RHS given u.
    """

    v = u[0]
    theta = u[1]
    return numpy.array([-g*sin(theta) - C_D/C_L*g/v_t**2*v**2,
                        -g*cos(theta)/v + g/v_t**2*v,
                        v*cos(theta),
                        v*sin(theta)])


def euler_step(u, f, dt):
    """Returns the solution at the next time-step using Euler's method.

    Parameters
    ----------
    u : array of float
        solution at the previous time-step.
    f : function
        function to compute the right hand-side of the system of equation.
    dt : float
        time-increment.

    Returns
    -------
    u_n_plus_1 : array of float
        approximate solution at the next time step.
    """

    return u + dt * f(u)


def dist(v0, theta0):
    # set initial conditions
    x0 = 0.0    # horizotal position is arbitrary
    y0 = 10.0   # initial altitude
    dt = 0.01   # time increment

    # initialize the array containing the solution for each time-step
    u = numpy.empty((2, 4))
    u[0] = numpy.array([v0, theta0, x0, y0])  # initial values
    u[1] = numpy.array([v0, theta0, x0, y0])  # solution

    # stop calculating distance when y <= 0
    y = u[1, 3]
    while y > 0.:
        # for n in range(10):
        u[1] = euler_step(u[1], f, dt)
        y = u[1, 3]

    return u


def invDist(v0, theta0):
    calcDist = dist(v0, theta0)
    return 1./calcDist[1, 2]


maxDist = numpy.empty((2, 4))
for i in range(-5, 6, 1):
    v0 = v_t + i
    for j in range(-1, 2, 1):
        theta0 = (j*numpy.pi)/16
        calcDist = dist(v0, theta0)
        x = calcDist[1, 2]
        if x > maxDist[1, 2]:
            maxDist = calcDist[:]

print maxDist

# now use the scipy optimizer
cons = ({'type': 'ineq', 'fun': lambda x:  x[0]},
        {'type': 'ineq', 'fun': lambda x:  invDist(x[0], x[1])},
        )

optResult = minimize(lambda x: invDist(x[0], x[1]), [v_t, 0.0],
                     method='COBYLA',
                     constraints=cons)

print optResult
print optResult.x

maxDist2 = dist(optResult.x[0], optResult.x[1])
print maxDist2
