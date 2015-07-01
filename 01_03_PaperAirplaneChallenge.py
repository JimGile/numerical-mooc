# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 21:40:48 2015

@author: Jim
"""

from math import sin, cos
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


# set initial conditions
v0 = v_t+0.      # start at the trim velocity (or add a delta)
# theta0 = numpy.pi/-9.0  # initial angle of trajectory
theta0 = numpy.pi  # initial angle of trajectory
x0 = 0.0      # horizotal position is arbitrary
y0 = 100.0   # initial altitude
dt = 0.01

print cos(theta0)

# initialize the array containing the solution for each time-step
u = numpy.empty((1, 4))
u[0] = numpy.array([v0, theta0, x0, y0])  # fill 1st element with initial vals

print u[0]

# time loop - Euler method
y = u[0, 3]
while y > 0.:
    # for n in range(10):
    u[0] = euler_step(u[0], f, dt)
    y = u[0, 3]

print u[0]
