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


def dist(v0, theta0, x0, y0, dt):
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

# set initial conditions
x0 = 0.0    # horizotal position is arbitrary
y0 = 10.0   # initial altitude
dt = 0.01   # time increment

maxDist = numpy.empty((2, 4))
for i in range(-5, 6, 1):
    v0 = v_t + i
    for j in range(-1, 2, 1):
        theta0 = (j*numpy.pi)/16
        calcDist = dist(v0, theta0, x0, y0, dt)
        x = calcDist[1, 2]
        if x > maxDist[1, 2]:
            maxDist = calcDist[:]

print maxDist

T = 15.0                           # final time
N = int(T/dt) + 1                  # number of time-steps
t = numpy.linspace(0.0, T, N)      # time discretization

# initialize the array containing the solution for each time-step
u = numpy.empty((N, 4))
u[0] = maxDist[0]  # fill 1st element with initial vals

print u[0]

# time loop - Euler method
for n in range(N-1):
    u[n+1] = euler_step(u[n], f, dt)

# get the glider's position with respect to the time
x = u[:, 2]
y = u[:, 3]

# visualization of the path
pyplot.figure(figsize=(8, 6))
pyplot.grid(True)
pyplot.xlabel(r'x', fontsize=18)
pyplot.ylabel(r'y', fontsize=18)
pyplot.suptitle('Distance = %.2f' % maxDist[1, 2], fontsize=18)
pyplot.title('Glider trajectory, flight time = %.2f' % T, fontsize=14)
pyplot.plot(x, y, 'k-', lw=2)
