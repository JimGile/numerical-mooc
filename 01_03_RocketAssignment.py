# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11, 2015

@author: Jim
"""

from math import sin, cos, log, ceil
import numpy
from matplotlib import pyplot
#%matplotlib inline
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

# model parameters:
ms = 50.     # in kg is the weight of the rocket shell
g = 9.81     # gravity in m s^{-2}
rho = 1.091  # in kgm3 is the avg air density (assumed constant for flight)
r = 0.5      # in meteres is the maximum cross sectional radius of the rocket
A = numpy.pi*r**2        # is the maximum cross sectional area of the rocket
ve = 325.    # in ms is the exhaust speed
Cd = 0.15    # is the drag coefficient

# set initial conditions
t0 = 0.      # in seconds is t=0
mp0 = 100.   # in kg at time t=0 is the initial weight of rocket propellant
br0 = 20.    # in kg/s at time t=0 is the initial burn rate
v0 = 0.0     # in m/s at time t=0 is the initial velocity
h0 = 0.0     # in m at time t=0 is the initial height

# set time parameters
T = 37.5     # final time
dt = 0.1     # time increment

def burnRate(t):
    """Returns the burn rate at the given time t.
    Assumes constant burn rate of br0 until all propellant mp0 has benn burned.

    Parameters
    ----------
    t : float
        time

    Returns
    -------
    brt : float
        the burn rate at the given time t.
    """
    brt = (br0 if (br0*(t+dt) < mp0) else 0.)
    return brt


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

    t = u[0]
    mp = u[1]
    v = u[2]
    h = u[3]
    brt = burnRate(t)
    return numpy.array([1,
                        brt*(-1),
                        ((brt*ve)/(ms+mp))-((rho*v*abs(v)*A*Cd)/(2*(ms+mp)))-g,
                        v])


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


N = int(T/dt) + 1                  # number of time-steps
# t = numpy.linspace(0.0, T, N)      # time discretization

# initialize the array containing the solution for each time-step
u = numpy.empty((N, 4))
u[0] = numpy.array([t0, mp0, v0, h0])  # fill 1st element with initial vals

# time loop - Euler method
for n in range(N-1):
    u[n+1] = euler_step(u[n], f, dt)

# print u
maxV = numpy.max(u[:, 2])    # Max velocity
maxH = numpy.max(u[:, 3])    # Max height
print u[numpy.where(u[:, 2] == maxV)]
print u[numpy.where(u[:, 3] == maxH)]
print u[numpy.where(u[:, 3] < 0)]

# get the rocket's position with respect to the time
time = u[:, 0]
height = u[:, 3]

# visualization of the path
pyplot.figure(figsize=(8, 6))
pyplot.grid(True)
pyplot.xlabel(r'time', fontsize=18)
pyplot.ylabel(r'height', fontsize=18)
pyplot.title('Rocket trajectory, flight time = %.2f' % T, fontsize=18)
pyplot.plot(time, height, 'k-', lw=2)
