import matplotlib.pyplot as plt
from har_bal import *
import scipy as sp
from scipy import pi,sin,cos
f = 2
omega = 2.*pi * f
numsteps = 11
t = sp.arange(0,1/omega*2*pi,1/omega*2*pi/numsteps)
x = sp.array([sin(omega*t)])
v = sp.array([omega*cos(omega*t)])
states = sp.append(x,v,axis = 0)
state_derives = harmonic_deriv(omega,states)
plt.plot(t,states.T,t,state_derives.T,'x')
# [<matplotlib.line...]
