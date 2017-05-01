import numpy as np
import sys
import random as rand
import matplotlib.pyplot as plt
import matplotlib as mpl
import gc
from math import *
from time import time
from numpy.linalg import norm as norm
import scipy.optimize as opt
HUGE = 1e10
TINY = 1e-10

PARALLEL = 1
COLLISION = 1
INPUT = 0
COUETTE = 0
# reflection_type = 'specular'
reflection_type = 'diffuse'

N_steps = 10  # Number of steps
Ne = int(1e5)  # Number of simulated molecules
alpha = 1  # temperature accommodation coef

# Constants (Diatomic Nitrogen)
###########################################
N_a = 6.022e23  # [mol^-1] Avagadro number
k = 1.38064852e-23  # [m^2*kg/s^2/K] Boltzmann constant
molar_mass = 28.0134e-3  # [kg/mol] molar mass
m = molar_mass / N_a  # [kg] mass
d = 4.11e-10  # [m] average diameter

# Size of the box
###########################################
if COLLISION:
    X_max = 0.1  # [m] box size
    Y_max = 0.1  # [m] box size
    N_x = 60
    N_y = 60
    n_subdiv = 2
    Kn = 0.01  # Knudsen number
else:
    X_max = 1  # [m] box size
    Y_max = 1  # [m] box size
    N_x = 100
    N_y = 100
    n_subdiv = 2
    Kn = 100  # Knudsen number

# Calculated size parameters
############################################
L_box = X_max  # [m] characteristic box size
cell_size = X_max / N_x
subcell_size = cell_size / n_subdiv
Z_max = subcell_size  # one layer of subcells in z direction
V = X_max * Y_max * Z_max
V_cell = cell_size ** 2 * Z_max
V_subcell = subcell_size ** 3
boundary = {'down': 0, 'up': X_max, 'left': 0, 'right': Y_max}

# Boundary conditions
############################################
T = 300  # [K] Temperature
T_w = 300  # [K] Wall Temperature

# Calculated parameters
############################################
beta = sqrt(m / 2 / k / T)  # initial beta
c_mean_initial = 2 / sqrt(pi) / beta  # mean speed (needed for time step calculation)
lambd = Kn * L_box  # mean free path
n_initial = 1 / (sqrt(2) * pi * lambd * d ** 2)  # number density in reservoir

c_r_max_initial = 200 # [m/s] initial value(can increase during simulation)

N = n_initial * V  # number of real molecules
Fn = N / Ne  # number of real molecules per simulated one
sigma = pi * d ** 2  # cross section (hard sphere model)

############################################
nu_init = 2 * pi * d ** 2 * n_initial * sqrt(2 * k * T / m)
t_coll_init = 1 / nu_init

N_in_init = int(Ne / N_x / N_y / n_subdiv ** 2)
print ("cell_size=", cell_size)
print ("n_subdiv=", n_subdiv)
print ("subcell_size=", subcell_size)
print ("V_subcell=", V_subcell)
N_coll_init = 1 / (2 * V_subcell) * N_in_init * (N_in_init - 1) * Fn * sigma * c_r_max_initial * t_coll_init
############################################

############################################
time_step = lambd / c_mean_initial  # total time of simulation

############################################
#INPUT
############################################
u_0 = 0
u_couette = 200
A = Z_max*Y_max
s = u_0*beta
N_input_real = time_step*A*n_initial/(2*sqrt(pi)*beta)*(exp(-s**2)+sqrt(pi)*s*(1+erf(s)))
N_output_real = time_step*A*n_initial/(2*sqrt(pi)*beta)*(exp(s**2)-sqrt(pi)*s*(1+erf(-s)))
N_input = int(round(N_input_real/Fn))
N_output = int(round(N_output_real/Fn))