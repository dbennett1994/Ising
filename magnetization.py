#!/usr/bin/env python

from __future__ import division #necessary to perform division of ints effectively cast into floats before division

import sys
import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand

import IsingClass #module written to carry out the required calculations (metropolis, neighbours, energy, magnetization, etc.)


def main(J):
    steps = 2000 #steps for equilibrium
    temp_steps = 40 #step size for increase in temperature [J/k]

    kT = np.linspace(1.0,10.0,temp_steps) #temperature values

    #arrays for future plotting
    magnetization = np.zeros(temp_steps) 
    Energy = np.zeros(temp_steps)
    susceptibility = np.zeros(temp_steps)
    specific_heat = np.zeros(temp_steps)
    mag_square = np.zeros(temp_steps) 
    e_square = np.zeros(temp_steps) 

    for l in range(len(kT)):
        temperature = kT[l]

        #initialise observables to zero
        M = 0
        E = 0
        M_sq = 0
        E_sq = 0

        #reset to random initial configuration
        lattice = model.spin_config(N,N)  

        #set out to reach equilibrium
        for t in range(steps):
            model.metropolis(lattice, temperature,N,N,J,h)

        #collect statistics
        stat_steps = int(steps/2) 
        for k in range(stat_steps):  

            model.metropolis(lattice, temperature,N,N,J,h) 

            mag = model.mag_per_spin(lattice,N,N)
            mag2 = model.mag2_per_spin(lattice,N,N)
            engy = model.energy(lattice,N,N,J)
            engy2 = model.energy2(lattice,N,N,J)

            M += mag
            E += engy
            M_sq += mag2
            E_sq += engy2

        #observables, after equilibrium is achieved, normalised by sweeps
        magnetization[l] = abs(M)/stat_steps 
        Energy[l] = E/stat_steps 
        mag_square[l] = M_sq/stat_steps
        e_square[l] = E_sq/stat_steps

    #perform manipulation directly on array instead of loop
    susceptibility = (mag_square - magnetization**2)/kT
    specific_heat = (e_square - Energy**2)/(kT**2)
    plt.plot(kT, magnetization, '--o', label = "J = %.2f"%J)
    plt.legend()    


plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.title("Average spin site magentization versus temperature")
plt.xlabel("$T$ $[J/k_B]$")
plt.ylabel("$|<M>|$ $[\mu]$")
plt.grid()



N = 16 #grid dimensions
h = 0.0

J_vals = []
args = len(sys.argv)

for i in range(2, args,1):
    J_vals.append(float(sys.argv[i]))

model = IsingClass.Ising() #constructor

for value in J_vals:
    J = value
    main(J)



plt.show()


