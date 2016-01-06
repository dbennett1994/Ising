from __future__ import division #necessary to perform division of ints effectively cast into floats before division

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand

#Ising model basic class
class Ising():
    
    #NOTE: decouple the N's. Nx and Ny to be specified explicitly to allow for application to 1D by setting Ny=1 or Nx=1

    def spin_config(self,Nx,Ny):
        lattice = rand.choice([1,-1],[Nx,Ny]) 
        return lattice
    
    #plot the spin configuration in the lattice 
    def lattice_plt(self, lattice):
        plt.matshow(lattice, cmap=plt.cm.gray)
        #LaTeX
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        plt.title("Spin configuration in the lattice")
        plt.axis('off')
        plt.show()
    
    def nb_sum(self,lattice,i,j,Nx,Ny):
        return lattice[(i+1) % Nx][j] + lattice[(i-1) % Nx][j] + lattice[i][(j+1) % Ny] +  lattice[i][(j-1) % Ny]
    
    #delta_energy with periodic BC imposed using the modulo operator, per spin site
    def dEnergy(self, lattice,i,j,Nx,Ny,J): 
        nbs = self.nb_sum(lattice, i, j,Nx,Ny)
        return 2*J*lattice[i][j]*nbs
    
    #metropolis with periodic BC
    def metropolis(self, lattice, kT,Nx,Ny,J):
        for i in range(Nx):
            for j in range(Ny):
                delta_E = self.dEnergy(lattice,i ,j,Nx,Ny,J)
                site = lattice[i][j]
                if delta_E<0:
                    site *= -1 #energy is lowered, accept
                elif np.exp(-delta_E/kT) > rand.rand(): #COMMMENT ON RAND(); uniform distribution in [0,1)
                    site *= -1 #accept according to Boltzman dist.
                lattice[i][j] = site
        return lattice
    
    #energy per site
    def energy(self, lattice,Nx,Ny,J):
        erg = 0
        for i in range(Nx):
            for j in range(Ny):
                Snb = self.nb_sum(lattice,i,j,Nx,Ny)
                erg += -J*lattice[i][j]*Snb
        norm_energy = erg/(2*(Nx*Ny)) #divide by twice number of spins to avoid overcounting
        return norm_energy
    
    #energy squared per site                       
    def energy2(self, lattice,Nx,Ny,J):
        erg2 = 0
        for i in range(Nx):
            for j in range(Ny):
                #sum of neighbouring spins
                sum_nb2 = self.nb_sum(lattice,i,j,Nx,Ny)
                #energy per site (not normalised)
                erg2 += (-J*lattice[i][j]*sum_nb2)**2
        norm_energy2 = erg2/(4*Nx*Ny) 
        return norm_energy2
    
    #magnetization per spin
    def mag_per_spin(self,lattice,Nx,Ny):
        Mag = 0
        for i in range(Nx):
            for j in range(Ny):
                Mag += lattice[i][j]
        return Mag/(Nx*Ny)
                       
    #magnetization squared per spin
    def mag2_per_spin(self,lattice,Nx,Ny):
        Mag2 = 0
        for i in range(Nx):
            for j in range(Ny):
                Mag2 += lattice[i][j]**2
        return Mag2/(Nx*Ny)
    
    #look at ground state spin energy with squared (1D line) lattice
    def square_energy_ground_state(self, Nx,Ny,J):
        #REFERENCE: http://physics.stackexchange.com/questions/133005/how-to-calculate-the-ground-state-energy-for-the-ising-model
        
        #square lattice
        lattice_groud_state = np.ones((Nx,Ny), dtype=np.int) #int, save memory - ground state is given by all spins collinear.
        #given the symmetry, simply choose all spins up

        lattice_nb_sum = []

        for i in range(N):
            for j in range(N):
                nbs_sum = lattice_groud_state[i][j]*self.nb_sum(lattice_ground_state,i,j,Nx, Ny) #projection of spin site onto sum of all edges (topological equivalence of graph and its dual in contribution to energy, Statistical Physics II)
                lattice_nb_sum.append(nbs_sum) #system ground state energy x 2

        ground_energy = (np.sum(lattice_nb_sum)*(-J))/(2.0*N*N) #spins can be up or down; energy per site

        E_ground = np.int(ground_energy)
        print "Ground state energy = %dNJ, for N the total number of spin sites"%E_ground 
        
        
    

  

   



