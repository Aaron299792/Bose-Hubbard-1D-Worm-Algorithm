import sys
import os

sys.path.append(os.path.abspath('../'))

from modules import Hamiltonian1D

#Params 
t = 1.0
U = 4.0
mu = 0.0
nmax = 5

ham = Hamiltonian1D(t, U, mu, nmax)
diag_energy = ham.onsite_energy(3)
nd_energy34 = ham.bosonic_matrix_element(3,4)
nd_energy32 = ham.bosonic_matrix_element(3, 2)

print('-----------------------------')
print('Test of the Hamiltonian class')
print('-----------------------------\n')
print(f'Parameters: \n t = {t} \n U = {U} \n mu = {mu} \n nmax = {nmax} \n')
print(f'Diagonal energy on site: {diag_energy:.3f}')
print(f'Creation 3 -> 4: {nd_energy34:.3f}')
print(f'Anihilation 3 -> 2: {nd_energy32:.3f}')
print(f'Creation 3 -> 5: {ham.bosonic_matrix_element(3,5)}')
print(f'Anihilation 3 -> 1: {ham.bosonic_matrix_element(3,1)}')
print(f'Creation 3 -> 6: {ham.bosonic_matrix_element(3,6)}')
print(f'Anihilation 3 -> 0: {ham.bosonic_matrix_element(3,0)}')
print(f'Anihilation 3 -> -1: {ham.bosonic_matrix_element(3,-1)}')





