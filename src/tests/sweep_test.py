import sys
import os

sys.path.append( os.path.abspath('../') )

from modules import TYPE_HOP, TYPE_WORM_HEAD, TYPE_WORM_TAIL, TYPE_WORM_DUMMY
from modules import Hamiltonian1D, Lattice1D, WormConfiguration, WormAlgorithm

#Params
L = 5
nmax = 5
t = 0.01
U = 0.0
mu = 0.0

beta = 10.0
#objects
latt = Lattice1D(L)
ham = Hamiltonian1D(t, U, mu, nmax)
w = WormAlgorithm(latt,ham, beta, n_max = 5, c_worm = 2.0, energy_off = 1.0, seed = 11, epsilon_time = 1.0e-3, kink_prob = 0.4, glue_prob = 0.3, move_prob = 0.3)

sweeps = 10000

for _ in range(sweeps):
    w.monte_carlo_sweep(updates_per_sweep = 100)

rates = w.get_acceptance_rates()


print('------------------')
print('Test of simulation')
print('------------------\n')
print(f'Parameters: \n L = {L} \n t = {t} \n U = {U} \n mu = {mu} \n nmax = {nmax} \n beta = {beta} \n site of test = {3} \n')
print(f'Statistics: \n Insert {rates['insert']} \n Glue {rates['glue']} \n Move {rates['move']} \n InsertKink {rates['insertkink']} \n DeleteKink{rates['deletekink']}')
