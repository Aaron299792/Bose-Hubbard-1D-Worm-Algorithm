import sys
import os

sys.path.append(os.path.abspath('../'))

from modules import Lattice1D

latt = Lattice1D(10)

print('---------------------------')
print('Test of the Lattice1D class')
print('---------------------------')
print('Number of sites on the lattice: ', latt.get_nsites())
print('Neighbors of each site:', latt._neighbor_table)

