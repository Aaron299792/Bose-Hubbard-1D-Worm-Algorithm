import sys
import os

sys.path.append( os.path.abspath('../') )

from modules import TYPE_HOP, TYPE_WORM_HEAD, TYPE_WORM_TAIL, TYPE_WORM_DUMMY
from modules import Hamiltonian1D, Lattice1D, WormConfiguration

#Params
L = 5
nmax = 5
t = 1.0
U = 4.0
mu = 0.5
beta = 2.0
site_test = 3
#objects
latt = Lattice1D(L)
ham = Hamiltonian1D(t, U, mu, nmax)
cfg = WormConfiguration(latt, ham, beta)
#times
time1 = 5.5
time2 = 3.2
time3 =5.502
neighbors = latt.get_neighbors(site_test)


print('-------------------------------')
print('Test of the Configuration class')
print('-------------------------------\n')
print(f'Parameters: \n L = {L} \n t = {t} \n U = {U} \n mu = {mu} \n nmax = {nmax} \n beta = {beta} \n site of test = {3} \n')
print(f'Test times: time_1 = {time1}, time_2 = {time2}, time_3 = {time3}')
print(f'Mapped time_1 [0, beta[ : {cfg._norm_time(time1):.3f}')
print(f'Mapped time_2 [0, beta[ : {cfg._norm_time(time2):.3f}')
print(f'Time distance in mapping t1 - t2: {cfg._time_distance(time1, time2):.3f}')
print(f'Time distance in mapping t2 - t1: {cfg._time_distance(time2, time1):.3f}')
print(f'Time Close between t1, t2: {cfg._time_close(time1, time2, epsilon = 0.01):.3f}')
print(f'Time Close between t1, t3: {cfg._time_close(time1, time3, epsilon = 0.01):.3f}')
events = cfg._event_times(site_test)
print(f'Time events in site {site_test} : {events}')
