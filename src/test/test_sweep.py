import os 
import sys

sys.path.append(os.path.abspath("../modules"))

from lattice import Lattice1D
from hamiltonian import Hamiltonian1D
from worm import WormAlgorithm
from configuration import WormConfiguration

#Parameters
L = 2
t = 1.0
U = 50.0
mu = 1.0
beta = 1.0
sweeps = 10000
updates_per_sweep = 100000


print("------------Constructing objects-----------")
print(f"L ={L}, t={t}, U={U}, mu={mu}, beta={beta}")
lattice = Lattice1D(L)
hamiltonian = Hamiltonian1D(t, U, mu, n_max=4)
worm = WormAlgorithm(lattice, hamiltonian, beta, n_max=4, c_worm= 1000.0, energy_off= 1.0, seed=11)

"""
for sweep in range(sweeps):
    worm.monte_carlo_sweep(updates_per_sweep)
    rates = worm.get_acceptance_rates()
    print(f"Sweep {sweep + 1}: "
          f"Z-sector={worm.config.in_z_sector}, "
          f"Insert={rates['insert']}, "
          f"Glue={rates['glue']}, "
          f"Move={rates['move']}, ")
"""
if worm.stats["glue_accepts"] > 0:
    print("Worm glued at sweep {sweep + 1}, total glued: {worm.stats['glue_accepts']}")
else: 
    print("No glues encounter")
"""   
for site in range(lattice.get_nsites()):
    events = worm.config.dump_site(site)
    print(f"\nSite {site} events: ")
    for event in events:
        print(f"  t={event['time']:.6f}, type={event['type']}, "
              f"occL={event['occ_left']}, occR={event['occ_right']}")
"""
