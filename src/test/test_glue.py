import os
import sys

sys.path.append(os.path.abspath("../modules"))

from lattice import Lattice1D
from hamiltonian import Hamiltonian1D
from worm import WormAlgorithm
from configuration import WormConfiguration

# Instantiate simple system
L = Lattice1D(1)
H = Hamiltonian1D(t=0.0, U=1.0, mu=0.0, n_max=5)
WA = WormAlgorithm(lattice=L, hamiltonian=H, beta=1.0,
                        n_max=5, c_worm=0.5, seed=123,
                        epsilon_time=1e-3, kink_prob=0.0,
                        glue_prob=0.3, move_prob=0.7)

# Run few sweeps with original glue_worm
for _ in range(200):
    WA.monte_carlo_sweep(50)
WA.get_acceptance_rates(), WA.stats
