import sys
import os

sys.path.insert(0, "../modules")
import math

from lattice import Lattice1D
from hamiltonian import Hamiltonian1D
from worm import WormAlgorithm
from configuration import (
    WormConfiguration,
    TYPE_WORM_HEAD,
    TYPE_WORM_TAIL,
    TYPE_HOP,
)

# Simple deterministic RNG stub so probabilities do not get in the way
class StubRng:
    def random(self):
        return 0.0  # always accept when allowed

    def choice(self, arr):
        return arr[0]  # always pick the first neighbour


def test_kinks():
    # 1D lattice with 2 sites
    L = Lattice1D(2)
    H = Hamiltonian1D(t=1.0, U=0.0, mu=0.0, n_max=5)
    beta = 2.0

    wa = WormAlgorithm(
        L, H, beta,
        n_max=5,
        c_worm=1.0,
        energy_off=1.0,
        seed=1,
        epsilon_time=1.0e-3,
        kink_prob=0.5,
        glue_prob=0.0,
        move_prob=0.0,
        verbose=False,
    )

    # Replace RNG with deterministic stub for this test
    wa.rng = StubRng()

    # Build configuration with occupation 2 on all sites
    wa.config = WormConfiguration(L, H, beta, initial_occupation=2)

    tau = 0.5

    # Insert worm head and tail manually on site 0
    idx_head = wa.config.insert_element(
        0, tau, TYPE_WORM_HEAD,
        occ_left=2, occ_right=3, linked_site=0
    )
    idx_tail = wa.config.insert_element(
        0, tau + 1e-3, TYPE_WORM_TAIL,
        occ_left=3, occ_right=2, linked_site=0
    )

    wa.config.worm_head_site = 0
    wa.config.worm_head_time = wa.config.events[0][idx_head]['time']
    wa.config.worm_head_wpm  = 0

    wa.config.worm_tail_site = 0
    wa.config.worm_tail_time = wa.config.events[0][idx_tail]['time']
    wa.config.worm_tail_wpm  = 0

    wa.config.in_z_sector = False

    print("Initial occupations at tau:",
          wa.config.get_occupation_at_time(0, tau),
          wa.config.get_occupation_at_time(1, tau))

    # ---- INSERT KINK ----
    ok_insert = wa.insert_kink()
    print("insert_kink returned:", ok_insert)
    print("insertkink_accepts:", wa.stats['insertkink_accepts'])

    n_hops = sum(
        1
        for s in range(wa.config.nsites)
        for ev in wa.config.events[s]
        if ev['type'] == TYPE_HOP
    )
    print("TYPE_HOP events after insert:", n_hops)

    print("worm head site/time after insert:",
          wa.config.worm_head_site, wa.config.worm_head_time)
    print("worm tail site/time after insert:",
          wa.config.worm_tail_site, wa.config.worm_tail_time)

    # ---- DELETE KINK ----
    ok_delete = wa.delete_kink()
    print("delete_kink returned:", ok_delete)
    print("deletekink_accepts:", wa.stats['deletekink_accepts'])

    n_hops2 = sum(
        1
        for s in range(wa.config.nsites)
        for ev in wa.config.events[s]
        if ev['type'] == TYPE_HOP
    )
    print("TYPE_HOP events after delete:", n_hops2)


if __name__ == "__main__":
    test_kinks()

