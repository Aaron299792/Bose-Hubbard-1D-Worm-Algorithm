# src/test/test_configuration.py
import os
import sys
import unittest

import numpy as np

MODULES_DIR = os.path.join(os.path.dirname(__file__), "..", "modules")
sys.path.append(MODULES_DIR)

import lattice
import hamiltonian
import configuration


class TestConfiguration(unittest.TestCase):
    def setUp(self):
        self.lat = lattice.Lattice1D(4)
        self.H = hamiltonian.Hamiltonian1D(t=1.0, U=2.0, mu=0.5, n_max=5)
        self.beta = 3.0

    def test_initial_configuration(self):
        conf = configuration.WormConfiguration(
            self.lat, self.H, self.beta, initial_occupation=2
        )

        # Single dummy event per site at time = beta
        self.assertEqual(conf.nsites, self.lat.get_nsites())
        for site in range(conf.nsites):
            evs = conf.events[site]
            self.assertEqual(len(evs), 1)
            ev = evs[0]
            self.assertAlmostEqual(ev["time"], self.beta)
            self.assertEqual(ev["occ_left"], 2)
            self.assertEqual(ev["occ_right"], 2)

        # Occupation is constant in imaginary time
        for t in np.linspace(0.0, self.beta, 5):
            for site in range(conf.nsites):
                self.assertEqual(conf.get_occupation_at_time(site, float(t)), 2)

        # Should not raise
        conf.validate_configuration()

    def test_insert_and_remove_element_preserve_order(self):
        conf = configuration.WormConfiguration(
            self.lat, self.H, self.beta, initial_occupation=1
        )
        site = 0

        # Insert a hop event at tau = beta/2
        idx = conf.insert_element(
            site,
            self.beta / 2.0,
            configuration.TYPE_HOP,
            occ_left=1,
            occ_right=2,
            linked_site=1,
        )

        # Times at this site must be sorted
        times = [e["time"] for e in conf.events[site]]
        self.assertEqual(times, sorted(times))

        # Removing the element returns us to the original config
        conf.remove_element(site, idx)
        self.assertEqual(len(conf.events[site]), 1)

        conf.validate_configuration()


if __name__ == "__main__":
    unittest.main()

