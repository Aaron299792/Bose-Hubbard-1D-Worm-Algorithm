# src/test/test_lattice_hamiltonian.py
import os
import sys
import math
import unittest

import numpy as np

# Add src/modules to path so we can `import lattice, hamiltonian, ...`
MODULES_DIR = os.path.join(os.path.dirname(__file__), "..", "modules")
sys.path.append(MODULES_DIR)

import lattice
import hamiltonian


class TestLatticeHamiltonian(unittest.TestCase):
    def test_lattice_neighbors_periodic(self):
        L = 4
        lat = lattice.Lattice1D(L)
        self.assertEqual(lat.get_nsites(), L)

        # Periodic 1D neighbors
        self.assertEqual(lat.get_neighbors(0), [L - 1, 1])
        self.assertEqual(lat.get_neighbors(1), [0, 2])
        self.assertEqual(lat.get_neighbors(3), [2, 0])

        # Out-of-range sites
        with self.assertRaises(IndexError):
            lat.get_neighbors(-1)
        with self.assertRaises(IndexError):
            lat.get_neighbors(L)

    def test_lattice_invalid_L(self):
        with self.assertRaises(ValueError):
            lattice.Lattice1D(0)
        with self.assertRaises(ValueError):
            lattice.Lattice1D(-5)

    def test_onsite_energy_bose_hubbard(self):
        t = 1.0
        U = 3.0
        mu = 0.7
        H = hamiltonian.Hamiltonian1D(t, U, mu, n_max=10)

        # Check Bose–Hubbard onsite energy:
        # E(n) = (U/2) n (n-1) - mu n
        for n in range(0, 11):
            e = H.onsite_energy(n)
            expected = 0.5 * U * n * (n - 1) - mu * n
            self.assertAlmostEqual(e, expected)

        # Occupation bounds enforced
        with self.assertRaises(ValueError):
            H.onsite_energy(-1)
        with self.assertRaises(ValueError):
            H.onsite_energy(11)

    def test_bosonic_matrix_element(self):
        H = hamiltonian.Hamiltonian1D(t=1.0, U=0.0, mu=0.0, n_max=5)

        # Check <n_to| b |n_from> and <n_to| b^\dagger |n_from>
        for n_from in range(0, 6):
            for n_to in range(0, 6):
                val = H.bosonic_matrix_element(n_from, n_to)
                if n_to == n_from - 1:
                    # annihilation
                    self.assertAlmostEqual(val, math.sqrt(n_from))
                elif n_to == n_from + 1:
                    # creation
                    self.assertAlmostEqual(val, math.sqrt(n_to))
                else:
                    self.assertAlmostEqual(val, 0.0)

        # Out of bounds -> 0
        self.assertEqual(H.bosonic_matrix_element(-1, 0), 0.0)
        self.assertEqual(H.bosonic_matrix_element(0, 6), 0.0)


if __name__ == "__main__":
    unittest.main()

