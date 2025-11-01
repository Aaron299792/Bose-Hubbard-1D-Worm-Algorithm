import numpy as np

class Hamiltonian1D:

    def __init__(self, t : float, U : float, mu : float, n_max : int = 50):
        if  t < 0:
            raise ValueError('t must be positive')
        self.t = t
        self.U = U
        self.mu = mu
        self.n_max = n_max

    def onsite_energy(self, n : int):
        """
        Energy eigenvalues corresponding to the diagonal part of the Hamiltonian.
        """

        if n < 0 or n > self.n_max:
            # the occupation value is delimited due to computational resources 
            raise ValueError("occupation out of range")
        
        return 0.5 * self.U * n * (n - 1) - self.mu * n
    
    def bosonic_matrix_element(self, n_from : int, n_to : int):
        """
        Matrix element for a single-site creation/annihilation operators
        """

        if n_from < 0 or n_from > self.n_max or n_to < 0 or n_to > self.n_max:
            return 0.0

        if n_to == n_from - 1:
            # b annihilation in that site
            return np.sqrt(n_from) 
        if n_to == n_from + 1:
            # b^\dagger creates a particle in that site
            return np.sqrt(n_to)
        
        return 0.0
