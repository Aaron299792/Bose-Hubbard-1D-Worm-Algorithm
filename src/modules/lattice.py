import numpy as np

class Lattice1D:
    """
    

    """
    
    def __init__(self, L):
        if L <= 0:
            raise ValueError("L must be a positive")

        self.L = int(L)
        self.nsites = self.L
        
        sites = np.arange(self.L)
        left = (sites - 1) % self.L
        right = (sites + 1) % self.L

        self._neighbor_table = [[ int(left[i]), int(right[i]) ] for i in range(self.L)]
        
    def get_nsites(self):
        return self.nsites

    def get_neighbors(self, site):
        if site < 0 or site >= self.L:
            raise IndexError('site out of range')
        return self._neighbor_table[site]

    
