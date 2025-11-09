import sys 
import os

sys.path.append(os.path.abspath("../modules"))

print("---------------------------------------------")
print("----------Importing test validation----------")
print("---------------------------------------------")
try:
    from lattice import Lattice1D
    print("--lattice.py imported correctly")
    
    from hamiltonian import Hamiltonian1D
    print("--hamiltonian.py imported correctly")

    from configuration import WormConfiguration
    print("--configuration.py imported correctly")

    from worm import WormAlgorithm
    print("--worm.py imported correctly")

    print("-------Importing was succesful----------")

except Exception:
    print("modules couldn't be imported correctly")
