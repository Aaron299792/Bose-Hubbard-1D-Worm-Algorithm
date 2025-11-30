import sys
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MODULES_PATH = os.path.abspath(os.path.join(BASE_DIR, "..", "modules"))

if MODULES_PATH not in sys.path:
    sys.path.insert(0, MODULES_PATH)

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
