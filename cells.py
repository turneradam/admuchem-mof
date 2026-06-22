import sys
from math import ceil
import numpy as np
from ase.io import read

rcut = 12.0
for path in sys.argv[1:]:
    c = read(path).cell.array
    v = abs(np.linalg.det(c))
    w = [v / np.linalg.norm(np.cross(c[(i+1)%3], c[(i+2)%3])) for i in range(3)]
    n = [ceil(2 * rcut / x) for x in w]
    print(f"{path.split('/')[-1]:32s} {n[0]} {n[1]} {n[2]}   {np.prod(n):>3d}x")
