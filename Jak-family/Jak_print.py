import sys
import time

from sympy import QQ

sys.path.insert(0, "../")
sys.path.insert(0, "./../../")
import parser
import clue
from sparse_polynomial import SparsePolynomial

system = parser.read_system("Barua.ode")

obs_sets = [
    ["aS000"],
    #["aS027"]
]

m = 0
for obs_set in obs_sets:
    m += 1
    print("===============================================")
    obs_polys = [SparsePolynomial.from_string(s, system['variables']) for s in obs_set]
    lumped = clue.do_lumping(system['equations'], obs_polys)
    subspaces = lumped["subspace"]
    polynomial = lumped["polynomials"]
    with open('obs'+str(m)+'/original_macrovariables.txt','w') as afile:
        for lines in subspaces:
            for i in lines:
                print(i,end=" ",file=afile)
            print(file=afile)
    with open('obs'+str(m)+'/original_ODE.txt','w') as afile:
        for lines in polynomial:
            print(lines,file = afile )
