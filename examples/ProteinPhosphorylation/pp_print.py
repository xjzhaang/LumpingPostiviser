import sys
import time

from sympy import QQ

sys.path.insert(0, "../")
sys.path.insert(0, "./../../")
import parser
import clue
from sparse_polynomial import SparsePolynomial


for n in range(2,6):

    system = parser.read_system(f"e{n}.ode")

    obs = [
        clue.SparsePolynomial.from_string("S0", system["variables"]),
        clue.SparsePolynomial.from_string("S1", system["variables"])
    ]

    #start = time.time()
    lumped = clue.do_lumping(system['equations'], obs)
    subspaces = lumped["subspace"]
    polynomial = lumped["polynomials"]
    with open('obs'+str(n)+'/original_macrovariables.txt','w') as afile:
        for lines in subspaces:
            for i in lines:
                print(i,end=" ",file=afile)
            print(file=afile)
    with open('obs'+str(n)+'/original_ODE.txt','w') as afile:
        for lines in polynomial:
            print(lines,file = afile )
            
