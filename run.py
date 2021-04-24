import sys
import os
import time

from sympy import QQ

sys.path.insert(0, "CLUE/")

import parser
import clue
from sparse_polynomial import SparsePolynomial

if len(sys.argv) < 3:
    print("Usage: python3 run.py path_to_ode_file path_to_obs_file [path_to_julia_binary]")
    exit()

julia_call = "julia"
if len(sys.argv) == 4:
    julia_call = sys.argv[3]

# Reading the system
system = parser.read_system(sys.argv[1])
model = system['name']

# Reading the observables
with open(sys.argv[2]) as f:
    lines = f.readlines()
    obs = [SparsePolynomial.from_string(line, system['variables']) for line in lines]

lumped = clue.do_lumping(system['equations'], obs, print_reduction=False)
subspaces = lumped["subspace"]
polynomial = lumped["polynomials"]
matrix_filename = model + '_CLUE_matrix.txt'
ode_filename = model + '_CLUE_ode.txt'

with open(matrix_filename,'w') as afile:
    for lines in subspaces:
        for i in lines:
            print(i,end=" ",file=afile)
        print(file=afile)

with open(ode_filename,'w') as afile:
    for lines in polynomial:
        print(lines,file = afile)

os.system(f"{julia_call} src/LumpPositive.jl {matrix_filename} {ode_filename} {model}")
print("Computed new matrix and ODE \n")

os.remove(matrix_filename)
os.remove(ode_filename)

with open(model + '_final_macrovariables.txt') as f:
    lines = f.readlines()

with open(model + '_final_macrovariables.txt', 'w') as afile:
    for index in range(len(lines)):
        string = 'y' + str(index) + ' = '
        line = lines[index].split(" ")[: -1]
        for l in range(len(line)):
            if line[l] == '1':
                string += system["variables"][l] + ' + '
            elif line[l] != '0':
                string += line[l]
                string += system["variables"][l] + ' + '
        string = string[: -3]
        print(string,file = afile)
print("Printed new macrovariables \n")
