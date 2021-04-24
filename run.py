import sys
import os
import time

from sympy import QQ

sys.path.insert(0, "CLUE/")
#sys.path.insert(0, "CLUE/examples/ProteinPhosphorylation/")

import parser
import clue
from sparse_polynomial import SparsePolynomial

model = sys.argv[1].split("/")[-1].strip(".ode")
system = parser.read_system(sys.argv[1])
f = open(sys.argv[2])
lines = f.readlines()
obs_sets = [eval(line) for line in lines]
f.close()


if len(obs_sets) == 1:
    for obs_set in obs_sets:
        print("===============================================")
        obs_polys = [SparsePolynomial.from_string(s, system['variables']) for s in obs_set]


        lumped = clue.do_lumping(system['equations'], obs_polys)
        subspaces = lumped["subspace"]
        polynomial = lumped["polynomials"]
        with open(model + '_CLUE_result.txt','w') as afile:
            for lines in subspaces:
                for i in lines:
                    print(i,end=" ",file=afile)
                print(file=afile)
        with open(model + '_CLUE_ode.txt','w') as afile:
            for lines in polynomial:
                print(lines,file = afile)
    os.chdir("src/")
    os.system("julia LumpPositive.jl " + str(model)) 
    print("Computed new matrix and ODE \n")
    os.chdir("../")


    os.remove(model + '_CLUE_ode.txt')
    f = open(model + '_final_macrovariables.txt')
    lines = f.readlines()
    f.close()
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
    
else:
    obs = 0
    for obs_set in obs_sets:
        obs += 1
        print("===============================================")
        obs_polys = [SparsePolynomial.from_string(s, system['variables']) for s in obs_set]

        if not os.path.exists(model + '_obs' + str(obs)):
            os.makedirs(model + '_obs' + str(obs))

        lumped = clue.do_lumping(system['equations'], obs_polys)
        subspaces = lumped["subspace"]
        polynomial = lumped["polynomials"]
        with open(model + '_obs' + str(obs) + '/' +  model + '_CLUE_result.txt','w') as afile:
            for lines in subspaces:
                for i in lines:
                    print(i,end=" ",file=afile)
                print(file=afile)
        with open(model + '_obs' + str(obs) + '/' +  model + '_CLUE_ode.txt','w') as afile:
            for lines in polynomial:
                print(lines,file = afile)

        os.chdir("src/")
        os.system("julia LumpPositive.jl " + str(model) + " obs" + str(obs))  
        print("Computed new matrix and ODE \n")
        os.chdir("../")


        os.remove(model + '_obs' + str(obs) + '/' + model + '_CLUE_ode.txt')
        f = open(model + '_obs' + str(obs) + '/' +  model + '_final_macrovariables.txt')
        lines = f.readlines()
        f.close()
        with open(model + '_obs' + str(obs) + '/' +  model + '_final_macrovariables.txt', 'w') as afile:
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

        