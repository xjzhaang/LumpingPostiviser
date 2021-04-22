import sys
import os
import time

from sympy import QQ

sys.path.insert(0, "CLUE/")
#sys.path.insert(0, "CLUE/examples/ProteinPhosphorylation/")

import parser
import clue
from sparse_polynomial import SparsePolynomial

#PROTEIN PHOSPHORYLATION CASE
if sys.argv[1] == 'm2.ode' or sys.argv[1] == 'm3.ode' or sys.argv[1] =='m4.ode' or sys.argv[1] == 'm5.ode':
    n = int(sys.argv[1][1])
    system = parser.read_system('examples/ProteinPhosphorylation/'f'm{n}/'f'm{n}.ode')

    obs = [
        clue.SparsePolynomial.from_string("S0", system["variables"]),
        clue.SparsePolynomial.from_string("S1", system["variables"])
    ]
    print(system["variables"])
    lumped = clue.do_lumping(system['equations'], obs)
    subspaces = lumped["subspace"]
    polynomial = lumped["polynomials"]
    with open('examples/ProteinPhosphorylation/m'+str(n)+'/CLUE_files/CLUE_matrix.txt','w') as afile:
        for lines in subspaces:
            for i in lines:
                print(i,end=" ",file=afile)
            print(file=afile)
    with open('examples/ProteinPhosphorylation/m'+str(n)+'/CLUE_files/CLUE_ode.txt','w') as afile:
        for lines in polynomial:
            print(lines,file = afile)
    #EXECUTE JULIA
    os.chdir("src/")
    os.system("julia LumpPositive.jl PP " + str(n))
    print("Computed new matrix and ODE \n")
    #Make matrix into macrovariables
    os.chdir("../")
    f = open('examples/ProteinPhosphorylation/'f'm{n}/Output/macrovariables.txt')
    lines = f.readlines()
    f.close()
    with open('examples/ProteinPhosphorylation/'f'm{n}/Output/macrovariables.txt', 'w') as afile:
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

#ANY OTHER MODEL
else:
    model = sys.argv[1].strip(".ode")

    if os.path.exists('examples/' + model + '/' + sys.argv[1]):
        system = parser.read_system('examples/' + model + '/' + sys.argv[1])
    elif os.path.exists('examples/' + sys.argv[1]):
        system = parser.read_system('examples/' + sys.argv[1])
        os.makedirs('examples/' + model + '/')
        os.rename('examples/' + sys.argv[1], 'examples/' + model + '/' + sys.argv[1])
    
    
    if os.path.exists('examples/' + model + '/' + sys.argv[2]):
        f = open('examples/' + model + '/' + sys.argv[2])
        lines = f.readlines()
        obs_sets = [eval(line) for line in lines]
        f.close()

    elif os.path.exists('examples/' + sys.argv[2]):
        f = open('examples/' + sys.argv[2])
        lines = f.readlines()
        obs_sets = [eval(line) for line in lines]
        f.close()
        os.rename('examples/' + sys.argv[2], 'examples/' + model + '/' + sys.argv[2])

    
    if len(obs_sets) == 1:
        for obs_set in obs_sets:
            print("===============================================")
            obs_polys = [SparsePolynomial.from_string(s, system['variables']) for s in obs_set]


            lumped = clue.do_lumping(system['equations'], obs_polys)
            subspaces = lumped["subspace"]
            polynomial = lumped["polynomials"]
            if not os.path.exists('examples/' + model + '/CLUE_files/'):
                os.makedirs('examples/' + model + '/CLUE_files/')
            if not os.path.exists('examples/' + model + '/Output/'):
                os.makedirs('examples/' + model + '/Output/')
            with open('examples/' + model + '/CLUE_files/CLUE_matrix.txt','w') as afile:
                for lines in subspaces:
                    for i in lines:
                        print(i,end=" ",file=afile)
                    print(file=afile)
            with open('examples/' + model + '/CLUE_files/CLUE_ode.txt','w') as afile:
                for lines in polynomial:
                    print(lines,file = afile)
        os.chdir("src/")
        os.system("julia LumpPositive.jl " + str(model)) 
        print("Computed new matrix and ODE \n")


        os.chdir("../")
        f = open('examples/' + model + '/Output/macrovariables.txt')
        lines = f.readlines()
        f.close()
        with open('examples/' + model + '/Output/macrovariables.txt', 'w') as afile:
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

            if not os.path.exists('examples/' + model + '/obs' + str(obs) + '/CLUE_files/'):
                os.makedirs('examples/' + model + '/obs' + str(obs) + '/CLUE_files/')
            if not os.path.exists('examples/' + model + '/obs' + str(obs) + '/Output/'):
                os.makedirs('examples/' + model + '/obs' + str(obs) + '/Output/')

            lumped = clue.do_lumping(system['equations'], obs_polys)
            subspaces = lumped["subspace"]
            polynomial = lumped["polynomials"]
            with open('examples/' + model + '/obs' + str(obs) + '/CLUE_files/CLUE_matrix.txt','w') as afile:
                for lines in subspaces:
                    for i in lines:
                        print(i,end=" ",file=afile)
                    print(file=afile)
            with open('examples/' + model + '/obs' + str(obs) + '/CLUE_files/CLUE_ode.txt','w') as afile:
                for lines in polynomial:
                    print(lines,file = afile)
    
        os.chdir("src/")
        os.system("julia LumpPositive.jl " + str(model) + " obs" + str(obs))  
        print("Computed new matrix and ODE \n")


        os.chdir("../")
        f = open('examples/' + model + '/obs' + str(obs) +'/Output/macrovariables.txt')
        lines = f.readlines()
        f.close()
        with open('examples/' + model + '/obs' + str(obs) +'/Output/macrovariables.txt', 'w') as afile:
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

        