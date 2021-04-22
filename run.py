import sys
import os
import time

from sympy import QQ

sys.path.insert(0, "CLUE/")
#sys.path.insert(0, "CLUE/examples/ProteinPhosphorylation/")

import parser
import clue
from sparse_polynomial import SparsePolynomial

if sys.argv[1] == "PP":
    n = int(sys.argv[2])
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
    os.system("julia LumpPositive.jl PP " + sys.argv[2])
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

elif sys.argv[1] == "Jak":
    system = parser.read_system('examples/Jak-family/Jak.ode')
    obs_sets = [
        ["S0"],
    ]
    for obs_set in obs_sets:
        print("===============================================")
        obs_polys = [SparsePolynomial.from_string(s, system['variables']) for s in obs_set]
        lumped = clue.do_lumping(system['equations'], obs_polys)
        subspaces = lumped["subspace"]
        polynomial = lumped["polynomials"]
        with open('examples/Jak-family/CLUE_files/CLUE_matrix.txt','w') as afile:
            for lines in subspaces:
                for i in lines:
                    print(i,end=" ",file=afile)
                print(file=afile)
        with open('examples/Jak-family/CLUE_files/CLUE_ode.txt','w') as afile:
            for lines in polynomial:
                print(lines,file = afile)

    os.chdir("src/")
    os.system("julia LumpPositive.jl Jak") 

    print("Computed new matrix and ODE \n")
    os.chdir("../")
    f = open('examples/Jak-family/Output/macrovariables.txt')
    lines = f.readlines()
    f.close()
    with open('examples/Jak-family/Output/macrovariables.txt', 'w') as afile:
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


elif sys.argv[1] == "fceri_ji":
    system = parser.read_system('examples/fceri_ji/fceri_ji.ode')

    obs_sets = [
        #["S0"],
        #["S2", "S178", "S267", "S77"],
        ["S2 + S178 + S267 + S77"],
        #["S7"],
        #["S1"]
    ]

    for obs_set in obs_sets:
        print("===============================================")
        obs_polys = [SparsePolynomial.from_string(s, system['variables']) for s in obs_set]


        lumped = clue.do_lumping(system['equations'], obs_polys)
        subspaces = lumped["subspace"]
        polynomial = lumped["polynomials"]
        with open('examples/fceri_ji/CLUE_files/CLUE_matrix.txt','w') as afile:
            for lines in subspaces:
                for i in lines:
                    print(i,end=" ",file=afile)
                print(file=afile)
        with open('examples/fceri_ji/CLUE_files/CLUE_ode.txt','w') as afile:
            for lines in polynomial:
                print(lines,file = afile)
    
    os.chdir("src/")
    os.system("julia LumpPositive.jl fceri_ji") 
    print("Computed new matrix and ODE \n")


    os.chdir("../")
    f = open('examples/fceri_ji/Output/macrovariables.txt')
    lines = f.readlines()
    f.close()
    with open('examples/fceri_ji/Output/macrovariables.txt', 'w') as afile:
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
    print("Bad argument. Possible arguments are: PP (number), fceri_ji, Jak")