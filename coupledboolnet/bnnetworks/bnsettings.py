"""
bnsettings 
"""

import os
import numpy as np

class ImportSavedData:
    
    path = os.getcwd()
    networkpath = path + "/coupledboolnet/data/example1/"
    importsaveddata = True
    
    if importsaveddata is True:
        ttable = np.loadtxt(networkpath +"tt.txt", dtype=int)
        initstate = np.loadtxt(networkpath +"init.txt", dtype=bool)
        varf = np.loadtxt(networkpath +"varf.txt", dtype=int)
        nv = np.loadtxt(networkpath +"nv.txt", dtype=int)
    else:
        ttable = None
        initstate = None
        varf = None
        nv = None

class RBNVariables:
    
    n = 3
    k = 3
    p = .5      

class PerturbationInputVariabale:
    
    perturbtype = "probabilistic"           
    
    if (perturbtype == "probabilistic"):
        booleanperturb = .1
        commnodeprob  = .4

        if (commnodeprob > 0):
            defaultnode = 0
            if (booleanperturb > 0):
                perturbperturbtype = "probabilistic_and_single_perturb"
            else:
                perturbperturbtype = "singleperturb_and_no_prob"
        else:
            perturbperturbtype = "probabilistic_and_no_perturb"

    elif (perturbtype == "periodic"):
        period= 0 
        defaultnode = 0

    elif (perturbtype == "stochastic"): 
        stochastictype = "Brownian"
    else: 
         raise Exception("perturbation type needs to be one of probabilistic, periodic, stochastic")

                
class GridVariables:
    
    numcells = 100
    if (numcells > 1):
        
        grid = True
        periodicboundary = True
        comm = "ising" 
        
        
        if (comm == "ising"):
            """
            
            PARAMETERS: 
            ___________
            J : theexchange energy between the spins 
            h : The magnetic moment is given by µ (external magnetic field)
            T_c : tempreature (J/k_B)
            """
              
            J = 1
            h = 0
            T_c = 4.4 #2.269 is critical
            kb = 8.617 * (10**-5) # eV/K

            initconfig = "antiferromagnetic"
            assert(initconfig in ["antiferromagnetic", "ferromagnetic", "disordered"])
            
        elif (comm == "linearthreshold"):
            linearthreshold = 2
        else:
            raise Exception("communication type must be one of ising, linearthreshold")

