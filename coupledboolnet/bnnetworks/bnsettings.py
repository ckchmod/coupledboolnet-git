"""
bnsettings
"""

import os
import numpy as np


class ImportSavedData:
    path = os.getcwd()
    networkpath = path + "/data/example_diversifying2/"
    #networkpath = path + "/data/example1/"
    #importsaveddata = False
    importsaveddata = True

    if importsaveddata is True:
        ttable = np.loadtxt(networkpath + "tt.txt", dtype=int)
        initstate = np.loadtxt(networkpath + "init.txt", dtype=bool)
        varf = np.loadtxt(networkpath + "varf.txt", dtype=int)
        nv = np.loadtxt(networkpath + "nv.txt", dtype=int)
    else:
        ttable = None
        initstate = None
        varf = None
        nv = None


class RBNVariables:
    """
    This gets overwritten, if importsaveddata is different (True)
    """
    n = 10
    k = 2
    p = .55


class PerturbationInputVariabale:
    perturbtype = "stochastic"

    assert (perturbtype in ["probabilistic", "periodic", "stochastic", "none"]), \
        "perturbation type needs to be one of probabilistic, periodic, stochastic, none"

    if (perturbtype == "probabilistic"):
        booleanperturb = 0.01
        commnodeprob = 0.0

        if (commnodeprob > 0):
            defaultnode = 0
            if (booleanperturb > 0):
                perturbperturbtype = "probabilistic_and_single_perturb"
            else:
                perturbperturbtype = "singleperturb_and_no_prob"
        else:
            perturbperturbtype = "probabilistic_and_no_perturb"

    elif (perturbtype == "periodic"):
        period = 0
        defaultnode = 0

    elif (perturbtype == "stochastic"):
        stochastictype = "Ising"
        booleanperturb = .01
        defaultnode = 0
        outputnode = 0

    elif(perturbtype == "none"):
        booleanperturb = 0
        commnodeprob = 0

class GridVariables:
    numcells = 10**2
    dT = 10 # Timestep delta for visualizations
    kldthreshold = 0

    if (numcells > 1):

        grid = True
        periodicboundary = True
        comm = "ising"

        assert (comm in ["ising"]), \
            "communication type needs to be one of Ising, ..."

        if (comm == "ising"):
            """
            PARAMETERS: 
            ___________
            J : the exchange energy between the spins 
            h : The magnetic moment is given by µ (external magnetic field)
            T_c : temperature (J/k_B)
            """

            J = -1 # Antiferromagnetic
            h = .001
            T_c = .001 # 2.269 is critical   `
            kb = 8.617 * (10 ** -5)  # eV/K

            initconfig = "antiferromagnetic"
            assert (initconfig in ["antiferromagnetic", "ferromagnetic", "disordered"])

        elif (comm == "linearthreshold"):
            linearthreshold = 2
        else:
            raise Exception("communication type must be one of ising, linearthreshold")

