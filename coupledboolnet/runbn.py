import matplotlib.pyplot as plt
import numpy as np
from math import ceil, sqrt
from numpy import save
from coupledboolnet.bnnetworks.bn import BooleanNetwork, inputtest, bitstoints
from coupledboolnet.bnnetworks.ising import Ising
from coupledboolnet.bnnetworks.bnsettings import RBNVariables, PerturbationInputVariabale, GridVariables, ImportSavedData
from coupledboolnet.visualization.steadystates import *
"""
Run RBNp
@author: chrisk
"""
rbnobj = RBNVariables()
perturbobj = PerturbationInputVariabale()
gridobj = GridVariables()
importeddata = ImportSavedData()

# Simulation Parameters
numiter = 1
timestep = 2**10

def runsteadystateeverything():
    #inputtest(rbnobj,importeddata, perturbobj, gridobj, timestep)
    rbnP = BooleanNetwork(rbnobj, perturbobj, gridobj, importeddata)
    rbnP.bool_next_all(timestep, gridobj)
    # save('data', rbnP.states)
    #statedistributionviz(rbnP.grid.numcells , rbnP.states, rbnP.n, perturbobj.booleanperturb)
    showanimation(rbnP.grid.numcells , rbnP.states, rbnP.n, perturbobj.defaultnode , gridobj.dT)
    # transitiondiagram(rbnP.state_transition())

def isingcheck():
    """
    TBD
    """
    isingmodel = Ising(GridVariables, timestep)
    print(isingmodel.isingrunall())
    pass

def checktime():
    """
    Just benchmark for running code
    """
    import timeit
    print(timeit.timeit(stmt=runsteadystateeverything, number=10))

def main():
    runsteadystateeverything()
    #checktime()
    #isingcheck()
    print("Simulation Finished.")

main()













