import numpy as np
import timeit
import os
from coupledboolnet.bnnetworks.bn import BooleanNetwork, inputtest
from coupledboolnet.bnnetworks.ising import Ising
from coupledboolnet.bnnetworks.bnsettings import RBNVariables, PerturbationInputVariabale, GridVariables, ImportSavedData
from coupledboolnet.bnnetworks.bnstatistics import steadystatesrobust
from coupledboolnet.visualization.steadystates import *
from coupledboolnet.bnnetworks.bnstatistics import kldpairwise, lyapunovexp
import pandas as pd
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
datasave = True

WORKING_PATH = os.getcwd()
SAVE_PATH = WORKING_PATH + "/data/output_dat"

def runsteadystateeverything():

    indexcount = 0
    simcount = 1

    for i in range(simcount):
        # inputtest(rbnobj,importeddata, perturbobj, gridobj, timestep)
        rbnP = BooleanNetwork(rbnobj, perturbobj, gridobj, importeddata)
        rbnP.bool_next_all(timestep, gridobj)

        steady, bins = steadystatesrobust(rbnP.states)

        kldmatrix = kldpairwise(steady)

        if datasave is True:

            diversity_threshold = np.mean(kldmatrix)

            dataCol = ["indexCount", "simcount", "k", "p","Lyapunov", "meanKLD", "medianKLD",
                       "varKLD", "t_final", "J", "h", "T_c"]
            dataRow = np.zeros((len(dataCol),1))
            lastRow = dataRow.shape[0]
            dataRow[lastRow, 0] = indexcount
            dataRow[lastRow, 1] = simcount
            dataRow[lastRow, 2] = rbnP.k
            dataRow[lastRow, 3] = rbnP.p
            dataRow[lastRow, 4] = lyapunovexp(rbnP.k,rbnP.p)
            dataRow[lastRow, 5] = diversity_threshold
            dataRow[lastRow, 6] = np.median(kldmatrix)
            dataRow[lastRow, 7] = np.var(kldmatrix)
            dataRow[lastRow, 8] = timestep
            dataRow[lastRow, 9] = gridobj.J
            dataRow[lastRow, 10] = gridobj.h
            dataRow[lastRow, 11] = gridobj.T_c

            if (diversity_threshold > 1):
                objname = str(indexcount) + "_" + str(simcount) + "_" + str(rbnP.k) + "_" + str(rbnP.p) + "_" + \
                           str(gridobj.h) + "_" + str(gridobj.T_c)+ ".npy"
                dataname = "data.npy"
                objfilename = os.path.join(SAVE_PATH, objname)
                datarowname = os.path.join(SAVE_PATH, dataname)
                np.save(objfilename, rbnP)
                np.save(datarowname, dataRow)

    #####     Visualizations
    # showanimation(rbnP.grid.numcells , rbnP.states, rbnP.n, perturbobj.defaultnode , gridobj)
    # statedistributionviz(rbnP.grid.numcells, rbnP.states, rbnP.n, perturbobj.booleanperturb, gridobj)
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
    print(timeit.timeit(stmt=runsteadystateeverything, number=2))

def main():
    runsteadystateeverything()
    #checktime()
    #isingcheck()
    print("Simulation Finished.")

main()













