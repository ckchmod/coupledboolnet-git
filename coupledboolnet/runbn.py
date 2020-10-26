import timeit, os
import pandas as pd, numpy as np, matplotlib.pyplot as plt, pickle as pkl
from coupledboolnet.bnnetworks.bn import BooleanNetwork, inputtest
from coupledboolnet.bnnetworks.ising import Ising
from coupledboolnet.bnnetworks.bnsettings import RBNVariables, PerturbationInputVariabale, GridVariables, ImportSavedData
from coupledboolnet.bnnetworks.bnstatistics import steadystatesrobust
from coupledboolnet.visualization.steadystates import *
from coupledboolnet.bnnetworks.bnstatistics import kldpairwise, lyapunovexp
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
SAVE_PATH = WORKING_PATH + "/data/output_data"

def runsteadystateeverything():

    indexcount = 0
    simcount = 1

    for i in range(simcount):
        # inputtest(rbnobj,importeddata, perturbobj, gridobj, timestep)
        rbnP = BooleanNetwork(rbnobj, perturbobj, gridobj, importeddata)
        rbnP.bool_next_all(timestep, gridobj)

        steady, bins = steadystatesrobust(rbnP.states)

        kldmatrix = kldpairwise(steady)
        #viewkldbar(kldmatrix)

        # need to load up df here
        if datasave is True:

            mean_KLD = np.mean(kldmatrix)
            dataCol = ["indexcount", "simcount", "k", "p","Lyapunov", "meanKLD", "medianKLD",
                       "varKLD", "t_final", "J", "h", "T_c"] #Turn this into pandas dataframe
            dfsave = pd.DataFrame(columns=dataCol)
            rowdf = pd.DataFrame([[indexcount, simcount, rbnP.k, rbnP.p, lyapunovexp(rbnP.k, rbnP.p), mean_KLD,
                       np.median(kldmatrix), np.var(kldmatrix), timestep, gridobj.J, gridobj.h, gridobj.T_c]],
                                 columns=dataCol)


            dfsave = dfsave.append(rowdf, ignore_index= True)
            dataname = "data.pkl"
            datarowname = os.path.join(SAVE_PATH, dataname)
            dfsave.to_pickle(datarowname)

            if (mean_KLD > gridobj.kldthreshold):
                objname = str(indexcount) + "_" + str(simcount) + "_" + str(rbnP.k) + "_" + str(rbnP.p) + "_" + \
                           str(gridobj.h) + "_" + str(gridobj.T_c)+ ".pkl"
                objfilename = os.path.join(SAVE_PATH, objname)
                pkl.dump(rbnP, open(objfilename, "wb"))

    #####     Visualizations
    # showanimation(rbnP.grid.numcells , rbnP.states, rbnP.n, perturbobj.defaultnode , gridobj)
    # statedistributionviz(rbnP.grid.numcells, rbnP.states, rbnP.n, perturbobj.booleanperturb, gridobj)
    # transitiondiagram(rbnP.state_transition())

def isingcheck():
    """
    Check all the variables are correct
    """
    isingmodel = Ising(GridVariables, timestep)
    print(isingmodel.isingrunall())
    pass

def checktime():
    """
    Benchmark time for running code
    """
    print(timeit.timeit(stmt=runsteadystateeverything, number=2))

def main():
    runsteadystateeverything()
    #checktime()
    #isingcheck()
    print("Simulation Finished.")

main()













