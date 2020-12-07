import os, warnings, time
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

def run_multi_sim(sys_arg, rbnobj, perturbobj, gridobj, importeddata, timestep, datasave, showviz):
    """
    Run Multi-Sim
    """
    if datasave is True:
        WORKING_PATH = os.getcwd()
        SAVE_PATH = WORKING_PATH + "/coupledboolnet/data/output_data"
        datafilename = "data" + str(sys_arg) + ".pkl"
        DATANAME = os.path.join(SAVE_PATH, datafilename)

    indexcount = 0
    simcount = 1
    sys_num = int(sys_arg)
    inc = 2
    k_range = np.linspace(2, 3, 2, dtype=int)
    p_range = np.linspace(0.55, 0.69, 2)
    h_range = np.linspace( sys_num - 1 + 0.1, sys_num + 0.1, 10)
    T_c_range = np.linspace( sys_num - 1 + 0.1, sys_num + 0.1, 10)
    data_col = ["indexcount", "simcount", "k", "p", "Lyapunov", "meanKLD", "medianKLD",
               "varKLD", "t_final", "J", "h", "T_c"]
    dfsave = pd.DataFrame(columns=data_col)

    for kk in k_range:
        for pp in p_range:
            for hh in h_range:
                for tc in T_c_range:
                    for i in range(simcount):
                        # inputtest(rbnobj,importeddata, perturbobj, gridobj, timestep)

                        rbnobj.k = kk
                        rbnobj.p = pp
                        gridobj.T_c = tc
                        gridobj.h = hh
                        rbnP = BooleanNetwork(rbnobj, perturbobj, gridobj, importeddata)
                        rbnP.bool_next_all(timestep, gridobj)

                        start_time = time.time()
                        steady, bins = steadystatesrobust(rbnP.states)
                        print("bn execution: --- %s seconds ---" % (time.time() - start_time))

                        kldmatrix = kldpairwise(steady)
                        mean_KLD = np.mean(kldmatrix)

                        rowdf = pd.DataFrame([[indexcount, simcount, rbnP.k, rbnP.p, lyapunovexp(rbnP.k, rbnP.p), mean_KLD,
                                               np.median(kldmatrix), np.var(kldmatrix), timestep, gridobj.J, gridobj.h,
                                               gridobj.T_c]],
                                             columns=data_col)
                        dfsave = dfsave.append(rowdf, ignore_index=True)
                        print("indexcount: " + str(indexcount) + " k: " + str(kk) + " p: " + str(pp) + " h: " + str(hh) + " T_c: " + str(tc) + " mean_KLD: " + str(mean_KLD))

                        # need to load up df here
                        if (datasave is True and (mean_KLD > gridobj.kldthreshold)):
                                objname = "s" + str(sys_arg) + "_" + str(indexcount) + "_" + str(simcount) + "_" + str(rbnP.k) + "_" + str(rbnP.p) + "_" + \
                                           str(round(gridobj.h,2)) + "_" + str(round(gridobj.T_c,2))+ ".pkl"
                                objfilename = os.path.join(SAVE_PATH, objname)
                                pkl.dump(rbnP, open(objfilename, "wb"))
                        indexcount = indexcount + 1

    if datasave is True:
        dfsave.to_pickle(DATANAME)

    if showviz is True:
        showanimation(rbnP.grid.numcells , rbnP.states, rbnP.n, perturbobj.defaultnode , gridobj)
        statedistributionviz(rbnP.grid.numcells, rbnP.states, rbnP.n, perturbobj.booleanperturb, gridobj)
        transitiondiagram(rbnP.state_transition())
        viewkldbar(kldmatrix)

def isingcheck(timestep):
    """
    Check all the variables are correct
    """
    isingmodel = Ising(GridVariables, timestep)
    print(isingmodel.isingrunall())
    pass

def main(sys_arg):

    print("Simulation Started.")

    # Global network parameters
    rbnobj = RBNVariables()
    perturbobj = PerturbationInputVariabale()
    gridobj = GridVariables()
    importeddata = ImportSavedData()

    # Simulation Parameters
    timestep = 5000
    datasave = True
    showviz = False

    warnings.filterwarnings("ignore", category=RuntimeWarning)
    run_multi_sim(sys_arg, rbnobj, perturbobj, gridobj, importeddata, timestep, datasave, showviz)
    #isingcheck(timestep)
    print("Simulation Finished.")

#main()













