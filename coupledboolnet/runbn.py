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

    # List Parameters
    indexcount = 0
    sys_num = int(sys_arg)
    k_range = np.linspace(sys_num, sys_num, 1, dtype=int)
    p_range = np.round(np.linspace(0.55, 0.69, 2), 2)
    T_c_range = np.round(np.linspace(0.1, 5, 2), 2)
    h_range = np.round(np.linspace(0.1, 5, 2), 2)
    simcount =np.linspace(1,1,1, dtype=int)

    # Simulation Begins
    data_col = ["indexcount", "simcount", "k", "p", "Lyapunov", "meanKLD", "medianKLD",
                "varKLD", "t_final", "J", "T_c", "h"]
    dfsave = pd.DataFrame(columns=data_col)

    # Check to see if the items are in
    path_exists = os.path.exists(DATANAME)

    if (path_exists is True):
        file = open(DATANAME, 'rb')
        dfsave = pkl.load(file)
        file.close()

        kk_index = np.where(k_range == dfsave["k"].iloc[-1])[0][0]
        pp_index = np.where(p_range == dfsave["p"].iloc[-1])[0][0]
        tc_index = np.where(T_c_range == dfsave["T_c"].iloc[-1])[0][0]
        hh_index = np.where(h_range == dfsave["h"].iloc[-1])[0][0]
        sc_index = np.where(simcount == dfsave["simcount"].iloc[-1])[0][0]
        indexcount = dfsave["indexcount"].iloc[-1]
        print(" --- Restarting from Checkpoint ---- ")
        print(str(kk_index) + "_" + str(pp_index) + "_" + str(tc_index) + "_"
              + str(hh_index) + "_" + str(sc_index))
    else:
        kk_index = 0
        pp_index = 0
        tc_index = 0
        hh_index = 0
        sc_index = 0
        indexcount = 0

    for kk in k_range[kk_index:]:
        for pp in p_range[pp_index:]:
            for tc in T_c_range[tc_index:]:
                for hh in h_range[hh_index:]:
                    for sc in simcount[sc_index:]:
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
                        lyap = np.round(lyapunovexp(rbnP.k, rbnP.p),5)
                        mean_KLD = np.round(np.mean(kldmatrix), 5)
                        med_KLD = np.round(np.median(kldmatrix), 5)
                        var_kld = np.round(np.var(kldmatrix), 5)

                        rowdf = pd.DataFrame([[indexcount, sc, rbnP.k, rbnP.p, lyap, mean_KLD,
                                               med_KLD, var_kld, timestep, gridobj.J, gridobj.T_c, gridobj.h]],
                                             columns=data_col)

                        dfsave = dfsave.append(rowdf, ignore_index=True)
                        print("ic: " + str(indexcount) + " sc: " + str(sc) + " k: " + str(kk) + " p: " + str(pp)  +
                              " T_c: " + str(tc) + " h: " + str(hh) + " mean_KLD: " + str(mean_KLD))

                        # need to load up df here
                        if (datasave is True and (mean_KLD > gridobj.kldthreshold)):
                            objname = "s" + str(sys_num) + "_" + str(indexcount) + "_" + str(sc) + "_" + str(rbnP.k) + \
                                      "_" + str(rbnP.p) + "_" + str(tc) + "_" + str(hh) + ".pkl"
                            objfilename = os.path.join(SAVE_PATH, objname)
                            pkl.dump(rbnP, open(objfilename, "wb"))

                        indexcount = indexcount + 1

                        if datasave is True:
                            dfsave.to_pickle(DATANAME)
                    sc_index = 0
                hh_index = 0
            tc_index = 0
        pp_index = 0
    kk_index = 0

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
    timestep = 12800
    datasave = True
    showviz = False

    warnings.filterwarnings("ignore", category=RuntimeWarning)
    run_multi_sim(sys_arg, rbnobj, perturbobj, gridobj, importeddata, timestep, datasave, showviz)
    #isingcheck(timestep)
    print("Simulation Finished.")

#main()













