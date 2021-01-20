from coupledboolnet.bnnetworks.bn import BooleanNetwork
import os
from numpy import arange
import pandas as pd
import matplotlib.pyplot as plt
import pickle as pkl
from coupledboolnet.visualization.steadystates import *
from coupledboolnet.bnnetworks.bnsettings import PerturbationInputVariabale, GridVariables
from scipy.optimize import curve_fit
import statsmodels.api as sm
def main():

    WORKING_PATH = os.getcwd()
    SAVE_PATH = WORKING_PATH + "/data/output_data/10_by_10_1_3"

    data_col = ["indexcount", "simcount", "k", "p", "Lyapunov", "meanKLD", "medianKLD",
                "varKLD", "t_final", "J", "T_c", "h"]
    df = pd.DataFrame(columns=data_col)
    for i in range(1, 6, 1):
        datafilename = "data" + str(i) + ".pkl"
        DATANAME = os.path.join(SAVE_PATH, datafilename)
        df = df.append(pd.read_pickle(DATANAME))

    print(df)


    #data_col = ["indexcount", "simcount", "k", "p", "Lyapunov", "meanKLD", "medianKLD",
    #            "varKLD", "t_final", "J", "T_c", "h"]
    #df1 = pd.DataFrame(columns=data_col)
    #df2 = pd.DataFrame(columns=data_col)


main()
