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

def viz_animation(DIVERSIFYING_PATH):
    ANI_PATH = DIVERSIFYING_PATH + "/9_1_2_0.55_0.1_1.1.pkl"
    file = open(ANI_PATH, 'rb')
    rbnP = pkl.load(file)
    file.close()
    perturbobj = PerturbationInputVariabale()
    gridobj = GridVariables()
    gridobj.h = .1
    gridobj.T_c = 1.1
    showanimation(rbnP.grid.numcells , rbnP.states, rbnP.n, perturbobj.defaultnode, gridobj)
    statedistributionviz(rbnP.grid.numcells, rbnP.states, rbnP.n, perturbobj.booleanperturb, gridobj)

def power_law(x, a, b):
    return a * np.power(x, b)

def vizsims10by10(DIVERSIFYING_PATH, CUM_PATH):
    
    data_col = ["indexcount", "simcount", "k", "p", "Lyapunov", "meanKLD", "medianKLD",
                "varKLD", "t_final", "J", "h", "T_c"]
    df1 = pd.DataFrame(columns=data_col)
    df2 = pd.DataFrame(columns=data_col)

    for i in range(1, 6, 1):
        datafilename = "data" + str(i) + ".pkl"
        DATANAME = os.path.join(DIVERSIFYING_PATH, datafilename)
        df1 = df1.append(pd.read_pickle(DATANAME))

    for i in range(1, 6, 1):
        datafilename = "data" + str(i) + ".pkl"
        DATANAME = os.path.join(CUM_PATH, datafilename)
        df2 = df2.append(pd.read_pickle(DATANAME))

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    fig.suptitle("10 by 10 Tissue",  fontsize=20)

    df1.plot(kind="scatter", x="T_c", y="meanKLD", alpha=0.25, ax=ax1)
    ax1.set_title('Single Diversifying Network')
    lowess = sm.nonparametric.lowess
    lowess = lowess(df1.meanKLD, df1.T_c, frac=.3)
    lowess_x = list(zip(*lowess))[0]
    lowess_y = list(zip(*lowess))[1]
    ax1.plot(lowess_x, lowess_y, '-', color='red')

    df2.plot(kind="scatter", x="T_c", y="meanKLD", alpha=0.25, ax=ax2)
    ax2.set_title('Multiple Randomized Networks, k=2~3, p=.55~69')
    lowess = sm.nonparametric.lowess
    lowess = lowess(df2.meanKLD, df2.T_c, frac=.3)
    lowess_x_T = list(zip(*lowess))[0]
    lowess_y_T = list(zip(*lowess))[1]
    ax2.plot(lowess_x_T, lowess_y_T, '-', color='red')

    df1.plot(kind="scatter", x="h", y="meanKLD", alpha=0.25, ax=ax3)
    ax3.set_title('Single Diversifying Network')
    # plt.plot(x_line, y_line, '--', color='green')
    lowess = sm.nonparametric.lowess
    lowess = lowess(df1.meanKLD, df1.h, frac=.3)
    lowess_x = list(zip(*lowess))[0]
    lowess_y = list(zip(*lowess))[1]
    ax3.plot(lowess_x, lowess_y, '-', color='red')

    df2.plot(kind="scatter", x="h", y="meanKLD", alpha=0.25, ax=ax4)
    ax4.set_title('Multiple Randomized Networks, k=2~3, p=.55~69')
    lowess = sm.nonparametric.lowess
    lowess = lowess(df2.meanKLD, df2.h, frac=.3)
    lowess_x_h = list(zip(*lowess))[0]
    lowess_y_h = list(zip(*lowess))[1]
    ax4.plot(lowess_x_h, lowess_y_h, '-', color='red')

    # 3D Figure
    fig = plt.figure()
    fig.suptitle("10 by 10 Tissue",  fontsize=20)
    ax1 = fig.add_subplot(121, projection='3d')
    ax1.scatter(df1.T_c, df1.h, df1.meanKLD, c='r', marker='o')
    ax1.set_title('Single Diversifying Network'), ax1.set_xlabel('Temperature'), ax1.set_ylabel(
        'Interaction h'), ax1.set_zlabel('Diversity')

    ax2 = fig.add_subplot(122, projection='3d')
    ax2.scatter(df2.T_c, df2.h, df2.meanKLD, c='r', marker='o')
    ax2.set_title("Multiple Randomized Network"), ax2.set_xlabel('Temperature'), ax2.set_ylabel(
        'Interaction h'), ax2.set_zlabel('Diversity')

    max_value_y_T = max(lowess_y_T)
    max_index = lowess_y_T.index(max_value_y_T)
    max_value_x_T = lowess_x_T[max_index]

    max_value_y_h = max(lowess_y_h)
    max_index = lowess_y_h.index(max_value_y_h)
    max_value_x_h = lowess_x_h[max_index]

    return([max_value_x_T, max_value_y_T], [max_value_x_h, max_value_y_h])

def vizsimstissues(TISSUE_PATH, num):
    data_col = ["indexcount", "simcount", "k", "p", "Lyapunov", "meanKLD", "medianKLD",
                "varKLD", "t_final", "J", "h", "T_c"]
    df1 = pd.DataFrame(columns=data_col)

    for i in range(1, 6, 1):
        datafilename = "data" + str(i) + ".pkl"
        DATANAME = os.path.join(TISSUE_PATH, datafilename)
        df1 = df1.append(pd.read_pickle(DATANAME))

    fig = plt.figure()
    ax1 = plt.subplot(1, 3, 1)
    plt.suptitle(str(num) + " by " + str(num) + " Tissue", fontsize=20)

    df1.plot(kind="scatter", x="T_c", y="meanKLD", alpha=0.25, ax=ax1)
    ax1.set_title('Single Diversifying Network')
    lowess = sm.nonparametric.lowess
    lowess = lowess(df1.meanKLD, df1.T_c, frac=.3)
    lowess_x_T = list(zip(*lowess))[0]
    lowess_y_T = list(zip(*lowess))[1]
    ax1.plot(lowess_x_T, lowess_y_T, '-', color='red')

    ax2 = plt.subplot(1, 3, 2)
    df1.plot(kind="scatter", x="h", y="meanKLD", alpha=0.25, ax=ax2)
    ax2.set_title('Single Diversifying Network')
    # plt.plot(x_line, y_line, '--', color='green')
    lowess = sm.nonparametric.lowess
    lowess = lowess(df1.meanKLD, df1.h, frac=.3)
    lowess_x_h = list(zip(*lowess))[0]
    lowess_y_h = list(zip(*lowess))[1]
    ax2.plot(lowess_x_h, lowess_y_h, '-', color='red')

    ax3 = fig.add_subplot(1, 3, 3, projection='3d')
    ax3.scatter(df1.T_c, df1.h, df1.meanKLD, c='r', marker='o')
    ax3.set_title('Single Diversifying Network'), ax3.set_xlabel('Temperature'), ax3.set_ylabel(
        'Interaction h'), ax3.set_zlabel('Diversity')

    max_value_y_T = max(lowess_y_T)
    max_index = lowess_y_T.index(max_value_y_T)
    max_value_x_T = lowess_x_T[max_index]

    max_value_y_h = max(lowess_y_h)
    max_index = lowess_y_h.index(max_value_y_h)
    max_value_x_h = lowess_x_h[max_index]

    return([max_value_x_T, max_value_y_T], [max_value_x_h, max_value_y_h])


def main():
    WORKING_PATH = os.getcwd()
    SAVE_PATH = WORKING_PATH + "/data/output_data"
    DIVERSIFYING_PATH = SAVE_PATH + "/example_diversifying_kamiak_ising"
    CUM_PATH = SAVE_PATH + "/example_cum_kamiak_ising"
    TISSUE_PATH02 = SAVE_PATH + "/multi_tissue_02"
    TISSUE_PATH04 = SAVE_PATH + "/multi_tissue_04"
    TISSUE_PATH08 = SAVE_PATH + "/multi_tissue_08"
    #TISSUE_PATH16 = SAVE_PATH + "/multi_tissue_16"

    T10, h10 = vizsims10by10(DIVERSIFYING_PATH, CUM_PATH)
    T02, h02 = vizsimstissues(TISSUE_PATH02, 2)
    T04, h04 = vizsimstissues(TISSUE_PATH04, 4)
    T08, h08 = vizsimstissues(TISSUE_PATH08, 8)
    # T16, h16 = vizsimstissues(TISSUE_PATH16, 16)

    gridsizes =np.array([4,16,64,100])
    Tdat = np.array([T02, T04, T08, T10], dtype=float)
    hdat =  np.array([h02, h04, h08, h10], dtype=float)


    print(Tdat)
    print(hdat)

    fig, (ax1, ax2) = plt.subplots(1, 2)
    plt.suptitle("Mutliple Tissues Power Law?", fontsize=20)

    ax1.scatter(gridsizes, Tdat[:,0], marker='o')
    pars, cov = curve_fit(f=power_law, xdata=gridsizes, ydata=Tdat[:,0], p0=[0, 0], bounds=(-np.inf, np.inf))
    x_line = arange(min(gridsizes), max(gridsizes),1)
    a, b = pars
    y_line = power_law(x_line, a, b)
    ax1.plot(x_line, y_line, '--', color='red')
    ax1.set_xlabel("Tissue Size"), ax1.set_ylabel("Critical Temperature")

    ax2.scatter(gridsizes, hdat[:,0], marker='o')
    pars, cov = curve_fit(f=power_law, xdata=gridsizes, ydata=hdat[:,0], p0=[0, 0], bounds=(-np.inf, np.inf))
    x_line = arange(min(gridsizes), max(gridsizes),1)
    a, b = pars
    y_line = power_law(x_line, a, b)
    ax2.plot(x_line, y_line, '--', color='red')
    ax2.set_xlabel("Tissue Size"), ax2.set_ylabel("Critical h")

    plt.show()

main()


# Question: What kind of analyses can be done on top of this?
# For paper
# 0) prbabilstic canalyzing function
# 1) find critical h_c and T_c (interporlation) and number of cells
# order parameter should be proportional to
# universitality class.

# 4x4, 8x8 Critical numbers (as a function of number of cells)

# Future
# with critical parametrs, develop Ladelle model(?)
