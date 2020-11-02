from coupledboolnet.bnnetworks.bn import BooleanNetwork
import os
import pandas as pd
import matplotlib.pyplot as plt
import pickle as pkl
from coupledboolnet.visualization.steadystates import *
from coupledboolnet.bnnetworks.bnsettings import PerturbationInputVariabale, GridVariables

WORKING_PATH = os.getcwd()
SAVE_PATH = WORKING_PATH + "/data/output_data"
DIVERSIFYING_PATH = SAVE_PATH + "/example_diversifying_kamiak_ising"
CUM_PATH = SAVE_PATH + "/example_cum_kamiak_ising"
ANI_PATH = DIVERSIFYING_PATH + "/9_1_2_0.55_0.1_1.1.pkl"

file = open(ANI_PATH,'rb')
rbnP = pkl.load(file)
file.close()
perturbobj = PerturbationInputVariabale()
gridobj = GridVariables()
gridobj.h = .1
gridobj.T_c = 1.1

#statedistributionviz(rbnP.grid.numcells, rbnP.states, rbnP.n, perturbobj.booleanperturb, gridobj)

data_col = ["indexcount", "simcount", "k", "p", "Lyapunov", "meanKLD", "medianKLD",
               "varKLD", "t_final", "J", "h", "T_c"]
df1 = pd.DataFrame(columns = data_col)
df2 = pd.DataFrame(columns = data_col)

for i in range(1, 6, 1):
    datafilename = "data" + str(i) + ".pkl"
    DATANAME = os.path.join(DIVERSIFYING_PATH, datafilename)
    df1 = df1.append(pd.read_pickle(DATANAME))

for i in range(1, 6, 1):
    datafilename = "data" + str(i) + ".pkl"
    DATANAME = os.path.join(CUM_PATH, datafilename)
    df2 = df2.append(pd.read_pickle(DATANAME))

#print(df1)
df1.plot(kind="scatter", x="h", y="meanKLD")
plt.title('Single Diversifying Network')
df2.plot(kind="scatter", x="h", y="meanKLD")
plt.title('Multiple Randomized Networks, k=2~3, p=.55~69')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(df1.T_c, df1.h, df1.meanKLD, c='r', marker='o')
ax.set_title('Single Diversifying Network'), ax.set_xlabel('Temperature'), ax.set_ylabel('Interaction h'), ax.set_zlabel('Diversity')
fig2 = plt.figure()
ax = fig2.add_subplot(111, projection='3d')
ax.scatter(df2.T_c, df2.h, df2.meanKLD, c='r', marker='o')
ax.set_title("Multiple Randomized Network"), ax.set_xlabel('Temperature'), ax.set_ylabel('Interaction h'), ax.set_zlabel('Diversity')
plt.show()

#showanimation(rbnP.grid.numcells , rbnP.states, rbnP.n, perturbobj.defaultnode, gridobj)

# Question: What kind of analyses can be done on top of this?
# For paper
# 0) prbabilstic canalyzing function
# 1) find critical h_c and T_c (interporlation) and number of cells
# order parameter should be proportional to
# universitality class.

# 4x4, 8x8 Critical numbers (as a function of number of cells)

# Future
# with critical parametrs, develop Ladelle model(?)