import matplotlib.pyplot as plt
from math import ceil, sqrt
from coupledboolnet.bnnetworks.bn import BooleanNetwork, inputtest
from coupledboolnet.bnnetworks.ising import Ising
from coupledboolnet.bnnetworks.bnsettings import RBNVariables, PerturbationInputVariabale, GridVariables, ImportSavedData
from coupledboolnet.visualization.statetransitiondiagram import transitiondiagram
from coupledboolnet.visualization.steadystates import stateevolutionviz,statedistributionviz
"""
Set all network Parameters

@author: chrisk
"""

rbnobj = RBNVariables()

perturbobj = PerturbationInputVariabale()

gridobj = GridVariables()

importeddata = ImportSavedData()

"""
Simulation Parameters
"""

numiter = 1
timestep = 2**12


def runsteadystateeverything():
    # Check all the inputs are correctinitiall. Incomplete.
    inputtest(rbnobj, importeddata, perturbobj, gridobj, timestep)

    # fig = plt.figure(0)    
    # for i in range(numiter):

    #     # make any incremental paramter change here
    #     perturbobj.booleanperturb = .00 + i*.1
    #     perturbobj.defaultnode = 0

    #     rbnP = BooleanNetwork(rbnobj, perturbobj, gridobj, importeddata)
    #     rbnP.bool_next_all(timestep, gridobj)
        
    #     plt.subplot(ceil(sqrt(numiter)), ceil(sqrt(numiter)), i+1)
    #     statedistributionviz(rbnP.states, rbnP.n, perturbobj.booleanperturb)
        


    rbnP = BooleanNetwork(rbnobj, perturbobj, gridobj, importeddata)
    rbnP.bool_next_all(timestep, gridobj)
    
    # from numpy import save
    # save('./coupledboolnet/data.npy', rbnP.states)

    #matlabsaving
    # import scipy.io
    # statetransitiontomatlab = rbnP.state_transition()
    # scipy.io.savemat('./coupledboolnet/MATLAB_TESTING/test.mat', {'allstates':rbnP.states , 'statetransition': statetransitiontomatlab})
    
    # # Visualize Probabilistic Boolean Network  
   
    # fig = plt.figure(0)
    # plt.subplot(ceil(sqrt(numiter)), ceil(sqrt(numiter)), 1)
    # statedistributionviz(rbnP.states, rbnP.n, perturbobj.booleanperturb)
         
    # fig.text(0.5, 0.99, 'Steady State Distributions: Single Node and Perturbation Node 0', ha='center', va='center')
    # fig.text(0.5, 0.04, 'states', ha='center', va='center')
    # fig.text(0.06, 0.5, 'counts', ha='center', va='center', rotation='vertical')
    # plt.show()
    


    # #Viz
    #fig = plt.figure(1)        
    #transitiondiagram(rbnP.state_transition())
    #plt.show()
    # fig = plt.figure(2)  
    #print(rbnP.states)      
    #stateevolutionviz(rbnP.states)
    # plt.show()

def isingcheck():
 
    isingmodel = Ising(GridVariables, timestep)
    print(isingmodel.isingrun())
    pass

def generatetissue():
    
    pass    
    
def main():
    #runsteadystateeverything()
    import timeit
    print(timeit.timeit(stmt = runsteadystateeverything, number = 100))
    #isingcheck()
    
main()













