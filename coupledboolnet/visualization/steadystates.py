"""
Steady-state Distribution
"""

from coupledboolnet.bnnetworks.bn import bitstoints, steadystates
import pandas as pd
import matplotlib.pyplot as plt


def stateevolutionviz(states):
    print("\n --- Time evolution of States --- ")  
    
    states = bitstoints(states)
    stateevolution = pd.DataFrame(states)
    stateevolution.plot(title="Time Evolution of BN")
    plt.xlabel('Time t')
    plt.ylabel('Decimal Encoding of Binary States')
    
def statedistributionviz(states, numgenes, stringinfo):
    print("\n --- Steady-state Distribution of States --- ")  
    
    ssdistribution, binnum = steadystates(states,numgenes)
    plt.bar(binnum, ssdistribution)
    plt.title('p = ' + "{:.2f}".format(stringinfo), fontsize=9)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) 
    
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False,
        labelleft=False) 
    plt.title('Steady-state Distribution of BN: p' + str(stringinfo), fontsize=9)
    plt.xlabel('Decimal Encoding of Binary States', fontsize=7)
    plt.ylabel('Counts', fontsize=7)
