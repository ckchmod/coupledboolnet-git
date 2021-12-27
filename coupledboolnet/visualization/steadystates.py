"""
Steady-state Distribution Visualizations
"""
from coupledboolnet.bnnetworks.bn import bitstoints,steadystates
from coupledboolnet.bnnetworks.bnstatistics import steadystatesrobust
import networkx as nx
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from math import ceil, sqrt
from numpy import linalg as LA

def transitiondiagram(transitionstates):
    print("\n --- Transition Diagrm of States --- ")

    G = nx.DiGraph()
    G.add_edges_from(transitionstates)
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True)
    plt.show()

def stateevolutionviz(states):
    """
    This needs to be fixed - or just use animation
    """
    print("\n --- Time evolution of States --- ")
    
    states = bitstoints(states)
    stateevolution = pd.DataFrame(states)
    stateevolution.plot(title="Time Evolution of BN")
    plt.xlabel('Time t')
    plt.ylabel('Decimal Encoding of Binary States')
    plt.show()
    
def statedistributionviz(numcells, states, numgenes, stringinfo, gridobj):
    print("\n --- Steady-state Distribution of States --- ")
    print("numcells: ", numcells)

    if (len(states.shape) == 3):
        numcellstodisplay = 9

        counter = 0
        tile =  ceil(sqrt(numcellstodisplay))
        fig = make_subplots(rows=tile, cols=tile, subplot_titles=("Cell 1", "Cell 2", "Cell 3", "Cell 4",
                                                                  "Cell 5", "Cell 6", "Cell 7", "Cell 8", "Cell 9"))
        for i in range(1, tile+1):
            for j in range(1, tile+1):

                ssdistribution, binnum = steadystatesrobust(np.transpose(states[counter, :, :]))
                fig.add_trace(
                    go.Bar(x = binnum, y = ssdistribution, marker_color='rgb(0, 0, 0)' ),
                    row = i, col = j
                )
                #fig.update_yaxes(type="log")
                counter = counter + 1
        fig.update_layout(title_text="Tissue of cell steady-state distributions. T_c={}, h={}".format(gridobj.T_c, gridobj.h),
                          showlegend=False)
        fig.show()
    else:
        ssdistribution, binnum = steadystatesrobust(states)
        fig = px.bar(x= binnum, y = ssdistribution)
        fig.show()

def pcolortest(numcells, states, numgenes):
    numcellline = int(np.sqrt(numcells))

    newmatrix = np.zeros((numcells, states.shape[2]))
    for t in range(states.shape[2]):
        tempstates1 = bitstoints(states[:,:,t]).reshape(numcellline * numcellline)/(2**numgenes)
        newmatrix[:,t] = tempstates1


    newvec = np.zeros(newmatrix.shape[1]-newmatrix.shape[0])
    for t in range(newmatrix.shape[1]-newmatrix.shape[0]):
        newvec[t] = LA.norm(newmatrix[:,t:t+newmatrix.shape[0]].shape)

    print(newvec)
    fig, ax0 = plt.subplots(1,1)
    c=ax0.pcolor(newmatrix[:, :20], vmin=0, vmax=1)
    ax0.set_title('s3_0_1_3_0.8_6.0_6.0')

    #plt.plot(newvec)
    #plt.title('Frobenius norm 100x100 of s3_0_1_3_0.8_5.0_3.0')
    #plt.xlabel('time snapshots')
    #plt.ylabel('Matrix norm')
    plt.show()

    return(newmatrix)


def showanimation(numcells, states, numgenes, defaultnode, gridobj):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11,6))
    cmap = matplotlib.colors.ListedColormap([i for i in range(2**numgenes)])

    gridobj.dT = 10
    plt.title("J={}".format(gridobj.J))
    ims = []
    numcellline = int(np.sqrt(numcells))

    if (len(states.shape) == 3):
        for t in range(states.shape[2]):
            if (t % gridobj.dT== 0):

                ax1.cla()
                tempstates1 = bitstoints(states[:, :, t]).reshape(numcellline, numcellline)
                ax1.imshow(tempstates1, animated=True)
                ax1.set_title("Overall Network: t={}".format(t))
                im1 = ax1.imshow(tempstates1, animated=True)

                ax2.cla()
                tempstates2 = states[:, defaultnode, t].reshape(numcellline, numcellline)
                ax2.imshow(tempstates2, animated=True)
                ax2.set_title("Comm Node: T_c={}".format(gridobj.T_c))
                im2 = ax2.imshow(tempstates2, animated=True)

                ax3.cla()
                tempstates3 = bitstoints(states[:, 1:, t]).reshape(numcellline, numcellline)
                ax3.imshow(tempstates3, animated=True)
                ax3.set_title("Non-comm Network: h={}".format(gridobj.h))
                im3 = ax3.imshow(tempstates3, animated=True)

                ax4.cla()
                tempstates4 = bitstoints(states[:, :, t]).reshape(numcellline, numcellline)
                ax4.imshow(tempstates4, animated=True)
                ax4.set_title("Average Network")
                im4 = ax4.imshow(tempstates3, animated=True)

                ims.append([im1, im2, im3, im4])

                plt.pause(.01)
                ani = animation.ArtistAnimation(fig, ims, blit=False)

            ani.save('ising-model.mp4')

def viewkldbar(kldmatrix):
    plt.bar(range(len(kldmatrix)), kldmatrix)
    plt.title('Pair-wise Kullback-Divergence of tissue network')
    plt.show()