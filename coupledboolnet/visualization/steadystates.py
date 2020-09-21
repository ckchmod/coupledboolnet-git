"""
Steady-state Distribution Visualizations
"""
from coupledboolnet.bnnetworks.bn import bitstoints, bittoint ,steadystates
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

def transitiondiagram(transitionstates):
    print("\n --- Transition Diagrm of States --- ")

    G = nx.DiGraph()
    G.add_edges_from(transitionstates)
    nx.draw(G, with_labels=True)
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
    
def statedistributionviz(numcells, states, numgenes, stringinfo):
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

                ssdistribution, binnum = steadystates(np.transpose(states[counter, :, :]),numgenes)
                fig.add_trace(
                    go.Bar(x = binnum, y = ssdistribution, marker_color='rgb(0, 0, 0)' ),
                    row = i, col = j
                )
                counter = counter + 1
        fig.update_layout(title_text="Tissue of cell steady-state distributions", showlegend=False)
        fig.show()
    else:
        ssdistribution, binnum = steadystates(states, numgenes)
        fig = px.bar(x= binnum, y = ssdistribution)
        fig.show()

def showanimation(numcells, states, numgenes, timestep, dt):
    fig, ax = plt.subplots()
    cmap = matplotlib.colors.ListedColormap([i for i in range(2**numgenes)])

    ims = []
    for t in range(timestep):

       if (t % dt == 0):

           ax.cla()
           states = bitstoints(states[:, :, t]).reshape(int(np.sqrt(numcells)),
                                                             int(np.sqrt(numcells)))
           ax.imshow(states, animated=True)
           ax.set_title("Timestep: {}".format(t))
           im = ax.imshow(states, animated=True)
           ims.append([im])
           plt.pause(.1)

           ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
                                           repeat_delay=1000)

    ani.save('ising-model.mp4')
    pass