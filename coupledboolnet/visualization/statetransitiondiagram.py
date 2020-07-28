"""
State-transition diagram
"""

import networkx as nx 
import matplotlib.pyplot as plt

def transitiondiagram(transitionstates):
   
    print("\n --- Transition Diagrm of States --- ")  

    G = nx.DiGraph()   
    G.add_edges_from(transitionstates)
    nx.draw(G, with_labels = True)
    plt.show()
