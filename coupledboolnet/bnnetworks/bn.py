"""
Boolean Network Class
"""

import numpy as np
from numpy import random

class BooleanNetwork():

    def __init__(self, rbnobj, perturbobj, gridobj, importobj = None):
        """
        Boolean Network Instantiation

        Parameters
        ----------
        """
        
        self.n = rbnobj.n
        self.k = rbnobj.k
        self.p = rbnobj.p
        self.booleanperturbobj = perturbobj
        self.grid = gridobj

        if importobj.initstate is None:
            self.initstate = random.randint(2, size=(self.grid.numcells, self.n), dtype=bool)
        else:
            self.initstate = importobj.initstate

        if importobj.ttable is None:
            # The bias p needs to be fixed
            #self.ttable = random.randint(2, size=(2**k, n))
            self.ttable = np.zeros((2**self.k, self.n), dtype=bool)
            for i in range(2**self.k):
                for j in range(self.n):
                    if random.random() < self.p:
                        self.ttable[i,j] = 1
                    else:
                        self.ttable[i,j] = 0
        else:
            self.ttable = importobj.ttable

        if importobj.nv is None:
            self.nv =  np.ones((self.k, self.n), dtype=int) * self.k
        else:
            self.nv = importobj.nv

        if importobj.varf is None:
            self.varf = np.zeros((self.k, self.n), dtype=int)
            for i in range (self.n):
                counter = 0 
                varftemp = np.array([],dtype=int)
                while (counter < self.k):
                    r=random.randint(1, self.n+1)
                    if (r not in varftemp): 
                        varftemp = np.append(varftemp, r)
                        counter += 1
                self.varf[:,i] = varftemp
            
        else:
            self.varf = importobj.varf

        self.interactingnode = 0


    def bool_next_all(self, timestep, grid):
        """
        Boolean update for all cells
        """

        numcells = self.grid.numcells

        self.states = np.zeros((timestep+1, numcells), dtype=int)
        self.states[0,:] = bittoint(self.initstate) 

        listallint = self.state_transition()
      
        if (self.booleanperturbobj.perturbtype == "probabilistic"):

            if (self.booleanperturbobj.booleanperturb != 0):
                for t in range(1, timestep+1): #Might need to come back and fix this
                    for i in range(numcells):
                        self.boolean_next_state(t)
                        self.boolean_perturbation(t)
            else:
                for t in range(1, timestep+1): #Might need to come back and fix this
                    for i in range(numcells):                       
                        self.states[t,:] = listallint[self.states[t-1] , 1]
            
        elif (self.booleanperturbobj.perturbtype == "periodic"):
            pass
        else:
            pass

    def boolean_next_state(self, currentime):
        """
        ***THIS METHOD IS DEPRECATED***
        Boolean next step: This needs to be replaced by like state-transition method
        
        """

        for genecol in range(self.n):
            tempvarf = self.varf[:,genecol] 
            tempvarf = tempvarf[tempvarf != -1]
            tempk = len(tempvarf) 
            wiringnode = 0
            
            for wire in range(len(tempvarf)):
                wiringnode = wiringnode + 2**(tempk -1 - wire)* self.states[currentime-1, tempvarf[wire]-1]

            self.states[currentime, genecol] = self.ttable[wiringnode, genecol]

    def boolean_perturbation(self, currentime):
        """
        ***Perturbaiton needs to be updated***

        Boolean perturbation
        """
        
        if (self.booleanperturbobj.perturbperturbtype == "probabilistic_and_no_perturb"):
            p = self.booleanperturbobj.booleanperturb
            if (random.random() < 1-(1-p)**(self.n)):
                for i in range(self.n):
                    if (random.random() < p):
                        self.states[currentime,i] = (self.states[currentime,i]+1)%2
                        
        elif(self.booleanperturbobj.perturbperturbtype == "singleperturb_and_no_prob"):
            p = self.booleanperturbobj.booleanperturb
            if (random.random() < p):
                i = self.booleanperturbobj.defaultnode
                self.states[currentime,i] = (self.states[currentime,i]+1)%2
                
        elif (self.booleanperturbobj.perturbperturbtype == "probabilistic_and_single_perturb"):
            p = self.booleanperturbobj.booleanperturb
            q = self.booleanperturbobj.commnodeprob

            if (random.random() < 1-(1-p)**(self.n)):
                for i in range(self.n):
                    if (random.random() < p):
                        self.states[currentime,i] = (self.states[currentime,i]+1)%2
            if (random.random() < q):
                commnode = self.booleanperturbobj.defaultnode
                self.states[currentime, commnode] = (self.states[currentime, commnode]+1)%2
        
        elif(self.booleanperturbobj.perturbtype == "periodic"):
            p = self.booleanperturbobj.booleanperturb            
            pass
        
        else:
            assert(self.booleanperturbobj.perturbtype == "stochastic")
            pass
           
        
    def state_transition(self):
        
        listallbin = inttobitlist(self.n)
        listallint = np.zeros((2**self.n, 2), dtype=int)
        listall = np.array([i for i in range(2**self.n)])
        listallint[:,0] = listall
        
        for i in range(len(listall)):
            nextstate = listallbin[i,:]
            newstate = np.zeros(self.n, dtype=bool)
            for genecol in range(self.n):
                tempvarf = self.varf[:,genecol]
                tempvarf = tempvarf[tempvarf != -1]
                reordered = np.zeros(len(tempvarf), dtype=bool)
                for wire in range(len(tempvarf)):
                    reordered[wire] = nextstate[tempvarf[wire]-1]
                newstate[genecol] = self.ttable[bittoint(reordered), genecol]
            listallint[i,1]=bittoint(newstate)    
 
        return(listallint)
    
        
    class Intercelluluar:
        """
        Interaction Rules
        """

        def __init__(self, h, J):
            pass

def inputtest(rbnobj, importedata, perturbobj, gridobj, timestep):
    """
    Class function to test all the inputs are correct
    """
    if rbnobj.k > rbnobj.n:
        raise Exception("numberof wirings (k) should not be greater than number of genes (n)")
    
    if (rbnobj.n != len(importedata.initstate)):
        #raise Exception("number of genes (n) should equal length of initstate")
        print("Overwritten: number of genes (n) should equal length of initstate")
        rbnobj.n = len(importedata.initstate)
    if (rbnobj.k != max(importedata.nv)):
        print("Overwritten: number of genes (n) should equal length of initstate")
        rbnobj.k = max(importedata.nv)
    
def inttobitlist(genes):

    nums = np.array([i for i in range(2**genes)], dtype=int)
    bin_nums = ((nums.reshape(-1,1) & (2**np.arange(genes))) != 0).astype(int)
    
    return(np.array(bin_nums[:,::-1]))

def bittoint(state):

    num = 0
    for i in range(len(state),0,-1):
        num = num + state[len(state)-i] * 2**(i-1)
    return(num)


def bitstoints(states):

    bitints = 2**np.arange(states.shape[1])[::-1] 
    return(states.dot(bitints) )

def steadystates(states, genes):
    states = bitstoints(states)
    return ( np.array([sum(states==i) for i in range(2**genes)]) / states.size, 
            np.array([i for i in range(2**genes)]))



