"""
Boolean Network Class
"""

import numpy as np
from coupledboolnet.bnnetworks.ising import isingsingle, isingsinglefastmetropolis
from numpy import random

class BooleanNetwork():

    def __init__(self, rbnobj, perturbobj, gridobj, importobj=None):
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


        if (gridobj.numcells == 1):
            if importobj.initstate is None:
                self.initstate = random.randint(2, size=(self.grid.numcells, self.n), dtype=bool)
            else:
                self.initstate = importobj.initstate
        elif (gridobj.numcells > 1):
            self.initstate = random.randint(2, size=(self.grid.numcells, self.n), dtype=bool)
            self.gridsize = int(np.sqrt(gridobj.numcells))
            if (gridobj.initconfig == "antiferromagnetic"):

                initconfig = np.zeros((self.gridsize, self.gridsize), dtype=int)
                initconfig[1::2, ::2] = 1
                initconfig[::2, 1::2] = 1
                self.initconfig = initconfig.reshape(gridobj.numcells)

            elif (gridobj.initconfig == "ferromagnetic"):
                initconfig = np.ones((self.gridsize, self.gridsize), dtype=int)
                self.initconfig = initconfig.reshape(gridobj.numcells)

            elif (gridobj.initconfig == "disordered"):
                initconfig = np.random.choice([0, 1], size=(self.gridsize, self.gridsize))
                self.initconfig = initconfig.reshape(gridobj.numcells)
            else:
                assert (gridobj.initconfig == "antiferromagnetic"), "initconfig not supported"
                pass
            self.initconfig = self.initconfig.astype(bool)
            self.initstate[:,0] = self.initconfig

            # self.initstate = random.randint(2, size=(self.grid.numcells, self.n), dtype=bool)
            # a = np.zeros((self.grid.numcells,))
            # a[::2] = True
            # a[1::2] = False
            # self.initstate[:, 0] = a

        if importobj.ttable is None:
            # The bias p needs to be fixed
            # self.ttable = random.randint(2, size=(2**k, n))
            self.ttable = np.zeros((2 ** self.k, self.n), dtype=bool)
            for i in range(2 ** self.k):
                for j in range(self.n):
                    if random.random() < self.p:
                        self.ttable[i, j] = 1
                    else:
                        self.ttable[i, j] = 0
        else:
            self.ttable = importobj.ttable

        if importobj.nv is None:
            self.nv = np.ones((self.k, self.n), dtype=int) * self.k
        else:
            self.nv = importobj.nv

        if importobj.varf is None:
            self.varf = np.zeros((self.k, self.n), dtype=int)
            for i in range(self.n):
                counter = 0
                varftemp = np.array([], dtype=int)
                while (counter < self.k):
                    r = random.randint(1, self.n + 1)
                    if (r not in varftemp):
                        varftemp = np.append(varftemp, r)
                        counter += 1
                self.varf[:, i] = varftemp

        else:
            self.varf = importobj.varf

        self.interactingnode = 0

    def bool_next_all(self, timestep, grid):
        """
        Boolean update for all cells
        """
        numcells = self.grid.numcells

        if (numcells == 1):

            self.states = np.zeros((timestep + 1, numcells * self.n), dtype=bool)
            self.states[0, :] = self.initstate

            # self.states = np.zeros((timestep+1, numcells), dtype=int)
            # self.states[0,:] = bittoint(self.initstate)
            # listallint = self.state_transition()

            if (self.booleanperturbobj.perturbtype == "probabilistic"):
                if (self.booleanperturbobj.booleanperturb != 0):
                    for t in range(1, timestep + 1):  # Might need to come back and fix this
                        for i in range(numcells):
                            self.boolean_next_state(t, i)
                            self.boolean_perturbation(t, i)
                else:
                    for t in range(1, timestep + 1):  # Might need to come back and fix this
                        for i in range(numcells):
                            self.boolean_next_state(t, i)
                            # Needs to be switched to this line eventually.
                            # self.states[t,:] = listallint[self.states[t-1] , 1]

            elif (self.booleanperturbobj.perturbtype == "periodic"):
                pass
            else:
                pass

        elif (numcells > 1):

            dimsize = int(np.sqrt(numcells))
            self.states = np.zeros((numcells, self.n, timestep + 1), dtype=bool)
            self.states[:, :, 0] = self.initstate

            if (self.booleanperturbobj.perturbtype == "probabilistic"):
                if (self.booleanperturbobj.booleanperturb != 0):
                    for t in range(1, timestep + 1):  # Might need to come back and fix this
                        for i in range(numcells):
                            self.boolean_next_state(t, i)
                            self.boolean_perturbation(t, i)
                else:
                    for t in range(1, timestep + 1):  # Might need to come back and fix this
                        for i in range(numcells):
                            self.boolean_next_state(t, i)
                            # Needs to be switched to this line eventually.
                            # self.states[t,:] = listallint[self.states[t-1] , 1]
            elif (self.booleanperturbobj.perturbtype == "stochastic"):
                if (self.booleanperturbobj.booleanperturb != 0):
                    for t in range(1, timestep + 1):  # Might need to come back and fix this
                        for i in range(numcells):
                            self.boolean_next_state(t, i)
                            self.boolean_perturbation(t, i)

                        # inputoutput same node (default)
                        if(self.booleanperturbobj.defaultnode == self.booleanperturbobj.outputnode):
                            nodestates = self.states[:, self.booleanperturbobj.defaultnode, t-1].reshape(dimsize, dimsize)
                            f_NN = self.states[:, self.booleanperturbobj.defaultnode, t].reshape(dimsize, dimsize)
                        else:
                            # inputoutput different node
                            nodestates = self.states[:, self.booleanperturbobj.outputnode, t-1].reshape(dimsize, dimsize)
                            f_NN = self.states[:, self.booleanperturbobj.defaultnode, t].reshape(dimsize, dimsize)

                        #newnodestates = isingsingle(nodestates, self.grid.J, self.grid.T_c, self.grid.h, f_NN)
                        newnodestates = isingsinglefastmetropolis(nodestates, self.grid.J, self.grid.T_c, self.grid.h, f_NN)

                        self.states[:, self.booleanperturbobj.defaultnode, t] = newnodestates.reshape(numcells)

    def boolean_next_state(self, currentime, cellid):
        """
        ***THIS METHOD IS DEPRECATED
        Boolean next step: This needs to be replaced by like state-transition method

        """
        if (self.grid.numcells == 1):

            for genecol in range(self.n):
                tempvarf = self.varf[:, genecol]
                tempvarf = tempvarf[tempvarf != -1]
                tempk = len(tempvarf)
                wiringnode = 0

                for wire in range(len(tempvarf)):
                    wiringnode = wiringnode + 2 ** (tempk - 1 - wire) * self.states[currentime - 1, tempvarf[wire] - 1]

                self.states[currentime, genecol] = self.ttable[wiringnode, genecol]

        elif (self.grid.numcells > 1):
            for genecol in range(self.n):
                tempvarf = self.varf[:, genecol]
                tempvarf = tempvarf[tempvarf != -1]
                tempk = len(tempvarf)
                wiringnode = 0

                # state(t, n_genes)
                # np.zeros((timestep+1, numcells * self.n), dtype=bool)
                # np.zeros((numcells, self.n, timestep+1)

                for wire in range(len(tempvarf)):
                    wiringnode = wiringnode + 2 ** (tempk - 1 - wire) * self.states[
                        cellid, tempvarf[wire] - 1, currentime - 1]

                self.states[cellid, genecol, currentime] = self.ttable[wiringnode, genecol]

    def boolean_perturbation(self, currentime, cellid):
        """
        Boolean perturbation
        """

        if (self.grid.numcells == 1):

            if (self.booleanperturbobj.perturbperturbtype == "probabilistic_and_no_perturb"):
                p = self.booleanperturbobj.booleanperturb
                if (random.random() < 1 - (1 - p) ** (self.n)):
                    for i in range(self.n):
                        if (random.random() < p):
                            self.states[currentime, i] = (self.states[currentime, i] + 1) % 2
                            # This needs to be changed
                            # newperturbed = (inttobitlist(self.states[currentime])[i]+1)%2
                            # self.states[currentime] = bittoint(newperturbed)

            elif (self.booleanperturbobj.perturbperturbtype == "singleperturb_and_no_prob"):
                p = self.booleanperturbobj.booleanperturb
                if (random.random() < p):
                    i = self.booleanperturbobj.defaultnode
                    self.states[currentime, i] = (self.states[currentime, i] + 1) % 2

            elif (self.booleanperturbobj.perturbperturbtype == "probabilistic_and_single_perturb"):
                p = self.booleanperturbobj.booleanperturb
                q = self.booleanperturbobj.commnodeprob

                if (random.random() < 1 - (1 - p) ** (self.n)):
                    for i in range(self.n):
                        if (random.random() < p):
                            self.states[currentime, i] = (self.states[currentime, i] + 1) % 2
                if (random.random() < q):
                    commnode = self.booleanperturbobj.defaultnode
                    self.states[currentime, commnode] = (self.states[currentime, commnode] + 1) % 2

            elif (self.booleanperturbobj.perturbtype == "periodic"):
                p = self.booleanperturbobj.booleanperturb
                pass

            elif (self.booleanperturbobj.perturbtype == "stochastic"):
                pass

        elif (self.grid.numcells > 1):

            if (self.booleanperturbobj.perturbtype == "probabilistic"):
                if (self.booleanperturbobj.perturbperturbtype == "probabilistic_and_no_perturb"):
                    p = self.booleanperturbobj.booleanperturb
                    if (random.random() < 1 - (1 - p) ** (self.n)):
                        for i in range(self.n):
                            if (random.random() < p):
                                self.states[cellid, i, currentime] = (self.states[cellid, i, currentime] + 1) % 2
                                # This needs to be changed
                                # newperturbed = (inttobitlist(self.states[currentime])[i]+1)%2
                                # self.states[currentime] = bittoint(newperturbed)

                elif (self.booleanperturbobj.perturbperturbtype == "singleperturb_and_no_prob"):
                    p = self.booleanperturbobj.booleanperturb
                    if (random.random() < p):
                        i = self.booleanperturbobj.defaultnode
                        self.states[cellid, i, currentime] = (self.states[cellid, i, currentime] + 1) % 2

                elif (self.booleanperturbobj.perturbperturbtype == "probabilistic_and_single_perturb"):
                    p = self.booleanperturbobj.booleanperturb
                    q = self.booleanperturbobj.commnodeprob

                    if (random.random() < 1 - (1 - p) ** (self.n)):
                        for i in range(self.n):
                            if (random.random() < p):
                                self.states[cellid, i, currentime] = (self.states[cellid, i, currentime] + 1) % 2
                    if (random.random() < q):
                        commnode = self.booleanperturbobj.defaultnode
                        self.states[cellid, commnode, currentime] = (self.states[cellid, commnode, currentime] + 1) % 2

            elif (self.booleanperturbobj.perturbtype == "periodic"):
                p = self.booleanperturbobj.booleanperturb
                pass

            elif (self.booleanperturbobj.perturbtype == "stochastic"):
                p = self.booleanperturbobj.booleanperturb
                if (random.random() < 1 - (1 - p) ** (self.n - 1)):
                    for i in range(self.n):
                        if (i != self.booleanperturbobj.defaultnode):
                            if (random.random() < p):
                                self.states[cellid, i, currentime] = (self.states[cellid, i, currentime] + 1) % 2

    def state_transition(self):

        listallbin = inttobitlist(self.n)
        listallint = np.zeros((2 ** self.n, 2), dtype=int)
        listall = np.array([i for i in range(2 ** self.n)])
        listallint[:, 0] = listall

        for i in range(len(listall)):
            nextstate = listallbin[i, :]
            newstate = np.zeros(self.n, dtype=bool)
            for genecol in range(self.n):
                tempvarf = self.varf[:, genecol]
                tempvarf = tempvarf[tempvarf != -1]
                reordered = np.zeros(len(tempvarf), dtype=bool)
                for wire in range(len(tempvarf)):
                    reordered[wire] = nextstate[tempvarf[wire] - 1]
                newstate[genecol] = self.ttable[bittoint(reordered), genecol]
            listallint[i, 1] = bittoint(newstate)

        return (listallint)

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
        # raise Exception("number of genes (n) should equal length of initstate")
        print("Overwritten: number of genes (n) should equal length of initstate")
        rbnobj.n = len(importedata.initstate)

    if (rbnobj.k != max(importedata.nv)):
        print("Overwritten: number of genes (n) should equal length of initstate")
        rbnobj.k = max(importedata.nv)

    if (gridobj.numcells < 1):
        raise Exception("numberof cells should be greater than or equal to one")


def inttobitlist(genes):
    nums = np.array([i for i in range(2 ** genes)], dtype=int)
    bin_nums = ((nums.reshape(-1, 1) & (2 ** np.arange(genes))) != 0).astype(int)

    return (np.array(bin_nums[:, ::-1]))


def bittoint(state):
    num = 0
    for i in range(len(state), 0, -1):
        num = num + state[len(state) - i] * 2 ** (i - 1)
    return (num)

def bitstoints(states):
    bitints = 2 ** np.arange(states.shape[1])[::-1]
    return (states.dot(bitints))

def bitstointsrobust(states):
    # if len(states.shape[2]) == 3:
    bitints = 2 ** np.arange(states.shape[0])[::-1]
    return (bitints.dot(states))

def steadystates(states, genes):
    states = bitstoints(states)
    return (np.array([sum(states == i) for i in range(2 ** genes)]) / states.size,
        np.array([i for i in range(2 ** genes)]))

