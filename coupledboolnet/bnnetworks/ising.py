"""
Standalone Module for simulating the Ising model
"""
import numpy as np

class Ising():

    def __init__(self, isingvariables, numsteps, initconfig=None):
        """
        Class used to instantiate the Ising model. This object can take on a grid/lattice of pre-configured spins or
        instantiate a random configuration. One can analyze the magnetization of the simulation, or use methods in conjunction
        with bn.py to form coupled Boolean networks by the Ising model. Supply Ising parameters in the GridVariables Class.

        ...
        Attributes
        ----------
        isingvariables : obj
            GridVariables Class

        numsteps : int
            number of timesteps for running the Ising model

        initconfig : float
            optional for supplying initial configuration of spins of the lattice
        """
        assert (isingvariables.numcells > 1), "Number of cells needs to be > 1"

        self.gridsize = int(np.sqrt(isingvariables.numcells))
        self.J = isingvariables.J
        self.h = isingvariables.h
        self.T_c = isingvariables.T_c
        self.hf = 0  # to be modified.
        self.dT = isingvariables.dT

        if (initconfig is None):

            if (isingvariables.initconfig == "antiferromagnetic"):

                initconfig = np.zeros((self.gridsize, self.gridsize), dtype=int)
                initconfig[1::2, ::2] = 1
                initconfig[::2, 1::2] = 1
                initconfig = np.where(initconfig == 0, -1, initconfig)
                self.initconfig = initconfig

            elif (isingvariables.initconfig == "ferromagnetic"):
                self.initconfig = np.ones((self.gridsize, self.gridsize), dtype=int)

            elif (isingvariables.initconfig == "disordered"):
                self.initconfig = np.random.choice([-1, 1], size=(self.gridsize, self.gridsize))

            else:
                assert (isingvariables.init == "antiferromagnetic"), "initconfig not supported"
                pass
        else:
            self.initconfig = initconfig

        self.numsteps = numsteps
        self.kb = 8.617 * (10 ** -5)  # eV/K

        self.dT = isingvariables.dT

    def isingsinglestep(self):
        mc_i = np.random.randint(self.initconfig.shape[0])
        mc_j = np.random.randint(self.initconfig.shape[1])
        N, S, W, E = boundaries(mc_i, mc_j, self.initconfig)

        s = self.J * self.initconfig[mc_i, mc_j]
        nb = (self.initconfig[S, mc_j] + self.initconfig[N, mc_j] + self.initconfig[mc_i, E] + self.initconfig[mc_i, W])

        dE = -2 * s * nb

        if dE < 0:
            s *= -1
        elif np.random.rand() < np.exp((-dE) / self.T_c):
            s *= -1
        self.initconfig[mc_i, mc_j] = s
        return (self.initconfig)

    def isingrunall(self):

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()

        for t in range(self.numsteps):

            mc_i = np.random.randint(self.initconfig.shape[0])
            mc_j = np.random.randint(self.initconfig.shape[1])
            N, S, W, E = boundaries(mc_i, mc_j, self.initconfig)

            s = self.initconfig[mc_i, mc_j]
            nb = (self.initconfig[S, mc_j] + self.initconfig[N, mc_j] + self.initconfig[mc_i, E] + self.initconfig[
                mc_i, W])

            NN_new = -self.J * -self.initconfig[mc_i,mc_j] * nb
            NN_old = -self.J * self.initconfig[mc_i,mc_j] * nb

            dE = NN_new - NN_old

            if dE < 0:
                s = -s
            elif np.random.rand() < np.exp((-dE) / self.T_c):
                s = -s
            self.initconfig[mc_i, mc_j] = s

            if (t % self.dT == 0):
                ax.cla()
                ax.imshow(self.initconfig)
                ax.set_title("Timestep: {}".format(t))
                plt.pause(0.1)


def isingsingle(initconfig, J, T_c, h, f_NN):
    initconfig = np.where(initconfig == 0, -1, initconfig)
    f_NN = np.where(f_NN == 0, -1, f_NN)

    mc_i = np.random.randint(initconfig.shape[0])
    mc_j = np.random.randint(initconfig.shape[1])
    N, S, W, E = boundaries(mc_i, mc_j, initconfig)

    s = initconfig[mc_i, mc_j]
    nb = (initconfig[S, mc_j] + initconfig[N, mc_j] + initconfig[mc_i, E] + initconfig[mc_i, W])

    NN_new = -J * -initconfig[mc_i, mc_j] * nb - h*f_NN[mc_i, mc_j]*s
    NN_old = -J * initconfig[mc_i, mc_j] * nb - h*f_NN[mc_i, mc_j]*s

    dE = NN_new - NN_old

    if dE < 0:
        s = -s
    elif np.random.rand() < np.exp((-dE) / T_c):
        s = -s
    initconfig[mc_i, mc_j] = s
    initconfig = np.where(initconfig == -1, False, initconfig)
    initconfig = initconfig.astype(bool)
    return (initconfig)

def energycomputation(i, j, initconfig, f_NN, J, T_c, h):
    """
    Compute the Hamiltonian of the local site

    ...

    Parameters
    ----------
     i : int
        ith location of where to compute local energy

     j : int
        jth location of where to compute local energy

    initconfig : ndarray
        ndarray of spin configurations at time t

    f_NN : BOOL
        spin value of the proposal

    J : int
        interaction constant of the Hamiltonian (+/-1)

    T_c : float
        temperature used for the Boltzmann distribution

    h : float
        strength of the external field
    """

    N, S, W, E = boundaries(i, j, initconfig )
    s = f_NN[i,j] #singrunall
    nb = initconfig[S, j] + initconfig[N, j] + initconfig[i, E] + initconfig[i, W]

    NN_new = J * s * nb + h * initconfig[i, j] * s
    NN_old = -J * s * nb - h * initconfig[i,j] * s

    dE = NN_new - NN_old

    if dE < 0:
        s = -s
    elif np.random.rand() < np.exp((-dE) / T_c):
        s = -s
    initconfig[i, j] = s

    return initconfig

def isingsinglefastmetropolis(initconfig, J, T_c, h, f_NN):
    # current model initconfig (t-1 interacting node)
    # f_nn (t current node)
    initconfig = np.where(initconfig == 0, -1, initconfig)
    f_NN = np.where(f_NN == 0, -1, f_NN)

    M = f_NN.shape[0]
    N = f_NN.shape[1]

    for i in range(M):
        for j in range(N):
            initconfig = energycomputation(i, j, initconfig, f_NN, J, T_c, h)

    initconfig = np.where(initconfig == -1, False, initconfig)
    initconfig = initconfig.astype(bool)

    return initconfig

def flipspin(spin):
    if (spin == 1):
        spin = -1
    elif (spin == -1):
        spin = 1
    else:
        raise Exception("spin wrong input")
    return (spin)


def boundaries(i, j, matrix):
    N = i - 1;
    S = i + 1;
    E = j + 1;
    W = j - 1;

    if (N == -1):
        N = matrix.shape[0] - 1

    if (S == matrix.shape[0]):
        S = 0

    if (W == -1):
        W = matrix.shape[1] - 1

    if (E == matrix.shape[1]):
        E = 0
    return (N, S, W, E)
