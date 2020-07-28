"""
Boolean Network Class
"""
import numpy as np

class Ising():
    
    def __init__(self, isingvariables, numsteps, initconfig = None):
        """
        Ising Method

        Parameters
        ----------
        """
        assert(isingvariables.numcells > 1), "Number of cells needs to be > 1"
        
        self.gridsize = int(np.sqrt(isingvariables.numcells))
        self.J = isingvariables.J
        self.h = isingvariables.h
        self.T_c = isingvariables.T_c
        self.hf = 0 #to be modified.
        
        if( initconfig is None):
        
            if (isingvariables.initconfig == "antiferromagnetic"):
                
                initconfig = np.zeros((self.gridsize, self.gridsize), dtype = int)
                initconfig[1::2, ::2] = 1
                initconfig[::2, 1::2] = 1
                initconfig = np.where(initconfig==0, -1, initconfig) 
                self.initconfig = initconfig
                
            elif (isingvariables.initconfig == "ferromagnetic"):
                self.initconfig = np.ones((self.gridsize, self.gridsize), dtype = int)
                
            elif (isingvariables.initconfig == "disordered"):
                self.initconfig = np.random.choice([-1, 1], size=(self.gridsize,self.gridsize))
                
            else:
                assert(isingvariables.init == "antiferromagnetic"), "initconfig not supported"
                pass
        else:
            self.initconfig = initconfig 
        
        self.numsteps = numsteps
        self.kb = 8.617 * (10**-5) # eV/K
    
    def isingrun(self):
        
        Energy_0 = np.zeros((self.gridsize, self.gridsize), dtype = bool)
        
        for i in range(self.initconfig.shape[0]):
            for j in range(self.initconfig.shape[1]):                
                N, S, W, E = boundaries(i,j,self.initconfig)
                Energy_0[i,j] =  -1 * self.initconfig[i,j] * \
                    (self.initconfig[S,j] + self.initconfig[N,j] + self.initconfig[i,E] + self.initconfig[i,W]);

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()


        for t in range(self.numsteps):

            mc_i = np.random.randint(self.initconfig.shape[0])
            mc_j = np.random.randint(self.initconfig.shape[1])
            N, S, W, E = boundaries(mc_i,mc_j,self.initconfig)
            
            #dE = energycomputation(self.initconfig[mc_i, mc_j], N, S, W, E, self.J, self.h, self.hf)
            s =  self.initconfig[mc_i, mc_j]
            nb = (self.initconfig[S,mc_j] + self.initconfig[N,mc_j] + self.initconfig[mc_i,E] + self.initconfig[mc_i,W])
            
            dE = 2*s*nb
            if dE < 0:
                s *= -1
            elif np.random.rand() < np.exp((-1*dE) /self.T_c):
                s *= -1
            self.initconfig[mc_i, mc_j] = s

            if (t % 100 == 0):
                """
                plot every 100 timesteps
                """
                ax.cla()
                ax.imshow(self.initconfig)
                ax.set_title("Timestep: {}".format(t))
                plt.pause(0.1)
    
def flipspin(spin):
    if (spin == 1):
        spin = -1
    elif (spin == -1):
        spin = 1
    else:
        raise Exception("spin wrong input")
    return(spin)
        
    
def boundaries(i,j, matrix):

    N = i-1; S = i+1; E = j+1; W = j-1;
    if (N == -1): 
        N = matrix.shape[0]-1
        
    if (S == matrix.shape[0]):
        S = 0
    
    if (W == -1): 
        W = matrix.shape[1]-1
    
    if (E == matrix.shape[1]):
        E = 0
    return(N, S, W, E)

         
def energycomputation(s_i , N_S_W_E, J, h_interac, h_f):
    
    N_S_W_E_temp = N_S_W_E;
    s_i_temp = s_i;
    h_f_temp = h_f; 
    
    # s =  self.initconfig[mc_i, mc_j]
    # nb = (self.initconfig[S,mc_j] + self.initconfig[N,mc_j] + self.initconfig[mc_i,E] + self.initconfig[mc_i,W])
          
    # (N_S_W_E[0] + N_S_W_E_temp[1] + N_S_W_E_temp[2] + N_S_W_E_temp[3]))
    # 2*s*nb
    
    
    energy = -2 * s_i * (N_S_W_E_temp[0] + N_S_W_E_temp[1] + N_S_W_E_temp[2] + N_S_W_E_temp[3])
    return energy
    
