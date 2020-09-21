# matlabsaving
import scipy.io
statetransitiontomatlab = rbnP.state_transition()
scipy.io.savemat('./coupledboolnet/MATLAB_TESTING/test.mat', {'allstates':rbnP.states , 'statetransition': statetransitiontomatlab})

Visualize Probabilistic Boolean Network

fig = plt.figure(0)
plt.subplot(ceil(sqrt(numiter)), ceil(sqrt(numiter)), 1)
statedistributionviz(rbnP.states, rbnP.n, perturbobj.booleanperturb)

fig.text(0.5, 0.99, 'Steady State Distributions: Single Node and Perturbation Node 0', ha='center', va='center')
fig.text(0.5, 0.04, 'states', ha='center', va='center')
fig.text(0.06, 0.5, 'counts', ha='center', va='center', rotation='vertical')
plt.show()
