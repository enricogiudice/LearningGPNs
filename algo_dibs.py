import sys
import numpy as np
import jax.random as random

from dibs.inference import JointDiBS
from dibs.models import DenseNonlinearGaussian, ScaleFreeDAGDistribution

key = random.PRNGKey(0)
key, subk = random.split(key)

n = int(sys.argv[2])
n_particles = int(sys.argv[3])
par = float(sys.argv[4])

data = sys.argv[1].split(',')  # import data and convert from string to array
data_lst = [float(i) for i in data]
data_arr = np.reshape(data_lst, (-1, n), order = 'F')

model = DenseNonlinearGaussian(graph_dist = ScaleFreeDAGDistribution(n), hidden_layers = [5], 
							   sig_param = par, obs_noise = 0.1*par)

post = JointDiBS(x = data_arr, interv_mask = None, inference_model = model)
gs, thetas = post.sample(key = subk, n_particles = n_particles, steps = 2000)  
dibs_mix = post.get_mixture(gs, thetas)

asvec = np.reshape(dibs_mix.g, -1)  # turn adjacency matrix into 1D string; 
print("".join(str(x) for x in asvec))
