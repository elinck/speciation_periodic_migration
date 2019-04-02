# load libraries
import matplotlib
matplotlib.use('PDF')
import moments
import pylab
import matplotlib.pyplot as plt
import numpy as np
from numpy import array
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D
import sys

# import data
dd = Misc.make_data_dict_vcf(vcf_filename="constant_sim_10005626.trees.vcf",popinfo_filename="id_pop_moments.txt")
pop_ids = ["pop1","pop2"]
projections = [20,20]
data = Spectrum.from_data_dict(dd, pop_ids, projections,polarized=False)
ns=data.sample_sizes
np.set_printoptions(precision=3)

# print info about sfs
print "projection", projections
print "sample sizes", data.sample_sizes
sfs_sum = np.around(data.S(), 2)
print "Sum of SFS = ", sfs_sum, '\n', '\n'

# define model
params=[1,1,1,1]

def IM(params, ns):
    nu1,nu2,T0,m = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T0, m = np.array([[0, m], [m, 0]]))
    return  fs

func=IM
upper_bound = [200,200,200,200]
lower_bound = [1e-3,1e-3,1e-3,1e-3]
params = moments.Misc.perturb_params(params, fold=2, upper_bound=upper_bound,lower_bound=lower_bound)

# run initial optimization
poptg = moments.Inference.optimize_log(params, data, func,lower_bound=lower_bound,
                                upper_bound=upper_bound,verbose=False, maxiter=10)

# run 25 optimizations with less pertubation
for i in range(25):
        poptg=moments.Misc.perturb_params(poptg, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)
        poptg=moments.Inference.optimize_log(poptg, data, func,lower_bound=lower_bound,
                                        upper_bound=upper_bound,verbose=False, maxiter=10)
        sample_sizes = ns
        model=func(poptg, sample_sizes)
        ll_model=moments.Inference.ll_multinom(model, data)
        print('Maximum log composite likelihood: {0}'.format(ll_model))
        theta = moments.Inference.optimal_sfs_scaling(model, data)
        fst=data.Fst()
        IMmod=open('IM_params.txt','a')
        IMmod.write(str(poptg[0])+'\t'+
            str(poptg[1])+'\t'+
            str(poptg[2])+'\t'+
            str(poptg[3])+'\t'+
            str(poptg[4])+'\t'+
            str(poptg[5])+'\t'+
            str(poptg[6])+'\t'+
            str(poptg[7])+'\t'+
            str(theta)+'\t'+
            str(fst)+'\t'+
            str(ll_model)+'\n')
        IMmod.close()
