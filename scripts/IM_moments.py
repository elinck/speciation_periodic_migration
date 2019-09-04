import matplotlib
matplotlib.use('PDF')
import moments
import pylab
import matplotlib.pyplot as plt
import numpy as np
from numpy import array
from moments import Misc,Spectrum,Numerics,Manips,Integration,Demographics1D,Demographics2D
import sys,os
import argparse

os.chdir("/home/cbattey2/speciation_cyclical_migration/")
parser=argparse.ArgumentParser()
parser.add_argument("--infile")
args=parser.parse_args()
#args=argparse.Namespace(infile="simulations/msprime/constant/constant_sim_10441991.trees.vcf")

# import data
print(args.infile)
dd = Misc.make_data_dict_vcf(vcf_filename=args.infile,popinfo_filename="id_pop_moments.txt")
pop_ids = ["1","2"]
projections = [40,40]
data = Spectrum.from_data_dict(dd, pop_ids, projections,polarized=False)
ns=data.sample_sizes
np.set_printoptions(precision=3)

# print info about sfs
print("projection", projections)
print("sample sizes", data.sample_sizes)
sfs_sum = np.around(data.S(), 2)
print("Sum of SFS = ", sfs_sum, '\n', '\n')

# define model
def IM(params, ns):
    nu1,nu2,T0,m = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T0, m = np.array([[0, m], [m, 0]]))
    return  fs
func=IM

upper_bound = [10,10,10,5]
lower_bound = [1e-4,1e-4,1e-4,0]

# run 10 optimizations from random starting parameters
for i in range(20):
    print("starting optimization "+str(i))
    poptg=[np.random.uniform(lower_bound[x],upper_bound[x]) for x in range(4)]
    poptg=moments.Inference.optimize_log(poptg, data, func,lower_bound=lower_bound,
                                    upper_bound=upper_bound,verbose=False, maxiter=10)
    sample_sizes = ns
    model=func(poptg, sample_sizes)
    ll_model=moments.Inference.ll_multinom(model, data)
    print('Maximum log composite likelihood: {0}'.format(ll_model))
    theta = moments.Inference.optimal_sfs_scaling(model, data)
    fst=data.Fst()
    IMmod=open('IM_params_v5.txt','a')
    IMmod.write(
        str(os.path.basename(args.infile))+'\t'+
        str(poptg[0])+'\t'+
        str(poptg[1])+'\t'+
        str(poptg[2])+'\t'+
        str(poptg[3])+'\t'+
        str(theta)+'\t'+
        str(fst)+'\t'+
        str(ll_model)+'\n')
    IMmod.close()

#run with >> parallel -j 30 ~/miniconda2/bin/python scripts/IM_moments.py --infile ~/speciation_cyclical_migration/simulations/msprime/periodic/{1} ::: `ls simulations/msprime/periodic/ | grep ".vcf"`
