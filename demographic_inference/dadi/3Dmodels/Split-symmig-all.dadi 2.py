3# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 2025
@author: Kaden Winspear. Body of script is modified from Annabel Beichman 2D dadi script.
"""
import matplotlib
matplotlib.use('Agg')
import sys
import argparse
import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum
from numpy import array
import datetime
todaysdate=datetime.datetime.today().strftime('%Y%m%d')


#modelName="3D.Split"
modelName="Split-symmig-all"

############### Parse input arguments ########################
parser = argparse.ArgumentParser(description='Infer a '+ modelName +' model from a 3D folded SFS in dadi')
parser.add_argument("--runNum",required=True,help="iteration number (e.g. 1-50)")
parser.add_argument("--pop",required=True,help="population identifier, e.g. 'CA'")
parser.add_argument("--mu",required=True,help="supply mutation rate in mutation/bp/gen")
parser.add_argument("--L",required=True,help="number of called neutral sites that went into making SFS (monomorphic+polymorphic)")
parser.add_argument("--sfs",required=True,help="path to FOLDED SFS in dadi format from easySFS (mask optional)")
parser.add_argument("--outdir",required=True,help="path to output directory")
# usage:
# python 1D.Bottleneck.dadi.py --runNum $i --pop CA --mu 8.64411385098638e-09 --L 4193488 --sfs [path to sfs] --outdir [path to outdir]
args = parser.parse_args()
runNum=str(args.runNum)
pop=str(args.pop)
mu=float(args.mu)
L=float(args.L)
outdir=str(args.outdir)
sfs=str(args.sfs)
maxiter=100
############### Input data ####################################
fs=dadi.Spectrum.from_file(sfs) # this is folded from easy SFS

# check if it's folded, if not folded, fold it
if fs.folded==False:
    fs=fs.fold()
else:
    fs=fs
############### Set up General Dadi Parameters ########################
ns = fs.sample_sizes
max_n = max(ns)
pts_l = [max_n + 5, max_n + 15, max_n + 25]
############### Set up Specific Model -- this will change from script to script #######################

# model adapted from Daniel Portik (daniel.portik@gmail.com) publicly available 3D model pipeline github. https://github.com/dportik

def split_symmig_all(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration is symmetrical between all population pairs (ie 1<->2, 2<->3, and 1<->3).
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    mA: Migration rate between population 1 and population (2,3)
    m1: Migration rate between populations 1 and 2 (2*Na*m) 
    m2: Migration rate between populations 2 and 3 
    m3: Migration rate between populations 1 and 3
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the split of pops 2 and 3 (in units of 2*Na generations).
    """
    #10 parameters
    nu1, nuA, nu2, nu3, mA, m1, m2, m3, T1, T2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)

    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=m3, m31=m3)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

param_names=("nu1","nuA","nu2","nu3","mA","m1","m2","m3","T1","T2")


# 6 parameters: (nu1, nuA, nu2, nu3, T1, T2)
lower_bound = [1e-3, 1e-3, 1e-3, 1e-3, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4]
upper_bound = [100,  100,  100,  100, 20.0, 20.0, 20.0, 20.0, 10.0, 10.0]

# starting guesses (can tweak based on your 2D runs / expectations)
p0 = [1.0, 1.0, 1.0, 1.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

func=split_symmig_all # set the function

############### Carry out optimization (same for any model) ########################
# Make extrapolation function:
func_ex = dadi.Numerics.make_extrap_log_func(func)
# perturb parameters
p0 = dadi.Misc.perturb_params(p0, fold=0.5, upper_bound=upper_bound, lower_bound=lower_bound)
# optimize:
print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(p0), maxiter=maxiter)
print('Finshed optimization **************************************************')

# Calculate the best-fit model AFS.
model = func_ex(popt, ns, pts_l)
# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, fs)
ll_data = dadi.Inference.ll_multinom(fs, fs)  # if you really want it

# calculate best fit theta
theta = dadi.Inference.optimal_sfs_scaling(model, fs)

###### model specific scaling of parameters (will depend on mu and L that you supply) #######
Nanc = theta / (4 * mu * L)

# Present-day and intermediate population sizes in diploids
nu1_scaled_dip = popt[0] * Nanc   # pop 1
nuA_scaled_dip = popt[1] * Nanc   # ancestral 2+3 at first split
nu2_scaled_dip = popt[2] * Nanc   # pop 2
nu3_scaled_dip = popt[3] * Nanc   # pop 3

# Times in generations
T1_scaled_gen = popt[8] * 2 * Nanc
T2_scaled_gen = popt[9] * 2 * Nanc

# Migration parameters:
# dadi uses m_ij = 2*Nanc*m (m = fraction of migrants per generation)
mA_2N = popt[4]
m1_2N = popt[5]
m2_2N = popt[6]
m3_2N = popt[7]

mA_pergen = mA_2N / (2.0 * Nanc)
m1_pergen = m1_2N / (2.0 * Nanc)
m2_pergen = m2_2N / (2.0 * Nanc)
m3_pergen = m3_2N / (2.0 * Nanc)

scaled_param_names = ("Nanc", "N1", "NA_23", "N2", "N3", "T1_gen", "T2_gen", "mA_pergen", "m1_pergen", "m2_pergen", "m3_pergen")

scaled_popt = (Nanc, nu1_scaled_dip, nuA_scaled_dip, nu2_scaled_dip, nu3_scaled_dip, T1_scaled_gen, T2_scaled_gen, mA_pergen, m1_pergen, m2_pergen, m3_pergen)

############### Write out output (same for any model) ########################
print('Writing out parameters **************************************************')

outputFile=open(str(outdir)+"/"+str(pop)+".dadi.inference."+str(modelName)+".runNum."+str(runNum)+".output","w")
# get all param names:
param_names_str='\t'.join(str(x) for x in param_names)
scaled_param_names_str = '\t'.join(str(x) for x in scaled_param_names)
header = (param_names_str + "\t" + scaled_param_names_str +
          "\ttheta\tLL_model\tLL_data\tmodelFunction\tmu\tL\tmaxiter\t" +
          "runNumber\trundate\tinitialParameters\tupper_bound\tlower_bound")

popt_str='\t'.join(str(x) for x in popt) # get opt'd parameters as a tab-delim string
scaled_popt_str='\t'.join(str(x) for x in scaled_popt)
# joint together all the output fields, tab-separated:
modelFunction = func.__name__  # 'split_nomig'

output = [popt_str, scaled_popt_str, theta, ll_model, ll_data, modelFunction, mu, L, maxiter, runNum, todaysdate, p0, upper_bound, lower_bound]
output='\t'.join(str(x) for x in output) # write out all the output fields
# this should result in a 2 row table that could be input into R / concatenated with other runs
outputFile.write(('{0}\n{1}\n').format(header,output))
outputFile.close()

############### Output SFS ########################
print('Writing out SFS **************************************************')

outputSFS=str(outdir)+"/"+str(pop)+".dadi.inference."+str(modelName)+".runNum."+str(runNum)+".expSFS"

model.to_file(outputSFS)


############### Output plot ########################
print('Making plots **************************************************')

import matplotlib.pyplot as plt

outputFigure = (outdir + "/" + pop +
                ".dadi.inference." + modelName +
                ".runNum." + runNum + ".3D_comp.png")

plt.figure(figsize=(15, 5))
dadi.Plotting.plot_3d_comp_multinom(model, fs)
plt.savefig(outputFigure)
plt.close()
###### exit #######
sys.exit()

