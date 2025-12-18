import matplotlib
matplotlib.use('Agg')  # so graphics show up on Hoffman
import sys
import argparse
import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum
from numpy import array
import datetime

todaysdate = datetime.datetime.today().strftime('%Y%m%d')
modelName = "1D.1Bottleneck.v3"

############### Parse input arguments ########################
parser = argparse.ArgumentParser(description='Infer a ' + modelName + ' model from a 1D folded SFS in dadi')
parser.add_argument("--runNum", required=True, help="iteration number (e.g. 1-50)")
parser.add_argument("--pop", required=True, help="population identifier, e.g. 'ENP/GOC'")
parser.add_argument("--mu", required=True, help="mutation rate in mutation/bp/gen")
parser.add_argument("--L", required=True, help="number of called neutral/synonymous sites used for SFS")
parser.add_argument("--sfs", required=True, help="path to folded SFS in dadi format from easySFS")
parser.add_argument("--outdir", required=True, help="output directory path")
args = parser.parse_args()

runNum = str(args.runNum)
pop = str(args.pop)
mu = float(args.mu)
L = float(args.L)
sfs = str(args.sfs)
outdir = str(args.outdir)
maxiter = 100

############### Input data ####################################
fs = dadi.Spectrum.from_file(sfs)
if fs.folded == False:
    fs = fs.fold()

############### General Dadi Parameters ########################
ns = fs.sample_sizes
pts_l = [ns[0] + 5, ns[0] + 15, ns[0] + 25]

############### Model Function (Recovery Fixed to 1 Gen Ago) ########################

# For Nanc = 15770, coalescent time for 1 generation = 1 / (2 * 15770) = 3.17e-5
tcur_1gen = 1 / (2 * 15100)  # ≈ 3.17e-5 = ENP, ESP ≈ 3.31e-5 

def bottleneck(params, ns, pts):
    nuB, nuF, TB = params
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = Integration.one_pop(phi, xx, TB, nuB)         # bottleneck
    phi = Integration.one_pop(phi, xx, tcur_1gen, nuF)  # fixed recovery time = 1 generation
    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

param_names = ("nuB", "nuF", "TB")
upper_bound = [10, 10, 1]
lower_bound = [1e-5, 1e-4, 1e-5]
p0 = [1.4, 0.02, 0.13]

func = bottleneck

############### Optimization ########################
func_ex = dadi.Numerics.make_extrap_log_func(func)
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)
print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(p0), maxiter=maxiter)
print('Finished optimization **************************************************')

model = func_ex(popt, ns, pts_l)
ll_model = dadi.Inference.ll_multinom(model, fs)
ll_data = dadi.Inference.ll_multinom(fs, fs)
theta = dadi.Inference.optimal_sfs_scaling(model, fs)

############### Parameter Scaling ########################
Nanc = theta / (4 * mu * L)
nuB_scaled_dip = popt[0] * Nanc
nuF_scaled_dip = popt[1] * Nanc
TB_scaled_gen = popt[2] * 2 * Nanc
TF_scaled_gen = tcur_1gen * 2 * Nanc  # Should = 1

scaled_param_names = ("Nanc_FromTheta_scaled_dip", "nuB_scaled_dip", "nuF_scaled_dip", "TB_scaled_gen", "TF_scaled_gen")
scaled_popt = (Nanc, nuB_scaled_dip, nuF_scaled_dip, TB_scaled_gen, TF_scaled_gen)

############### Output ########################
print('Writing out parameters **************************************************')
outputFile = open(f"{outdir}/{pop}.dadi.inference.{modelName}.runNum.{runNum}.output", "w")

param_names_str = '\t'.join(param_names)
scaled_param_names_str = '\t'.join(scaled_param_names)
header = f"{param_names_str}\t{scaled_param_names_str}\ttheta\tLL\tLL_data\tmodelFunction\tmu\tL\tmaxiter\trunNumber\trundate\tinitialParameters\tupper_bound\tlower_bound"
popt_str = '\t'.join(str(x) for x in popt)
scaled_popt_str = '\t'.join(str(x) for x in scaled_popt)

output = f"{popt_str}\t{scaled_popt_str}\t{theta}\t{ll_model}\t{ll_data}\t{func.__name__}\t{mu}\t{L}\t{maxiter}\t{runNum}\t{todaysdate}\t{p0}\t{upper_bound}\t{lower_bound}"
outputFile.write(f"{header}\n{output}\n")
outputFile.close()

############### Output SFS ########################
print('Writing out SFS **************************************************')
outputSFS = f"{outdir}/{pop}.dadi.inference.{modelName}.runNum.{runNum}.{todaysdate}.expSFS"
model.to_file(outputSFS)

############### Plot ########################
print('Making plots **************************************************')
import matplotlib.pyplot as plt
fig = plt.figure(1)
outputFigure = f"{outdir}/{pop}.dadi.inference.{modelName}.runNum.{runNum}.{todaysdate}.figure.png"
dadi.Plotting.plot_1d_comp_multinom(model, fs)
plt.savefig(outputFigure)

sys.exit()
