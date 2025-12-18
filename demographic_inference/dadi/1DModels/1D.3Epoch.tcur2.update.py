import matplotlib
matplotlib.use('Agg')  # So plots can render on remote servers
import sys
import argparse
import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum
from numpy import array
import datetime

todaysdate = datetime.datetime.today().strftime('%Y%m%d')
modelName = "1D.1Bottleneck.v3"

# === Parse arguments ===
parser = argparse.ArgumentParser(description='Infer a ' + modelName + ' model from a 1D folded SFS in dadi')
parser.add_argument("--runNum", required=True)
parser.add_argument("--pop", required=True)
parser.add_argument("--mu", required=True, help="Mutation rate (mut/bp/gen)")
parser.add_argument("--L", required=True, help="Number of called neutral sites (monomorphic + polymorphic)")
parser.add_argument("--sfs", required=True, help="Path to FOLDED SFS from easySFS")
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

runNum = str(args.runNum)
pop = str(args.pop)
mu = float(args.mu)
L = float(args.L)
sfs = str(args.sfs)
outdir = str(args.outdir)
maxiter = 100

# === Load folded SFS ===
fs = dadi.Spectrum.from_file(sfs)
if not fs.folded:
    fs = fs.fold()

# === Dadi parameters ===
ns = fs.sample_sizes
pts_l = [ns[0] + 5, ns[0] + 15, ns[0] + 25]

# === Model function with tcur = 2 generations ago ===
Nanc_fixed = 15100
tcur = 2 / (2 * Nanc_fixed)  # ~ 6.34e-5 = ENP ~ 6.62e-5 

def bottleneck(params, ns, pts):
    nuB, nuF, TB = params
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = Integration.one_pop(phi, xx, TB, nuB)      # bottleneck
    phi = Integration.one_pop(phi, xx, tcur, nuF)    # recovery exactly 2 generations ago
    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

param_names = ("nuB", "nuF", "TB")
upper_bound = [10, 10, 1]
lower_bound = [1e-5, 1e-4, 1e-5]
p0 = [1.4, 0.02, 0.13]

func = bottleneck
func_ex = dadi.Numerics.make_extrap_log_func(func)
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)

# === Optimization ===
print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(p0), maxiter=maxiter)
print('Finished optimization **************************************************')

# === Post-processing ===
model = func_ex(popt, ns, pts_l)
ll_model = dadi.Inference.ll_multinom(model, fs)
ll_data = dadi.Inference.ll_multinom(fs, fs)
theta = dadi.Inference.optimal_sfs_scaling(model, fs)

Nanc = theta / (4 * mu * L)
nuB_scaled_dip = popt[0] * Nanc
nuF_scaled_dip = popt[1] * Nanc
TB_scaled_gen = popt[2] * 2 * Nanc
TF_scaled_gen = tcur * 2 * Nanc  # should = 2 generations

scaled_param_names = ("Nanc_FromTheta_scaled_dip", "nuB_scaled_dip", "nuF_scaled_dip", "TB_scaled_gen", "TF_scaled_gen")
scaled_popt = (Nanc, nuB_scaled_dip, nuF_scaled_dip, TB_scaled_gen, TF_scaled_gen)

# === Output parameter file ===
print('Writing out parameters **************************************************')
outputFile = open(f"{outdir}/{pop}.dadi.inference.{modelName}.runNum.{runNum}.output", "w")

header = '\t'.join(param_names) + '\t' + \
         '\t'.join(scaled_param_names) + '\ttheta\tLL\tLL_data\tmodelFunction\tmu\tL\tmaxiter\trunNumber\trundate\tinitialParameters\tupper_bound\tlower_bound'

output = '\t'.join(map(str, popt)) + '\t' + \
         '\t'.join(map(str, scaled_popt)) + f"\t{theta}\t{ll_model}\t{ll_data}\t{func.__name__}\t{mu}\t{L}\t{maxiter}\t{runNum}\t{todaysdate}\t{p0}\t{upper_bound}\t{lower_bound}"

outputFile.write(f"{header}\n{output}\n")
outputFile.close()

# === Output expected SFS ===
print('Writing out SFS **************************************************')
model.to_file(f"{outdir}/{pop}.dadi.inference.{modelName}.runNum.{runNum}.{todaysdate}.expSFS")

# === Plot SFS ===
print('Making plots **************************************************')
import matplotlib.pyplot as plt
fig = plt.figure(1)
dadi.Plotting.plot_1d_comp_multinom(model, fs)
plt.savefig(f"{outdir}/{pop}.dadi.inference.{modelName}.runNum.{runNum}.{todaysdate}.figure.png")

# === Exit ===
sys.exit()
