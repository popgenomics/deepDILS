import os
import random
import argparse
from math import isnan

# pathway to the bin directory
binpath = '/home/croux/Programmes/deepDILS/scripts/dev'

# threshold of the frequency of the selected allele. <threshold: restart simulation; else: retain the simulation
threshold = 0.9

# Create the parser
parser = argparse.ArgumentParser(description='Simulation of one chromosome following random demographic and selective parameters from prior distributions specified by the user.')

parser.add_argument('--outfile', type=str, default='exp_0', help='Root name of the output files')
parser.add_argument('--model', type=str, default='EXP', help='Demographic model to simulate')
parser.add_argument('-n', type=int, default=40, help='Number of sampled gametes')
parser.add_argument('-S', type=int, default=2000, help='Number of segregating sites')
parser.add_argument('-L', type=int, default=100000, help='Length of the chromosome in nucleotides')
parser.add_argument('-r', type=float, default=1e-8, help='Local recombination rate /bp/gen')
parser.add_argument('-fA', type=float, default=0.05, help='Frequency of the selected allele when selection begins')
parser.add_argument('-m', type=float, default=0.4, help='Migration rate from pop1 to pop2 backward in time (4.N.m)')
parser.add_argument('--width', type=float, default=0.05, help='Width of sliding window')
parser.add_argument('--step', type=float, default=0.025, help='Sliding window step')

parser.add_argument('--min_N', type=int, default=500, help='Minimum number of individuals (PRIOR_min)')
parser.add_argument('--max_N', type=int, default=50000, help='Maximum number of individuals (PRIOR_max)')
parser.add_argument('--min_T', type=int, default=1000, help='Minimum time of split in generations(PRIOR_min)')
parser.add_argument('--max_T', type=int, default=500000, help='Maximum time of split in generations(PRIOR_max)')
parser.add_argument('--min_Ns', type=int, default=10, help='Minimum value for N.s (PRIOR_min)')
parser.add_argument('--max_Ns', type=int, default=100, help='Maximum value for N.s (PRIOR_max)')
parser.add_argument('--scalar_N_min', type=float, default=0.02, help='Minimum value for the ratio smallest Ne over bigger Ne (PRIOR_min)')
parser.add_argument('--scalar_N_max', type=float, default=0.2, help='Maximum value for the ratio smallest Ne over bigger Ne (PRIOR_max)')
parser.add_argument('--scalar_Tsplit_min', type=float, default=0.1, help='Minimum value for the ratio Tsplit over Ne_ancestral (PRIOR_min)')
parser.add_argument('--scalar_Tsplit_max', type=float, default=10, help='Maximum value for the ratio Tsplit over Ne_ancestral (PRIOR_max)')
parser.add_argument('--scalar_Tevent_min', type=float, default=0.2, help="Minimum value for the ratio Tevent over Tsplit (PRIOR_min). The event is 'demographic change' or 'selection start'.")
parser.add_argument('--scalar_Tevent_max', type=float, default=0.5, help="Maximum value for the ratio Tevent over Tsplit (PRIOR_max). The event is 'demographic change' or 'selection start'.")

# Get the arguments
args = parser.parse_args()

def test_trajectory(outfile, threshold, N):
	tested_trajectory = 0
	T_25 = float('nan')
	T_50 = float('nan')
	T_75 = float('nan')
	T_99 = float('nan')

	with open('{outfile}.trajectory'.format(outfile=outfile), "r") as infile:
		header = infile.readline()

		for line in infile:
			line = line.strip().split('\t')
			generation = float(line[0]) * 4 * N
			fA_i = float(line[2])

			if isnan(T_25) and fA_i >= 0.25:
				T_25 = generation
			if isnan(T_50) and fA_i >= 0.5:
				T_50 = generation
			if isnan(T_75) and fA_i >= 0.75:
				T_75 = generation
			if isnan(T_99) and fA_i >= 0.99:
				T_99 = generation

			if fA_i >= threshold:
				tested_trajectory = 1

	return tested_trajectory, T_25, T_50, T_75, T_99

def define_parameters(model, min_N, max_N, min_Ns, max_Ns, scalar_Tsplit_min, scalar_Tsplit_max, scalar_Tevent_min, scalar_Tevent_max, L, scalar_N_min, scalar_N_max):
	N_sampled = random.randint(min_N, max_N)
	Ns_sampled = random.uniform(min_Ns, max_Ns)
	scalar_T_split = random.uniform(scalar_Tsplit_min, scalar_Tsplit_max)
	scalar_T_dem = random.uniform(scalar_Tevent_min, scalar_Tevent_max)
	scalar_T_selection = random.uniform(scalar_Tevent_min, scalar_Tevent_max)

	# Position of the selected site in nucleotides between 1 and L
	Sp = random.randint(1, L)

	if model == 'CST':
		Nc = N_sampled
		Na = Nc
		Tsplit = int(scalar_T_split * Na)
		Tdem = int(scalar_T_dem * Tsplit)
		Tselection = int(scalar_T_selection * Tsplit)
		# Selective coefficient of the allele A (1+2s; 1+s; 1)
		s = Ns_sampled / Nc

	if model == 'BTL':
		scalar_N = random.uniform(scalar_N_min, scalar_N_max)
		Na = N_sampled
		Nc = int(Na * scalar_N)
		Tsplit = int(scalar_T_split * Na)
		Tdem = int(scalar_T_dem * Tsplit)
		Tselection = int(scalar_T_selection * Tsplit)
		# Selective coefficient of the allele A (1+2s; 1+s; 1)
		s = Ns_sampled / Nc

	if model == 'EXP':
		scalar_N = random.uniform(scalar_N_min, scalar_N_max)
		Nc = N_sampled
		Na = int(Nc * scalar_N)
		Tsplit = int(scalar_T_split * Na)
		Tdem = int(scalar_T_dem * Tsplit)
		Tselection = int(scalar_T_selection * Tsplit)
		# Selective coefficient of the allele A (1+2s; 1+s; 1)
		s = Ns_sampled / Na
	
	return(Nc, Na, Tsplit, Tdem, Tselection, s, Sp)

# conditions to stop the program
return_code = 1
tested_trajectory = 0

# loop until we have a successful fixation
sim_tmp = 0
while return_code != 0 or tested_trajectory != 1:
	sim_tmp += 1
	print('Try #{sim_tmp}'.format(sim_tmp=sim_tmp))
	
	# get parameters for simulations
	Nc, Na, Tsplit, Tdem, Tselection, s, Sp = define_parameters(args.model, args.min_N, args.max_N, args.min_Ns, args.max_Ns, args.scalar_Tsplit_min, args.scalar_Tsplit_max, args.scalar_Tevent_min, args.scalar_Tevent_max, args.L, args.scalar_N_min, args.scalar_N_max)

	commande = 'python3 {binpath}/pipeline_msms.py  --outfile {outfile} --n {n} --S {S} --r {r} --L {L} --Sp {Sp} --s {s} --fA {fA} --Nc {Nc} --Na {Na} --Ts {Tsplit} --Td {Tdem} --Tsel {Tselection} --m {m} --width {width} --step {step}'.format(
		binpath=binpath, outfile=args.outfile,
		n=args.n, S=args.S, r=args.r, L=args.L, fA=args.fA,
		m=args.m,
		Nc=Nc, Na=Na, Tsplit=Tsplit, Tdem=Tdem, Tselection=Tselection, s=s, Sp=Sp, 
		width=args.width, step=args.step)
	return_code = os.system(commande)
	if return_code == 0:
		tested_trajectory, T_25, T_50, T_75, T_99 = test_trajectory(args.outfile, threshold, Nc)

