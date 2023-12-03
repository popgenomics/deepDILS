#!/bin/path/python3
import os
import argparse

binpath = '/home/croux/Programmes/deepDILS/scripts/dev'

def run_simulation(outfile, n, S, r, L, Sp, s, fA, Nc, Na, m, Ts, Td, Tsel, width, step):
	# simulations
	command = 'python3 {binpath}/submit_msms.py  -n {n} -S {S} -r {r} -L {L} -Sp {Sp} -s {s} -fA {fA} -Nc {Nc} -Na {Na} -m {m} -Ts {Ts} -Td {Td} -Tsel {Tsel} -output {outfile}.msms'.format(
		binpath=binpath, n=n, S=S, r=r, L=L, Sp=Sp, s=s, fA=fA, Nc=Nc, Na=Na, m=m, Ts=Ts, Td=Td, Tsel=Tsel, outfile=outfile)
	print('\tSimulations:')
	print('\t' + command)
	os.system(command)

	# get trajectory of the advantageous allele
	command = 'python3 {binpath}/get_msms_trajectory.py --infile {outfile}.msms --outfile {outfile}.trajectory'.format(
		binpath=binpath, outfile=outfile)
	print('\n\tGet trajectory of selected mutation:')
	print('\t' + command)
	os.system(command)

	# compute summary statistics
	command = 'python3 {binpath}/msmscalc_onePop_sort_haplo.py --infile {outfile}.msms --nIndiv {n} --nCombParam 1 --regionSize {L} --width {width} --step {step} --nRep 1 --root {outfile}'.format(
		binpath=binpath, outfile=outfile, n=n, L=L, width=width, step=step)
	print('\n\tCompute summary statistics:')
	print('\t' + command)
	os.system(command)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Launches a simulation pipeline whose outputs are: 1) the output of msms (exp_0.msms); 2) the trajectory of the advantageous allele since its introduction Tsel generations ago (exp_0.trajectory); 3) the haplotypes sorted from top to bottom by decreasing frequency (exp_0_0_sorted_rows.txt); 4) the central positions of the sliding windows (exp_0_positions.txt); 5) the summary statistics for each window (exp_0_sumStats.txt).')
	
	parser.add_argument('--outfile', type=str, default='exp_0', help='Root name of the output files (default: %(default)s)')
	parser.add_argument('--n', type=int, default=20, help='Number of sampled gametes (default: %(default)s)')
	parser.add_argument('--S', type=int, default=2000, help='Number of segregating sites (default: %(default)s)')
	parser.add_argument('--r', type=float, default=1e-8, help='Local recombination rate /bp/gen (default: %(default)s)')
	parser.add_argument('--L', type=int, default=100000, help='Length of the chromosome in nucleotides (default: %(default)s)')
	parser.add_argument('--Sp', type=int, default=20000, help='Position of the selected site in nucleotides between 1 and L (default: %(default)s)')
	parser.add_argument('--s', type=float, default=0.001, help='Selective coefficient of the allele A (1+2s; 1+s; 1) (default: %(default)s)')
	parser.add_argument('--fA', type=float, default=0.01, help='Frequency of advantageous allele A when selection starts (default: %(default)s)')
	parser.add_argument('--Nc', type=int, default=100000, help='Number of current individuals (default: %(default)s)')
	parser.add_argument('--Na', type=int, default=10000, help='Number of ancestral individuals (default: %(default)s)')
	parser.add_argument('--Ts', type=int, default=100000, help='Time of speciation (default: %(default)s)')
	parser.add_argument('--Td', type=int, default=25000, help='Time of demographic change (default: %(default)s)')
	parser.add_argument('--Tsel', type=int, default=20000, help='Time at which selection begins (default: %(default)s)')
	parser.add_argument('--m', type=float, default=0, help='Migration rate from pop1 to pop2 backward in time (4.N.m) (default: %(default)s)')
	parser.add_argument('--width', type=float, default=0.05, help='Width of sliding window (in proportion between 0 and 1 of the chromosome length) (default: %(default)s)')
	parser.add_argument('--step', type=float, default=0.05, help='Sliding window step (in proportion between 0 and 1 of the chromosome length) (default: %(default)s)')
	
	args = parser.parse_args()

	# Execute the simulation using arguments from the command line
	run_simulation(args.outfile, args.n, args.S, args.r, args.L, args.Sp, args.s, args.fA, args.Nc, args.Na, args.m, args.Ts, args.Td, args.Tsel, args.width, args.step)

