import os
import argparse
import math 

def main(args):
	description = 'Script that launches msms for a given combination of parameters.\nSimulates a focal population whose size can vary over time, and can receive migrants from a sister population.\nDirectional selection can also be involved.\n'
	print(args)
	

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Script that launches msms for a given combination of parameters. Simulates a focal population whose size can vary over time, and can receive migrants from a sister population. Directional selection can also be involved.')

	# arguments with default values 
	parser.add_argument('-n', '--number_sampled_gametes', type=int, default=20, help='Number of sampled gametes (default: %(default)s)')
	parser.add_argument('-S', '--number_segregating_sites', type=int, default=100, help='Number of segregating sites (default: %(default)s)')
	parser.add_argument('-r', '--recombination_rate', type=float, default=1e-8, help='Local recombination rate (default: %(default)s)')
	parser.add_argument('-L', '--chromosome_length', type=int, default=100000, help='Length of the chromosome in nucleotides (default: %(default)s)')
	parser.add_argument('-Sp', '--selected_site_position', type=int, default=20000, help='Position of the selected site in nucleotides (default: %(default)s)')
	parser.add_argument('-s', '--selective_coefficient', type=float, default=0.001, help='Selective coefficient of the allele A (default: %(default)s)')
	parser.add_argument('-fA', '--frequency_adv', type=float, default=0.01, help='Frequency of advantageous allele A when selection starts (default: %(default)s)')
	parser.add_argument('-Nc', '--current_individuals', type=int, default=100000, help='Number of current individuals (default: %(default)s)')
	parser.add_argument('-Na', '--ancestral_individuals', type=int, default=10000, help='Number of ancestral individuals (default: %(default)s)')
	parser.add_argument('-m', '--migration_rate', type=float, default=0.0, help='Migration rate from pop1 to pop2 backward in time (default: %(default)s)')
	parser.add_argument('-Ts', '--time_speciation', type=int, default=100000, help='Time of speciation (default: %(default)s)')
	parser.add_argument('-Td', '--time_demographic_change', type=int, default=25000, help='Time of demographic change (default: %(default)s)')
	parser.add_argument('-Tsel', '--time_selection_starts', type=int, default=20000, help='Time at which selection begins (default: %(default)s)')
	parser.add_argument('-output', '--outputfile_name', type=str, default='output.msms', help='Name of the output file where simulations are recorded (default: %(default)s)')

	args = parser.parse_args()
	main(args)

# outputfile name
outfile = args.outputfile_name

# sampling
n = args.number_sampled_gametes

# genomic
S = args.number_segregating_sites
r = args.recombination_rate
L = args.chromosome_length
Sp = args.selected_site_position
s = args.selective_coefficient
freq_adv = args.frequency_adv

# demography
N_cur = args.current_individuals
N_anc = args.ancestral_individuals
m_12 = args.migration_rate

# events
T_speciation = args.time_speciation
T_demography = args.time_demographic_change
T_selection = args.time_selection_starts

# rescaling
s_aa = 2*N_cur*0
s_Aa = 2*N_cur*s
s_AA = 2*N_cur*2*s

# getting the format of positions
order_magnitude = math.floor(math.log10(L)) + 1
format_pos = 10**(-order_magnitude)
format_pos ="{:.{}f}".format(format_pos, order_magnitude)

# command line
# if selection starts before the demographic change (forward in time)
if T_selection > T_demography:
	# the rescale 2.N.s in order to keep 's' constant over time
	command_line = 'msms -N {N} -ms {n} 1 -s {S} -r {rho} {L} -I 2 {n} 0 0 -m 1 2 {m_12} -m 2 1 0 -en {T_demography} 1 {N_anc} -ej {T_speciation} 2 1 -SI {T_selection} 2 {freq_adv} 0 -Sc 0 1 {SAA} {SAa} {Saa} -Sc {T_demography} 1 {SAA_d} {SAa_d} {Saa_d} -Sc 0 2 0 0 0 -Sp {Sp} -oFormat {format_pos} -oTrace > {outfile}'.format(n=n, S=S, rho=4*N_cur*r*L, L=L, m_12=m_12, T_selection=T_selection/(4*N_cur), T_demography=T_demography/(4*N_cur), N_anc=N_anc/N_cur, T_speciation=T_speciation/(4*N_cur), freq_adv=freq_adv, SAA=s_AA, SAa=s_Aa, Saa=s_aa, SAA_d=s_AA*N_anc/N_cur, SAa_d=s_Aa*N_anc/N_cur, Saa_d=s_aa*N_anc/N_cur, Sp=Sp/L, N=N_cur, format_pos=format_pos, outfile=outfile)
else:
	command_line = 'msms -N {N} -ms {n} 1 -s {S} -r {rho} {L} -I 2 {n} 0 0 -m 1 2 {m_12} -m 2 1 0 -en {T_demography} 1 {N_anc} -ej {T_speciation} 2 1 -SI {T_selection} 2 {freq_adv} 0 -Sc 0 1 {SAA} {SAa} {Saa} -Sc 0 2 0 0 0 -Sp {Sp} -oFormat {format_pos} -oTrace > {outfile}'.format(n=n, S=S, rho=4*N_cur*r*L, L=L, m_12=m_12, T_selection=T_selection/(4*N_cur), T_demography=T_demography/(4*N_cur), N_anc=N_anc/N_cur, T_speciation=T_speciation/(4*N_cur), freq_adv=freq_adv, SAA=s_AA, SAa=s_Aa, Saa=s_aa, Sp=Sp/L, N=N_cur, format_pos=format_pos, outfile=outfile)
print('\n{command_line}\n'.format(command_line=command_line))


# perform simulations
max_attempts = 10
attempt = 0

while attempt < max_attempts:
	print('\tattempt #{0}'.format(attempt+1))
	return_code = os.system(command_line)

	if return_code == 0:
		print('\t\tsimulation successfully performed')
		break
	else:
		attempt += 1
		print('\t\tsimulation failed')

if attempt == max_attempts:
	print('\n\tSimulation failed after 10 attempts')

