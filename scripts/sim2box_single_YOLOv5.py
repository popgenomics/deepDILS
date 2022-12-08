from matplotlib import pyplot
from matplotlib import cm
import os
import numpy as np
from numpy import nanmean
from numpy import nan
from numpy import isnan
from pandas import isnull
from numpy import nanmin
from numpy import nanmax
import cv2
import gc
import sys

epsilon = 0.00001

scalar = 1 # if scalar = 0.5, then the box is defined by 0.5 theta; if scalar = 1 then the box is defined by theta
binary = cm.get_cmap('binary', 512)

help_message = '\tpython scripts that works after the SLiM simulations to produce PNG files representing the data along the simulated chromosome\n'
help_message += '\tto run this script for instance on the simulation producing the following files: 123_parameters.txt 123_positions.txt 123_sumStats.txt 123_trees.txt\n'
help_message += '\t\tpython3 sim2box_single.py dpi=300 datapath=/home/croux/Programmes/yolo_box/simulations simulation=123 object=posSelection\n'

# version using a loop over files

## ARGUMENTS
dpi=nan
datapath=nan
simulation_target=nan
object_to_recognize=nan
for arg in sys.argv:
	arg = arg.split('=')
	if arg[0] == 'dpi':
		dpi=int(arg[1]) # dpi = 300
	if arg[0] == 'datapath':
		datapath=arg[1] # datapath='/home/croux/Programmes/yolo_box/simulations'
	if arg[0] == 'simulation':
		simulation_target=arg[1] # simulation_target = '5'
	if arg[0] == 'object':
		object_to_recognize=arg[1] # neutral; posSelection
	if arg[0] == 'theta':
		theta=int(arg[1]) # 0 (don't use theta=4.N.u for the expected theta, but use the neutral simulations); 1 (use theta=4.N.u)
	if arg[0] == 'phasing':
		phasing=int(arg[1]) # 0: don't use the statistics associated to LD; 1: use the LD statistics
	if arg[0] == 'plotStats':
		plotStats=int(arg[1]) # 0: don't plot individual stats (pi.png, theta.png, etc....); 1: plots the individuals stats
	if arg[0] == 'getAllData':
		getAllData=int(arg[1]) # 0: in order to get the total range of stat. values, don't read all simulations but a reference table; 1: read all simulations and produce a reference table for standardization

if isnan(dpi)==True:
	print('\n\ta value of dpi has to be specified\n')
	print(help_message)
	sys.exit(-1)
if isnull(datapath)==True:
	print('\n\ta datapath has to be specified\n')
	print(help_message)
	sys.exit(-1)
if isnull(simulation_target)==True:
	print("\n\ta simulation's ID has to be specified\n")
	print(help_message)
	sys.exit(-1)
if isnull(object_to_recognize)==True:
	print("\n\tthe object to recognize has to be specified\n")
	print(help_message)
	sys.exit(-1)

if object_to_recognize=='neutral':
	object_to_recognize=0
else:
	if object_to_recognize=='posSelection':
		object_to_recognize=1	
## END OF TMP

## START OF FUNCTIONS
def getRanges():
	ranges_stats = {}
	infile = open('range_stats.txt', 'r')
	for line in infile:
		tmp = line.strip().split(':')
		ranges_stats[tmp[0]] = float(tmp[1].replace(' ', ''))
	return(ranges_stats)

def getSNP(snp, positions):
	# function usefull for the rawdata
	# it 1. gets the position of 'snp' within the vector 'positions'; 2. converts it in a float in [0,1]
	# snp is a float readen from the ms output (positions) (in [0,1])
	# positions is a vector of positions (all in [0,1])
	distances = [ abs(i-snp) for i in positions ]
	closest_snp = distances.index(min(distances))
	
	return(closest_snp / (1.0*len(positions)))

def plot_globalPic(data, colormaps, vmin, vmax, iteration, model):
	"""
	Helper function to plot data with associated colormap.
	"""
#	data = np.random.randn(30, 30)
	n = len(colormaps)
	fig, axs = pyplot.subplots(1, n, figsize=(n * 2 + 2, 3), constrained_layout=True, squeeze=False)
	for [ax, cmap] in zip(axs.flat, colormaps):
		psm = ax.pcolormesh(data, cmap=cmap, rasterized=True, vmin=vmin, vmax=vmax)
#		fig.colorbar(psm, ax=ax)
	pyplot.axis("off")
	pyplot.savefig("{iteration}_{model}_globalPic.jpg".format(iteration=iteration, model=model), bbox_inches='tight', pad_inches = 0, dpi=dpi)
	#pyplot.close()
	gc.collect()
	im = cv2.imread("{iteration}_{model}_globalPic.jpg".format(iteration=iteration, model=model))
	res={}
	res['width'] = im.shape[1]
	res['height'] = im.shape[0]
	return(res)
#	pyplot.show()

def plot_genome(x, y, xlim, ylim, iteration, root):
	# function that plot x and y in a png file of kind 123_pi.png ({iteration}_{root}.png})
	# x=[0,1,2,3]
	# y=[1,2,1,2]
	# xlim=[-1, 4]
	# ylim=[0, 3]
	fig = pyplot.figure(num=1)
	pyplot.scatter(x, y,
		       marker = 's', edgecolors = 'none')
	pyplot.xlim(xlim[0], xlim[1])
	pyplot.ylim(ylim[0], ylim[1])
	pyplot.axis('off')
	pyplot.tight_layout()
	pyplot.savefig('{0}_{1}.jpg'.format(iteration, root), bbox_inches='tight', pad_inches = 0, dpi=dpi)
	fig.clf()
	pyplot.close('all')
	gc.collect()
	#pyplot.show()

def readTrees(simulation_target):
	res = {}
	infile = open('{0}_trees.txt'.format(simulation_target), 'r')
	test = 0
	for line in infile:
		if 'positions' in line:
			test = 1
			positions = line.strip().split(': ')[1]
			positions = [ float(i) for i in positions.split(' ') ]
			
			seq = []
		if test!=0 and 'positions' not in line:
			seq.append(line.strip())
	infile.close()
	res['positions'] = positions
	res['sequences'] = seq
	del positions
	del seq
	return(res)


def readMS(simulation_target):
	res = {'sweep':{}, 'neutral':{}}
	for model_tmp in ['sweep', 'neutral']:
		infile = open('{simulation}_{model}.ms'.format(simulation=simulation_target, model=model_tmp), 'r')
		test = 0
		for line in infile:
			if 'positions' in line:
				test = 1
				positions = line.strip().split(': ')[1]
				positions = [ float(i) for i in positions.split(' ') ]
				
				seq = []
			if test!=0 and 'positions' not in line:
				seq.append(line.strip())
		infile.close()
		res[model_tmp]['positions'] = positions
		res[model_tmp]['sequences'] = seq
		del positions
		del seq
	return(res)


def get_distances(positions):
	# computes the distance between 2 SNPs.
	# for SNPs Sa Sb Sc Sd, the function returns a list [abs(Sa-0)+abs(Sa-Sb) ; abs(Sb-Sa)+abs(Sb-Sc); abs(Sc-Sb)+abs(Sc-Sd); abs(Sd-1)+abs(Sd-Sc) ]
	S = len(positions)
	distances = [0]*S
	for pos in range(S):
		if pos==0:
			distances[pos] = abs(positions[pos]-0) + abs(positions[pos]-positions[pos+1])
		else:
			if pos==(S-1):
				distances[pos] = abs(positions[pos]-1) + abs(positions[pos]-positions[pos-1])
			else:
				distances[pos] = abs(positions[pos]-positions[pos-1]) + abs(positions[pos]-positions[pos+1])
	return(distances)
			
	

	
def plotRawData_fullLength(tree, L, simulation_target, colormaps):
	# tree: output of readTrees; a dictionnary containing tree['positions'] and tree['sequences']
	# L: length of the simulated chromosome
	# simulation_target: ID of the simulation {simulation_target]_sumStats.txt for instance
	sequences = tree['sequences']
	positions = tree['positions']
	positions = [ int(round(i*L,0)) for i in positions ]
	nSam = len(sequences)
	data = np.zeros((nSam, L), dtype=int)
	for sam in range(nSam):
		for pos in range(L):
			if pos in positions:
				data[sam][pos]=sequences[sam][positions.index(pos)]
	
	vmin=0
	vmax=1
	n = len(colormaps)
	fig, axs = pyplot.subplots(1, n, figsize=(n * 2 + 2, 3), constrained_layout=True, squeeze=False)
	for [ax, cmap] in zip(axs.flat, colormaps):
		psm = ax.pcolormesh(data, cmap=cmap, rasterized=True, vmin=vmin, vmax=vmax)
#		fig.colorbar(psm, ax=ax)
	pyplot.axis("off")
	pyplot.savefig('{0}_rawData.jpg'.format(simulation_target), bbox_inches='tight', pad_inches = 0, dpi=dpi)
#	pyplot.show()
	pyplot.close()
	gc.collect()
	
	
def plotRawData_onlySNPs(tree, L, simulation_target, distances_target, colormaps, model):
	# tree: output of readTrees; a dictionnary containing tree['positions'] and tree['sequences']
	# L: length of the simulated chromosome
	# simulation_target: ID of the simulation {simulation_target]_sumStats.txt for instance
	sequences = tree['sequences']
	positions = tree['positions']
	nSam = len(sequences)
	data = np.zeros((nSam+2, len(positions)), dtype=float)
	for sam in range(nSam):
		for pos in range(len(positions)):
			data[sam+1][pos]=sequences[sam][pos]

	for pos in range(len(positions)-1):
		data[0][pos] = distances_target[pos]
		data[nSam+1][pos] = distances_target[pos]
	
	vmin=0
	vmax=1
	n = len(colormaps)
	fig, axs = pyplot.subplots(1, n, figsize=(n * 2 + 2, 3), constrained_layout=True, squeeze=False)
	for [ax, cmap] in zip(axs.flat, colormaps):
		psm = ax.pcolormesh(data, cmap=cmap, rasterized=True, vmin=vmin, vmax=vmax)
#		fig.colorbar(psm, ax=ax)
	pyplot.axis("off")
	pyplot.savefig('{iteration}_{model}_rawData.jpg'.format(iteration=simulation_target, model=model), bbox_inches='tight', pad_inches = 0, dpi=dpi)
#	pyplot.show()
#	pyplot.close()
	gc.collect()
	im = cv2.imread("{iteration}_{model}_rawData.jpg".format(iteration=simulation_target, model=model))
	res={}
	res['width'] = im.shape[1]
	res['height'] = im.shape[0]
	return(res)


def readStats(simulations, L):
	res = {'sweep':{}, 'neutral':{}}
	distances = {'sweep':{}, 'neutral':{}}
	tree = {}
	for iteration_tmp in range(len(simulations)):
		for model_tmp in ['sweep', 'neutral']:
			# summary statistics
			infile = open('{simulation}_{model}_sumStats.txt'.format(simulation=simulations[iteration_tmp], model=model_tmp), 'r')
			line_tmp = infile.readline().strip().split('\t')
			
			# get the header
			if iteration_tmp==0:
				header = {}
				for i in range(len(line_tmp)):
					res[model_tmp][i] = [nan]*len(simulations)
					header[i] = line_tmp[i]
			
			# get the stat values
			line_tmp = infile.readline().strip().split('\t')
			for i in range(len(line_tmp)):
				res[model_tmp][i][iteration_tmp] = float(line_tmp[i])
			infile.close()
			
			# positions of SNPs
#			tree = readTrees(simulation_target=simulations[iteration_tmp])

		tree[simulations[iteration_tmp]] = readMS(simulation_target=simulations[iteration_tmp])
		for model_tmp in ['sweep', 'neutral']:
#			distances[model_tmp][iteration_tmp] = get_distances(tree_tmp[model_tmp]['positions'])
			distances[model_tmp][iteration_tmp] = [ tree[simulations[iteration_tmp]][model_tmp]['positions'][i]-tree[simulations[iteration_tmp]][model_tmp]['positions'][i-1] for i in range(1, len(tree[simulations[iteration_tmp]][model_tmp]['positions']), 1) ]
	
	output = {}
	output['header'] = header
	output['stats'] = res
	output['distances'] = distances
	output['MS'] = tree
	del res
	del header
	return(output)


def getSelectedPosition(simulations, L, simulation_target):
	res = {}
	for iteration_tmp in range(len(simulations)):
		if simulations[iteration_tmp]==str(simulation_target):
			infile = open('{0}_parameters.txt'.format(simulations[iteration_tmp]), 'r')
			for line in infile:
				if 'iteration_tmp' not in line:
					line = line.strip().split('\t')
					res[simulations[iteration_tmp]] = float(line[6])/L
			infile.close()
	return(res)

def getParameters(simulations):
	# simulations = [ list of IDs of all simulations within the datapath/ ]
	list_non_numerical_parameters = ['sim_id', 'outcome']
	list_non_numerical_values = ['NF']

	res = {}
	
	for iteration_tmp in simulations:
		res[iteration_tmp] = {}
		infile = open('{0}_sweep_parameters.txt'.format(iteration_tmp), 'r')
		
		line = infile.readline()
		header = {}
		cnt = 0
		for header_tmp in line.strip().split('\t'):
			header[cnt] = header_tmp
			cnt += 1
		
		line = infile.readline()
		line = line.strip().split('\t')
		for i in range(cnt):
			if header[i] not in list_non_numerical_parameters and line[i] not in list_non_numerical_values:
				value_tmp = float(line[i])
			else:
				value_tmp = line[i]
			res[iteration_tmp][header[i]] = value_tmp

		infile.close()
	
	return(res)
	

def getBoundaries(selected_pos, all_pos, pi_obs, pi_exp):
	# selected_pos = target # in [0, 1]
	# all_pos = x_tmp_pi # list of pos # in [0, 1]
	# pi_obs = y_tmp_pi # list of obs pi # in [0, 1]
	# pi_exp = 4.N.mu # 0.002 p/ex
	
	# get the window the closest to the target of selection
	closest_pos = [ abs(i-selected_pos) for i in all_pos ]
	closest_pos = [ i for i in range(len(closest_pos)) if closest_pos[i]==nanmin(closest_pos) ][0] # position of the window which is the closest to the selected target
	
	# get the x_min value
	x_min = 0
	for pos in range(closest_pos, 0, -1):
		if pi_obs[pos]>=pi_exp:
			x_min = pos
			break
		else:
			x_min = pos

	# get the x_max value
	x_max = len(all_pos)
	for pos in range(closest_pos, len(all_pos), 1):
		if pi_obs[pos]>=pi_exp:
			x_max = pos
			break
		else:
			x_max = pos
	res = {}
	
	if isnan(all_pos[x_min])==False and isnan(all_pos[x_max])==False:
#		res['x_min'] = int(round(all_pos[x_min] * width, 0))
#		res['x_max'] = int(round(all_pos[x_max] * width, 0))
		res['x_min'] = all_pos[x_min]
		res['x_max'] = all_pos[x_max]
	else:
		res['x_min'] = nan
		res['x_max'] = nan
	return(res)


##### END OF FUNCTIONS ######

##### START TREATMENT OF DATA #####
if getAllData == 1:
	simulations = [ int(i.split('_')[0]) for i in os.listdir('./') if '.ms' in i ]
	simulations.sort()
	simulations = [ str(i) for i in simulations ]
	simulations = list(set(simulations))
else:
	simulations = [simulation_target]

param = getParameters(simulations)

L = param[simulation_target]['L']
Ne = param[simulation_target]['NeA']
mu = param[simulation_target]['mu']

# get positions: positions are produced by mscalc, corresponds to the mid SNP within each window
positions = {'sweep':[], 'neutral':[]}
for model_tmp in ['sweep', 'neutral']:
	infile = open('{simulation}_{model}_positions.txt'.format(simulation=simulation_target, model=model_tmp), 'r')
	tmp = infile.readline()
	tmp = infile.readline().strip().split('\t')
	positions[model_tmp] = [ float(i) for i in tmp ]
	infile.close()

if object_to_recognize == 1:
#	selected_positions = getSelectedPosition(simulations=simulations, L=L, simulation_target=simulation_target)
	selected_positions = param[simulation_target]['mut_pos']

# get stats
stats = readStats(simulations, L)

# get all distances
distances_all = [ val for model in stats['distances'] for i in stats['distances'][model] for val in stats['distances'][model][i] ]

# positions of statistics in the table
pi = [ i for i in stats['header'].keys() if 'pi_avg' in stats['header'][i] ]
pistd = [ i for i in stats['header'].keys() if 'pi_std' in stats['header'][i] ]
theta = [ i for i in stats['header'].keys() if 'thetaW_avg' in stats['header'][i] ]
tajD = [ i for i in stats['header'].keys() if 'tajD' in stats['header'][i] ]
achazY = [ i for i in stats['header'].keys() if 'achazY' in stats['header'][i] ]
pearson_pi_r = [ i for i in stats['header'].keys() if 'pearson_r' in stats['header'][i] ]
pearson_pi_pval = [ i for i in stats['header'].keys() if 'pearson_pval' in stats['header'][i] ]
nHaplo = [ i for i in stats['header'].keys() if 'nHaplo' in stats['header'][i] ]
H1 = [ i for i in stats['header'].keys() if stats['header'][i][0:3]=='H1_' ]
H2 = [ i for i in stats['header'].keys() if stats['header'][i][0:3]=='H2_' ]
H12 = [ i for i in stats['header'].keys() if stats['header'][i][0:4]=='H12_' ]
H2overH1 = [ i for i in stats['header'].keys() if 'H2overH1' in stats['header'][i] ]
D = [ i for i in stats['header'].keys() if stats['header'][i][0:2]=='D_' ]
r2 = [ i for i in stats['header'].keys() if 'r2' in stats['header'][i] ]

# get min and max values
# pi_range = [ item for i in pi for item in stats['stats'][i] ]
# pistd_range = [ item for i in pistd for item in stats['stats'][i] ]
# theta_range = [ item for i in theta for item in stats['stats'][i] ]
# tajD_range = [ item for i in tajD for item in stats['stats'][i] ]
# achaz_range = [ item for i in achazY for item in stats['stats'][i] ]
# pearson_r_range = [ item for i in pearson_pi_r for item in stats['stats'][i] ]
# pearson_pval_range = [ item for i in pearson_pi_pval for item in stats['stats'][i] ]

pi_range = [ val for i in pi for model in stats['stats'] for val in stats['stats'][model][i] ]
pistd_range = [ val for i in pistd for model in stats['stats'] for val in stats['stats'][model][i] ]
theta_range = [ val for i in theta for model in stats['stats'] for val in stats['stats'][model][i] ]
tajD_range = [ val for i in tajD for model in stats['stats'] for val in stats['stats'][model][i] ]
achaz_range = [ val for i in achazY for model in stats['stats'] for val in stats['stats'][model][i] ]
pearson_r_range = [ val for i in pearson_pi_r for model in stats['stats'] for val in stats['stats'][model][i] ]
pearson_pval_range = [ val for i in pearson_pi_pval for model in stats['stats'] for val in stats['stats'][model][i] ]
nHaplo_range = [ val for i in nHaplo for model in stats['stats'] for val in stats['stats'][model][i] ]
H1_range = [ val for i in H1 for model in stats['stats'] for val in stats['stats'][model][i] ]
H2_range = [ val for i in H2 for model in stats['stats'] for val in stats['stats'][model][i] ]
H12_range = [ val for i in H12 for model in stats['stats'] for val in stats['stats'][model][i] ]
H2overH1_range = [ val for i in H2overH1 for model in stats['stats'] for val in stats['stats'][model][i] ]
D_range = [ val for i in D for model in stats['stats'] for val in stats['stats'][model][i] ]
r2_range = [ val for i in r2 for model in stats['stats'] for val in stats['stats'][model][i] ]


# range (min, max) of summary statistics over all simulations
if getAllData==1:
	min_pi = nanmin(pi_range)
	max_pi = nanmax(pi_range)
	min_pistd = nanmin(pistd_range)
	max_pistd = nanmax(pistd_range)
	min_theta = nanmin(theta_range)
	max_theta = nanmax(theta_range)
	min_tajD = nanmin(tajD_range)
	max_tajD = nanmax(tajD_range)
	min_achaz = nanmin(achaz_range)
	max_achaz = nanmax(achaz_range)
	min_pearsonR = nanmin(pearson_r_range)
	max_pearsonR = nanmax(pearson_r_range)
	min_pearsonP = nanmin(pearson_pval_range)
	max_pearsonP = nanmax(pearson_pval_range)
	min_nHaplo = nanmin(nHaplo_range)
	max_nHaplo = nanmax(nHaplo_range)
	min_H1 = nanmin(H1_range)
	max_H1 = nanmax(H1_range)
	min_H2 = nanmin(H2_range)
	max_H2 = nanmax(H2_range)
	min_H12 = nanmin(H12_range)
	max_H12 = nanmax(H12_range)
	min_H2overH1 = nanmin(H2overH1_range)
	max_H2overH1 = nanmax(H2overH1_range)
	min_D = nanmin(D_range)
	max_D = nanmax(D_range)
	min_r2 = nanmin(r2_range)
	max_r2 = nanmax(r2_range)
	
	outfile_range = open('range_stats.txt', 'w')
	outfile_range.write('min_pi: {stat}\n'.format(stat=min_pi))
	outfile_range.write('min_pistd: {stat}\n'.format(stat=min_pistd))
	outfile_range.write('min_theta: {stat}\n'.format(stat=min_theta))
	outfile_range.write('min_tajD: {stat}\n'.format(stat=min_tajD))
	outfile_range.write('min_achaz: {stat}\n'.format(stat=min_achaz))
	outfile_range.write('min_pearsonR: {stat}\n'.format(stat=min_pearsonR))
	outfile_range.write('min_pearsonP: {stat}\n'.format(stat=min_pearsonP))
	outfile_range.write('min_nHaplo: {stat}\n'.format(stat=min_nHaplo))
	outfile_range.write('min_H1: {stat}\n'.format(stat=min_H1))
	outfile_range.write('min_H2: {stat}\n'.format(stat=min_H2))
	outfile_range.write('min_H12: {stat}\n'.format(stat=min_H12))
	outfile_range.write('min_H2overH1: {stat}\n'.format(stat=min_H2overH1))
	outfile_range.write('min_D: {stat}\n'.format(stat=min_D))
	outfile_range.write('min_r2: {stat}\n'.format(stat=min_r2))
	outfile_range.write('max_pi: {stat}\n'.format(stat=max_pi))
	outfile_range.write('max_pistd: {stat}\n'.format(stat=max_pistd))
	outfile_range.write('max_theta: {stat}\n'.format(stat=max_theta))
	outfile_range.write('max_tajD: {stat}\n'.format(stat=max_tajD))
	outfile_range.write('max_achaz: {stat}\n'.format(stat=max_achaz))
	outfile_range.write('max_pearsonR: {stat}\n'.format(stat=max_pearsonR))
	outfile_range.write('max_pearsonP: {stat}\n'.format(stat=max_pearsonP))
	outfile_range.write('max_nHaplo: {stat}\n'.format(stat=max_nHaplo))
	outfile_range.write('max_H1: {stat}\n'.format(stat=max_H1))
	outfile_range.write('max_H2: {stat}\n'.format(stat=max_H2))
	outfile_range.write('max_H12: {stat}\n'.format(stat=max_H12))
	outfile_range.write('max_H2overH1: {stat}\n'.format(stat=max_H2overH1))
	outfile_range.write('max_D: {stat}\n'.format(stat=max_D))
	outfile_range.write('max_r2: {stat}\n'.format(stat=max_r2))
	outfile_range.close()
else:
	ranges = getRanges()
	min_pi = ranges['min_pi']
	min_pistd = ranges['min_pistd']
	min_theta = ranges['min_theta']
	min_tajD = ranges['min_tajD']
	min_achaz = ranges['min_achaz']
	min_pearsonR = ranges['min_pearsonR']
	min_pearsonP = ranges['min_pearsonP']
	min_nHaplo = ranges['min_nHaplo']
	min_H1 = ranges['min_H1']
	min_H2 = ranges['min_H2']
	min_H12 = ranges['min_H12']
	min_H2overH1 = ranges['min_H2overH1']
	min_D = ranges['min_D']
	min_r2 = ranges['min_r2']
	max_pi = ranges['max_pi']
	max_pistd = ranges['max_pistd']
	max_theta = ranges['max_theta']
	max_tajD = ranges['max_tajD']
	max_achaz = ranges['max_achaz']
	max_pearsonR = ranges['max_pearsonR']
	max_pearsonP = ranges['max_pearsonP']
	max_nHaplo = ranges['max_nHaplo']
	max_H1 = ranges['max_H1']
	max_H2 = ranges['max_H2']
	max_H12 = ranges['max_H12']
	max_H2overH1 = ranges['max_H2overH1']
	max_D = ranges['max_D']
	max_r2 = ranges['max_r2']

# range (min, max) of distances between 3 SNPs over all simulations
min_distance = nanmin(distances_all)
max_distance = nanmax(distances_all)


# iteration : corresponds to the iterration [0, 1, 2, 3, 4] associated to one simulation_target, for instance,
# among ['123', '124', '125', '130', '131'] if runs 123, 124, 125, 130 and 131 are the one that have been performed
iteration = [ i for i in range(len(simulations)) if simulations[i]==simulation_target ][0]

# plot stats
#x_tmp = positions['sweep'][iteration] # positions

#########
# PLOTS #
#########
distances_target = {}
models = ['sweep', 'neutral']

for model in models:
	x_tmp = positions[model]

	if plotStats == 1:
		# pi
		y_tmp_pi = [ stats['stats'][model][i][iteration] for i in pi ] # values
		plot_genome(x=x_tmp, y=y_tmp_pi, xlim=[0,1], ylim=[min_pi, max_pi], iteration=simulations[iteration], root='{model}_pi'.format(model=model))

		# pi std
		y_tmp_pistd = [ stats['stats'][model][i][iteration] for i in pistd ] # values
		plot_genome(x=x_tmp, y=y_tmp_pistd, xlim=[0,1], ylim=[min_pistd, max_pistd], iteration=simulations[iteration], root='{model}_pistd'.format(model=model))

		# theta
		y_tmp_theta = [ stats['stats'][model][i][iteration] for i in theta ] # values
		plot_genome(x=x_tmp, y=y_tmp_theta, xlim=[0,1], ylim=[min_theta, max_theta], iteration=simulations[iteration], root='{model}_theta'.format(model=model))

		# Tajima's D
		y_tmp_tajD = [ stats['stats'][model][i][iteration] for i in tajD ] # values
		plot_genome(x=x_tmp, y=y_tmp_tajD, xlim=[0,1], ylim=[min_tajD, max_tajD], iteration=simulations[iteration], root='{model}_tajimaD'.format(model=model))

		# Achaz
		y_tmp_achazY = [ stats['stats'][model][i][iteration] for i in achazY] # values
		plot_genome(x=x_tmp, y=y_tmp_achazY, xlim=[0,1], ylim=[min_achaz, max_achaz], iteration=simulations[iteration], root='{model}_achaz'.format(model=model))

		# Pearson R
		y_tmp_pearsonR = [ stats['stats'][model][i][iteration] for i in pearson_pi_r] # values
		plot_genome(x=x_tmp, y=y_tmp_pearsonR, xlim=[0,1], ylim=[min_pearsonR, max_pearsonR], iteration=simulations[iteration], root='{model}_pearsonR'.format(model=model))

		# Pearson P
		y_tmp_pearsonP = [ stats['stats'][model][i][iteration] for i in pearson_pi_pval] # values
		plot_genome(x=x_tmp, y=y_tmp_pearsonP, xlim=[0,1], ylim=[min_pearsonP, max_pearsonP], iteration=simulations[iteration], root='{model}_pearsonP'.format(model=model))

		# nHaplo
		y_tmp_nHaplo = [ stats['stats'][model][i][iteration] for i in nHaplo] # values
		plot_genome(x=x_tmp, y=y_tmp_nHaplo, xlim=[0,1], ylim=[min_nHaplo, max_nHaplo], iteration=simulations[iteration], root='{model}_nHaplo'.format(model=model))

		# H1
		y_tmp_H1 = [ stats['stats'][model][i][iteration] for i in H1] # values
		plot_genome(x=x_tmp, y=y_tmp_H1, xlim=[0,1], ylim=[min_H1, max_H1], iteration=simulations[iteration], root='{model}_H1'.format(model=model))

		# H2
		y_tmp_H2 = [ stats['stats'][model][i][iteration] for i in H2] # values
		plot_genome(x=x_tmp, y=y_tmp_H2, xlim=[0,1], ylim=[min_H2, max_H2], iteration=simulations[iteration], root='{model}_H2'.format(model=model))

		# H12
		y_tmp_H12 = [ stats['stats'][model][i][iteration] for i in H12] # values
		plot_genome(x=x_tmp, y=y_tmp_H12, xlim=[0,1], ylim=[min_H12, max_H12], iteration=simulations[iteration], root='{model}_H12'.format(model=model))

		# H2overH1
		y_tmp_H2overH1 = [ stats['stats'][model][i][iteration] for i in H2overH1] # values
		plot_genome(x=x_tmp, y=y_tmp_H2overH1, xlim=[0,1], ylim=[min_H2overH1, max_H2overH1], iteration=simulations[iteration], root='{model}_H2overH1'.format(model=model))

		# D
		y_tmp_D = [ stats['stats'][model][i][iteration] for i in D] # values
		plot_genome(x=x_tmp, y=y_tmp_D, xlim=[0,1], ylim=[min_D, max_D], iteration=simulations[iteration], root='{model}_D'.format(model=model))

		# r2
		y_tmp_r2 = [ stats['stats'][model][i][iteration] for i in r2] # values
		plot_genome(x=x_tmp, y=y_tmp_r2, xlim=[0,1], ylim=[min_r2, max_r2], iteration=simulations[iteration], root='{model}_r2'.format(model=model))
	else:
		y_tmp_pi = [ stats['stats'][model][i][iteration] for i in pi ] # values
	
	# observed distances between SNPs
	distances_target[model] = [ (dst-min_distance)/(max_distance-min_distance) for dst in stats['distances'][model][iteration] ]


	# GLOBAL PICTURE
	if phasing == 0: # if unphased data
	#	data = np.zeros((7, len(positions[iteration])), dtype=float)
		data = np.zeros((7, len(x_tmp)), dtype=float)
	else: # if phased data
		data = np.zeros((14, len(x_tmp)), dtype=float)


	data[0] = [ (stats['stats'][model][i][iteration]-min_pi)/(max_pi-min_pi) for i in pi ] # pi
	data[1] = [ (stats['stats'][model][i][iteration]-min_pistd)/(max_pistd-min_pistd) for i in pistd ] # pi
	data[2] = [ (stats['stats'][model][i][iteration]-min_theta)/(max_theta-min_theta) for i in theta ] # theta
	data[3] = [ (stats['stats'][model][i][iteration]-min_tajD)/(max_tajD-min_tajD) for i in tajD ] # tajD
	data[4] = [ (stats['stats'][model][i][iteration]-min_achaz)/(max_achaz-min_achaz) for i in achazY ] # achaz
	data[5] = [ (stats['stats'][model][i][iteration]-min_pearsonR)/(max_pearsonR-min_pearsonR) for i in pearson_pi_r ] # achaz
	data[6] = [ (stats['stats'][model][i][iteration]-min_pearsonP)/(max_pearsonP-min_pearsonP) for i in pearson_pi_pval ] # achaz

	if phasing == 1:
		data[7] = [ (stats['stats'][model][i][iteration]-min_nHaplo)/(max_nHaplo-min_nHaplo) for i in nHaplo ] # achaz
		data[8] = [ (stats['stats'][model][i][iteration]-min_H1)/(max_H1-min_H1) for i in H1 ] # achaz
		data[9] = [ (stats['stats'][model][i][iteration]-min_H2)/(max_H2-min_H2) for i in H2 ] # achaz
		data[10] = [ (stats['stats'][model][i][iteration]-min_H12)/(max_H12-min_H12) for i in H12 ] # achaz
		data[11] = [ (stats['stats'][model][i][iteration]-min_H2overH1)/(max_H2overH1-min_H2overH1) for i in H2overH1 ] # achaz
		data[12] = [ (stats['stats'][model][i][iteration]-min_D)/(max_D-min_D) for i in D ] # achaz
		data[13] = [ (stats['stats'][model][i][iteration]-min_r2)/(max_r2-min_r2) for i in r2 ] # achaz

	dimensions_globalPic = plot_globalPic(data, [binary], 0, 1, simulations[iteration], model=model) # [width; height] in pixels
	del data

# get coordinates in pixel
for model_tmp in ['sweep', 'neutral']:
	x_tmp = positions[model_tmp]
	
#	if object_to_recognize==1:
	if model_tmp == 'sweep':
		# if selected target
	#	target = selected_positions[simulations[iteration]] # position of the selected target in 0,1
		target = selected_positions/L
		if theta==1:
			pi_exp = 4*Ne*mu*scalar
		else:
			pi_exp = nanmean([ stats['stats']['neutral'][i][iteration] for i in pi ])

		y_tmp_pi = [ stats['stats'][model_tmp][i][iteration] for i in pi ] # values

		xmin_xmax = getBoundaries(selected_pos = target, all_pos = x_tmp, pi_obs = y_tmp_pi, pi_exp = pi_exp)
		center_global = target
		width_global = max([abs(center_global-xmin_xmax['x_min']), abs(center_global-xmin_xmax['x_max'])])*2 # center_global=center_global; width=twice the distance between the center_global and the most distant boundary (x_min or x_max)
		if center_global <0.5:
			if center_global-width_global/2.0 < 0:
				width_global = 2.0*center_global - epsilon
		else:
			if center_global+width_global/2.0 > 1:
				width_global = (1.0-center_global)*2 - epsilon

		output_path = '{object} {x_center} {y_center} {width} {height}\n'.format(x_center=center_global, y_center=0.5, width=width_global, height=1.0, object=1) # object = 1 means the object is a "sweep"
		outfile = open('{0}/{1}_sweep_globalPic.txt'.format(datapath, simulation_target), 'w')
		outfile.write(output_path)
		outfile.close()


	else:
		# if neutral (for the moment)
		center_global=0.5
		xmin_xmax = {}
		xmin_xmax['x_min']=0
		xmin_xmax['x_max']=1
		width_global = 1

		output_path = '{object} {x_center} {y_center} {width} {height}\n'.format(x_center=center_global, y_center=0.5, width=width_global, height=1.0, object=0) # object = 0 means the object is neutrality
		outfile = open('{0}/{1}_neutral_globalPic.txt'.format(datapath, simulation_target), 'w')
		outfile.write(output_path)
		outfile.close()

# #output_path = '{object} {x_center} {y_center} {width} {height}\n'.format(x_center=(x_max+x_min)/2.0, y_center=0.5, width=x_max-x_min, height=1.0, object=object_to_recognize)
# output_path = '{object} {x_center} {y_center} {width} {height}\n'.format(x_center=center_global, y_center=0.5, width=width_global, height=1.0, object=object_to_recognize)
# outfile = open('{0}/{1}_train_globalPic.txt'.format(datapath, simulation_target), 'w')
# outfile.write(output_path)
# outfile.close()


# RAW DATA
for model_tmp in ['sweep', 'neutral']:
#	tree = readMS(simulation_target=simulation_target) 
	tree = stats['MS'][simulations[iteration]][model_tmp]
#	dimensions_rawData = plotRawData_onlySNPs(tree=tree, L=L, simulation_target=simulation_target, distances_target=distances_target, colormaps=[binary])
	dimensions_rawData = plotRawData_onlySNPs(tree=tree, L=L, simulation_target=simulation_target, distances_target=distances_target[model_tmp], colormaps=[binary], model=model_tmp)
#	if object_to_recognize == 1:
	if model_tmp == 'sweep':
		x_tmp = positions[model_tmp]
		y_tmp_pi = [ stats['stats'][model_tmp][i][iteration] for i in pi ] # values
		if theta==1:
			pi_exp = 4*Ne*mu*scalar
		else:
			pi_exp = nanmean([ stats['stats']['neutral'][i][iteration] for i in pi ])
		
		xmin_xmax = getBoundaries(selected_pos = target, all_pos = x_tmp, pi_obs = y_tmp_pi, pi_exp = pi_exp)
		center_global = selected_positions/L
		center_raw = getSNP(center_global, tree['positions'])
		xmin_raw = getSNP(xmin_xmax['x_min'], tree['positions'])
		xmax_raw = getSNP(xmin_xmax['x_max'], tree['positions'])
		width_raw = max([abs(center_raw-xmin_raw), abs(center_raw-xmax_raw)])*2 
		if center_raw <0.5:
			if center_raw-width_raw/2.0 < 0:
				width_raw = 2.0*center_raw - epsilon
		else:
			if center_raw+width_raw/2.0 > 1:
				width_raw = (1.0-center_raw)*2 - epsilon
		output_path = '{object} {x_center} {y_center} {width} {height}\n'.format(x_center=center_raw, y_center=0.5, width=width_raw, height=1.0, object=1) # object=1 for sweep, 0 for neutral
	else:
		# if neutral (for the moment)
		center_raw = 0.5
		xmin_xmax = {}
		xmin_xmax['x_min'] = 0
		xmin_xmax['x_max'] = 1
		width_raw = 1
#	output_path = '{object} {x_center} {y_center} {width} {height}\n'.format(x_center=center_raw, y_center=0.5, width=width_raw, height=1.0, object=object_to_recognize)
		output_path = '{object} {x_center} {y_center} {width} {height}\n'.format(x_center=center_raw, y_center=0.5, width=width_raw, height=1.0, object=0) # object=1 for sweep, 0 for neutral
	#output_path = '{object} {x_min} {y_min} {x_max} {y_max}\n'.format(path=datapath, simulation=simulations[iteration], x_min=x_min, y_min=y_min, x_max=x_max, y_max=y_max, object=object_to_recognize)
	outfile = open('{datapath}/{iteration}_{model}_rawData.txt'.format(datapath=datapath, iteration=simulation_target, model=model_tmp), 'w')
	outfile.write(output_path)
	outfile.close()

