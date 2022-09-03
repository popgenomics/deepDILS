from matplotlib import pyplot
from matplotlib import cm
import os
import numpy as np
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

help_message = '\tpython scripts that works after the SLiM simulations to produce PNG files representing the data along the simulated chromosome\n'
help_message += '\tto run this script for instance on the simulation producing the following files: 123_parameters.txt 123_positions.txt 123_sumStats.txt 123_trees.txt\n'
help_message += '\t\tpython3 sim2box_single.py dpi=300 datapath=/home/croux/Programmes/yolo_box/simulations simulation=123 object=posSelection\n'
binary = cm.get_cmap('binary', 512)

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
def getSNP(snp, positions):
	# function usefull for the rawdata
	# it 1. gets the position of 'snp' within the vector 'positions'; 2. converts it in a float in [0,1]
	# snp is a float readen from the ms output (positions) (in [0,1])
	# positions is a vector of positions (all in [0,1])
	distances = [ abs(i-snp) for i in positions ]
	closest_snp = distances.index(min(distances))
	
	return(closest_snp / (1.0*len(positions)))

def plot_globalPic(data, colormaps, vmin, vmax, iteration):
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
	pyplot.savefig("{0}_globalPic.png".format(iteration), bbox_inches='tight', pad_inches = 0, dpi=dpi)
	pyplot.close()
	gc.collect()
	im = cv2.imread("{0}_globalPic.png".format(iteration))
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
	pyplot.savefig('{0}_{1}.png'.format(iteration, root), bbox_inches='tight', pad_inches = 0, dpi=dpi)
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


def get_distances(positions, L):
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
	pyplot.savefig('{0}_rawData.png'.format(simulation_target), bbox_inches='tight', pad_inches = 0, dpi=dpi)
#	pyplot.show()
	pyplot.close()
	gc.collect()
	
	
def plotRawData_onlySNPs(tree, L, simulation_target, distances_target, colormaps):
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

	for pos in range(len(positions)):
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
	pyplot.savefig('{0}_rawData.png'.format(simulation_target), bbox_inches='tight', pad_inches = 0, dpi=dpi)
#	pyplot.show()
	pyplot.close()
	gc.collect()
	im = cv2.imread("{0}_rawData.png".format(simulation_target))
	res={}
	res['width'] = im.shape[1]
	res['height'] = im.shape[0]
	return(res)


def readStats(simulations, L):
	res = {}
	distances = {}
	for iteration in range(len(simulations)):
		# summary statistics
		infile = open('{0}_sumStats.txt'.format(simulations[iteration]), 'r')
		tmp = infile.readline().strip().split('\t')
		
		if iteration==0:
			header = {}
			for i in range(len(tmp)):
				res[i] = [nan]*len(simulations)
				header[i] = tmp[i]
		
		tmp = infile.readline().strip().split('\t')
		for i in range(len(tmp)):
			res[i][iteration] = float(tmp[i])
		
		infile.close()
		
		# positions of SNPs
		tree = readTrees(simulation_target=simulations[iteration])
		distances[iteration] = get_distances(tree['positions'], L)
	
	output = {}
	output['header'] = header
	output['stats'] = res
	output['distances'] = distances
	del res
	del header
	return(output)


def getSelectedPosition(simulations, L, simulation_target):
	res = {}
	for iteration in range(len(simulations)):
		if simulations[iteration]==str(simulation_target):
			infile = open('{0}_parameters.txt'.format(simulations[iteration]), 'r')
			for line in infile:
				if 'iteration' not in line:
					line = line.strip().split('\t')
					res[simulations[iteration]] = float(line[6])/L
			infile.close()
	return(res)

def getParameters(simulations):
	# simulations = [ list of IDs of all simulations within the datapath/ ]
	res = {}
	iteration=0 # assumes that all simulations share the same Ne, L and mu
	infile = open('{0}_parameters.txt'.format(simulations[iteration]), 'r')
	for line in infile:
		if 'iteration' not in line:
			line = line.strip().split('\t')
			res['L'] = float(line[5])
			res['Ne'] = float(line[1])
			res['mu'] = float(line[4])
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
simulations = [ int(i.split('_')[0]) for i in os.listdir('./') if 'trees.txt' in i ]
simulations.sort()
simulations = [ str(i) for i in simulations ]

param = getParameters(simulations)
L = param['L']
Ne = param['Ne']
mu = param['mu']

# get positions: positions are produced by mscalc, corresponds to the mid SNP within each window
positions = []
for iteration in range(len(simulations)):
	infile = open('{0}_positions.txt'.format(simulations[iteration]), 'r')
	tmp = infile.readline()
	tmp = infile.readline().strip().split('\t')
	positions.append( [ float(i) for i in tmp ] )
	infile.close()

if object_to_recognize == 1:
	selected_positions = getSelectedPosition(simulations=simulations, L=L, simulation_target=simulation_target)

# get stats
stats = readStats(simulations, L)
distances_all = [ j for i in stats['distances'] for j in stats['distances'][i] ]

# positions of statistics in the table
pi = [ i for i in stats['header'].keys() if 'pi_avg' in stats['header'][i] ]
pistd = [ i for i in stats['header'].keys() if 'pi_std' in stats['header'][i] ]
theta = [ i for i in stats['header'].keys() if 'thetaW_avg' in stats['header'][i] ]
tajD = [ i for i in stats['header'].keys() if 'tajD' in stats['header'][i] ]
achazY = [ i for i in stats['header'].keys() if 'achazY' in stats['header'][i] ]
pearson_pi_r = [ i for i in stats['header'].keys() if 'pearson_r' in stats['header'][i] ]
pearson_pi_pval = [ i for i in stats['header'].keys() if 'pearson_pval' in stats['header'][i] ]


# get min and max values
pi_range = [ item for i in pi for item in stats['stats'][i] ]
pistd_range = [ item for i in pistd for item in stats['stats'][i] ]
theta_range = [ item for i in theta for item in stats['stats'][i] ]
tajD_range = [ item for i in tajD for item in stats['stats'][i] ]
achaz_range = [ item for i in achazY for item in stats['stats'][i] ]
pearson_r_range = [ item for i in pearson_pi_r for item in stats['stats'][i] ]
pearson_pval_range = [ item for i in pearson_pi_pval for item in stats['stats'][i] ]

# range (min, max) of summary statistics over all simulations
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

# range (min, max) of distances between 3 SNPs over all simulations
min_distance = nanmin(distances_all)
max_distance = nanmax(distances_all)

# plot stats
x_tmp = positions[iteration] # positions

# iteration : corresponds to the iterration [0, 1, 2, 3, 4] associated to one simulation_target, for instance,
# among ['123', '124', '125', '130', '131'] if runs 123, 124, 125, 130 and 131 are the one that have been performed
iteration = [ i for i in range(len(simulations)) if simulations[i]==simulation_target ][0]

# pi
y_tmp_pi = [ stats['stats'][i][iteration] for i in pi ] # values
plot_genome(x=x_tmp, y=y_tmp_pi, xlim=[0,1], ylim=[min_pi, max_pi], iteration=simulations[iteration], root='pi')

# pi std
y_tmp_pistd = [ stats['stats'][i][iteration] for i in pistd ] # values
plot_genome(x=x_tmp, y=y_tmp_pistd, xlim=[0,1], ylim=[min_pistd, max_pistd], iteration=simulations[iteration], root='pistd')

# theta
y_tmp_theta = [ stats['stats'][i][iteration] for i in theta ] # values
plot_genome(x=x_tmp, y=y_tmp_theta, xlim=[0,1], ylim=[min_theta, max_theta], iteration=simulations[iteration], root='theta')

# Tajima's D
y_tmp_tajD = [ stats['stats'][i][iteration] for i in tajD ] # values
plot_genome(x=x_tmp, y=y_tmp_tajD, xlim=[0,1], ylim=[min_tajD, max_tajD], iteration=simulations[iteration], root='tajimaD')

# Achaz
y_tmp_achazY = [ stats['stats'][i][iteration] for i in achazY] # values
plot_genome(x=x_tmp, y=y_tmp_achazY, xlim=[0,1], ylim=[min_achaz, max_achaz], iteration=simulations[iteration], root='achaz')

# Pearson R
y_tmp_pearsonR = [ stats['stats'][i][iteration] for i in pearson_pi_r] # values
plot_genome(x=x_tmp, y=y_tmp_pearsonR, xlim=[0,1], ylim=[min_pearsonR, max_pearsonR], iteration=simulations[iteration], root='pearsonR')

# Pearson P
y_tmp_pearsonP = [ stats['stats'][i][iteration] for i in pearson_pi_pval] # values
plot_genome(x=x_tmp, y=y_tmp_pearsonP, xlim=[0,1], ylim=[min_pearsonP, max_pearsonP], iteration=simulations[iteration], root='pearsonP')

# observed distances between SNPs
distances_target = [ (dst-min_distance)/(max_distance-min_distance) for dst in stats['distances'][iteration] ]


# GLOBAL PICTURE
data = np.zeros((7, len(positions[iteration])), dtype=float)
data[0] = [ (stats['stats'][i][iteration]-min_pi)/(max_pi-min_pi) for i in pi ] # pi
data[1] = [ (stats['stats'][i][iteration]-min_pistd)/(max_pistd-min_pistd) for i in pistd ] # pi
data[2] = [ (stats['stats'][i][iteration]-min_theta)/(max_theta-min_theta) for i in theta ] # theta
data[3] = [ (stats['stats'][i][iteration]-min_tajD)/(max_tajD-min_tajD) for i in tajD ] # tajD
data[4] = [ (stats['stats'][i][iteration]-min_achaz)/(max_achaz-min_achaz) for i in achazY ] # achaz
data[5] = [ (stats['stats'][i][iteration]-min_pearsonR)/(max_pearsonR-min_pearsonR) for i in pearson_pi_r ] # achaz
data[6] = [ (stats['stats'][i][iteration]-min_pearsonP)/(max_pearsonP-min_pearsonP) for i in pearson_pi_pval ] # achaz

dimensions_globalPic = plot_globalPic(data, [binary], 0, 1, simulations[iteration]) # [width; height] in pixels
del data

# get coordinates in pixel
if object_to_recognize==1:
	# if selected target
	target = selected_positions[simulations[iteration]] # position of the selected target in 0,1
	xmin_xmax = getBoundaries(selected_pos = target, all_pos = x_tmp, pi_obs = y_tmp_pi, pi_exp = 4*Ne*mu*scalar)
	center_global = target
	width_global = max([abs(center_global-xmin_xmax['x_min']), abs(center_global-xmin_xmax['x_max'])])*2 # center_global=center_global; width=twice the distance between the center_global and the most distant boundary (x_min or x_max)
	if center_global <0.5:
		if center_global-width_global/2.0 < 0:
			width_global = 2.0*center_global - epsilon
	else:
		if center_global+width_global/2.0 > 1:
			width_global = (1.0-center_global)*2 - epsilon
else:
	# if neutral (for the moment)
	center_global=0.5
	xmin_xmax = {}
	xmin_xmax['x_min']=0
	xmin_xmax['x_max']=1
	width_global = 1

#output_path = '{object} {x_center} {y_center} {width} {height}\n'.format(x_center=(x_max+x_min)/2.0, y_center=0.5, width=x_max-x_min, height=1.0, object=object_to_recognize)
output_path = '{object} {x_center} {y_center} {width} {height}\n'.format(x_center=center_global, y_center=0.5, width=width_global, height=1.0, object=object_to_recognize)
outfile = open('{0}/{1}_train_globalPic.txt'.format(datapath, simulation_target), 'w')
outfile.write(output_path)
outfile.close()


# RAW DATA
tree = readTrees(simulation_target=simulation_target) 
dimensions_rawData = plotRawData_onlySNPs(tree=tree, L=L, simulation_target=simulation_target, distances_target=distances_target, colormaps=[binary])
if object_to_recognize == 1:
	center_raw = getSNP(center_global, positions[iteration])
	xmin_raw = getSNP(xmin_xmax['x_min'], positions[iteration])
	xmax_raw = getSNP(xmin_xmax['x_max'], positions[iteration])
	width_raw = max([abs(center_raw-xmin_raw), abs(center_raw-xmax_raw)])*2 
	if center_raw <0.5:
		if center_raw-width_raw/2.0 < 0:
			width_raw = 2.0*center_raw - epsilon
	else:
		if center_raw+width_raw/2.0 > 1:
			width_raw = (1.0-center_raw)*2 - epsilon
else:
	# if neutral (for the moment)
	center_raw = 0.5
	xmin_xmax = {}
	xmin_xmax['x_min'] = 0
	xmin_xmax['x_max'] = 1
	width_raw = 1
output_path = '{object} {x_center} {y_center} {width} {height}\n'.format(x_center=center_raw, y_center=0.5, width=width_raw, height=1.0, object=object_to_recognize)
#output_path = '{object} {x_min} {y_min} {x_max} {y_max}\n'.format(path=datapath, simulation=simulations[iteration], x_min=x_min, y_min=y_min, x_max=x_max, y_max=y_max, object=object_to_recognize)
outfile = open('{0}/{1}_train_rawData.txt'.format(datapath, simulation_target), 'w')
outfile.write(output_path)
outfile.close()

