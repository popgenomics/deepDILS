#!/usr/bin/python

# test: msms 10 10 -s tbs -SAA 200 -SaA tbs -SF 1e-2 -Smu 0
# cat prior.txt | msms 50 10 -s tbs -r 1 10000 -SAA tbs -SaA 100 -SF 0.01 -N 100000 >output_test.msms
import sys
import time
from numpy import nan
from numpy import mean
from numpy import nanmean
from numpy import std 
from numpy import sum as somme
from numpy import sqrt
from numpy import power
from scipy.stats import pearsonr
from decimal import Decimal
from collections import Counter
from typing import Dict, List, Any
import argparse

minNumbSNP = 0 # statistics (pi, tajD, pearson's, etc ...) are only computed if a number of SNPs>= minNumbSNP is present within the studied fragment

parser = argparse.ArgumentParser(description='Takes the output of ms(ms) as input. Calculates summary statistics on sliding windows (width+step). Can calculate these statistics over several nRep replicates, for each of the nCombParam parameter combinations. Output files have root name root.')

parser.add_argument('--infile', type=str, default='simulation.msms', help='Name of the ms file (default: %(default)s)')
parser.add_argument('--nIndiv', type=int, default=20, help='Number of simulated sampled gametes (default: %(default)s)')
parser.add_argument('--nCombParam', type=int, default=1, help='Number of combination of parameters in the file (default: %(default)s)')
parser.add_argument('--regionSize', type=int, default=100000, help='Length of the chromosome (default: %(default)s)')
parser.add_argument('--width', type=float, default=0.1, help='Width of the sliding window (default: %(default)s)')
parser.add_argument('--step', type=float, default=0.05, help='Step of the sliding window (default: %(default)s)')
parser.add_argument('--nRep', type=int, default=1, help='Number of replicates of the same combination of parameters (default: %(default)s)')
parser.add_argument('--root', type=str, default='CST_0', help='Name root for outputs (default: %(default)s)')

args = parser.parse_args()

# Access the arguments using their names
inputFileName = args.infile
nIndiv = args.nIndiv
nCombParam = args.nCombParam
regionSize = args.regionSize
width = args.width
step = args.step
nRep = args.nRep
rootOutputFileName = args.root

##################
# don't touch    #
# nRep = 1       #
##################

def stdCustom(liste, moyenne, longueurRegion):
	# 'moyenne' = mean of ['liste' (size L) + longueurRegion x 0]
	res = 0.0
	for i in liste:
		res += power((i - moyenne), 2)
	res /= longueurRegion
	return(sqrt(res))


def getParams(x):
	res = {}
	x = x[3:].split("\t")
	for i in range(len(x)):
		if x[i][0]=="-":
			param = x[i][1:]
			value = float(x[i+1])
			res[param] = value
	return(res)

def parse_msms_v2(file_path: str, nIndiv: int, nCombParam: int, nRep: int) -> Dict[int, Dict[str, Any]]:
	"""
	Parses the msms output file.

	Args:
	- file_path (str): Path to the msms output file
	- nIndiv (int): Number of individuals
	- nCombParam (int): Combined parameter count
	- nRep (int): Repetition count

	Returns:
	- Dict[int, Dict[str, Any]]: Simulated alignments
	"""
	compteur = nCombParam * nRep
	res: Dict[int, Dict[str, Any]] = {}
	
	with open(file_path, "r") as infile:
		nLocus = -1
		nIndTmp = -1
		test = -9999
		
		for line in infile:
			line = line.strip()
			
			if "//" in line:
				test = 1
				nLocus += 1
				res[nLocus] = {}
				continue
			
			if "segsites" in line and nLocus >= 0:
				res[nLocus]["segsites"] = int(line.split(" ")[1])
				if res[nLocus]["segsites"] == 0:
					continue
				continue
			
			if "positions" in line and nLocus >= 0:
				res[nLocus]["positions"] = [float(pos) for pos in line.split(" ")[1:]]
				res[nLocus]["haplotypes"] = []
				nIndTmp = -1
				continue
			
			if line != "" and test == 1 and '\t' not in line and 'Frequency' not in line:
				nIndTmp += 1
				res[nLocus]["haplotypes"].append(line)
				
				if nIndTmp == (nIndiv - 1):
					tmp = comp_differences(res[nLocus]["haplotypes"], res[nLocus]["segsites"])
					res[nLocus]['kxy'] = tmp['kxy_allSNP']
					res[nLocus]['kxy_singletons'] = tmp['kxy_singletons']
	
	return res


def parse_msms(x, nIndiv, nCombParam, nRep):
	# returns the simulated alignments
	# x = msms output file
	compteur = nCombParam * nRep
	res = {}
	infile = open(x, "r")
	nLocus = -1
	nIndTmp = -1
	test = -9999
	for i in infile:
		i = i.strip()
		if "//" in i:
			test = 1
			nLocus += 1
#			if nLocus%(compteur/100) == 0:
#				print("simulation {0} over a total of {1}: {2}".format(nLocus, compteur, time.strftime("%H:%M:%S")))
			res[nLocus] = {}
			continue
		if "segsites" in i and nLocus >= 0:
			res[nLocus]["segsites"] = int(i.split(" ")[1])
			if res[nLocus]["segsites"] == 0:
				continue
			continue
		if "positions" in i and nLocus >= 0:
			res[nLocus]["positions"] = [ float(pos) for pos in i.split(" ")[1::] ]
			res[nLocus]["haplotypes"] = []
			nIndTmp = -1
			continue
		if i != "" and test == 1:
			nIndTmp += 1
#			print("{0}: {1}".format(nIndTmp, i))
			res[nLocus]["haplotypes"].append(i)
			if nIndTmp == (nIndiv-1):
#				print("{0}\n{1}".format(nLocus, "\n".join(res[nLocus]["haplotypes"])))
				tmp = comp_differences(res[nLocus]["haplotypes"], res[nLocus]["segsites"])
#				print(a["kxy_allSNP"])
				res[nLocus]['kxy'] = tmp['kxy_allSNP'] 
				res[nLocus]['kxy_singletons'] = tmp['kxy_singletons'] 
#				del res[nLocus]['haplotypes']	
	infile.close()
	return(res)


def comp_differences(alignement, nSegSites):
	n = len(alignement) # number of individuals in the alignement
	res = {}
	res["kxy_allSNP"] = [] # vector containing pi for different SNPs, the probability to sample 2 different alleles
	res["kxy_singletons"] = [] # vector containing pi for different SNPs for singletons only
	for i in range(nSegSites):
		n0 = 0
		n1 = 0
		for j in alignement:
			if j[i] == '0':
				n0 += 1
			if j[i] == '1':
				n1 += 1
		pi = n0*n1
		res["kxy_allSNP"].append(pi)
		if n0==1 or n1==1: # if only one derived allele or one ancestral allele
			res["kxy_singletons"].append(pi)
	return(res)


def comp_pi(alignement, nSegSites):
	n = len(alignement) # number of individuals in the alignement
	res = {}
	res["pi_allSNP"] = [] # vector containing pi for different SNPs, the probability to sample 2 different alleles
	res["pi_singletons"] = [] # vector containing pi for different SNPs for singletons only
	for i in range(nSegSites):
		n0 = 0
		n1 = 0
		nComparaisons = n*(n-1.0) / 2.0
		for j in alignement:
			if j[i] == '0':
				n0 += 1
			if j[i] == '1':
				n1 += 1
		pi = n0*n1/nComparaisons
		res["pi_allSNP"].append(pi)
		if n0==1 or n1==1: # if only one derived allele or one ancestral allele
			res["pi_singletons"].append(pi)
	res["nSingletons"] = len(res["pi_singletons"])
	res["nSNP"] = len(res["pi_allSNP"])
	return(res)


def window(width, step):
	width=Decimal(width)
	step=Decimal(step)
	bins = {}
	cnt = 0
	minB = 0
	maxB = minB + width
	bins[cnt] = {}
	bins[cnt]["min"] = float(minB)
	bins[cnt]["max"] = float(maxB)
	while((minB+width) <= 1 and maxB < 1):
		cnt += 1
		minB += step
		maxB += step
		if maxB > 1:
			maxB = 1
		bins[cnt] = {}
		bins[cnt]["min"] = float(minB)
		bins[cnt]["max"] = float(maxB)
	return(bins)


def tajimaD(nInd, pi, nS):
	# nInd = number of individuals in the alignment
	# pi = sum(Kxy over SNPs and pairwise comparisons) / number_of_pairwise_comparisons
	# nS = number of SNPs within the alignment
	# a1 and a2
	if nS < minNumbSNP:
#		return([nan, nan])
		return([0, 0])
	else:
		a1, a2 = 0.0, 0.0
		for i in range(nInd-1):
			a1 += 1.0/(i+1.0)
			a2 += 1.0/power(i+1, 2)
		# b1 and b2
		b1 = (nInd + 1.0) / (3.0 * (nInd - 1.0))
		b2 = 2.0 * (power(nInd, 2) + nInd + 3.0) / (9.0*nInd * (nInd - 1.0))
		# c1 and c2
		c1 = b1 - 1.0 / a1
		c2 = b2 - (nInd + 2.0) / (a1 * nInd) + a2/(power(a1, 2))
		# e1 and e2
		e1 = c1/a1
		e2 = c2/(power(a1, 2) + a2)
		# pi is assumed to already be: sum(pi over SNPs)/nCombination, let's compute thetaW
		thetaW = nS / a1
		# denominateur
		denominateur = e1*nS + e2*nS*(nS - 1.0)
		denominateur = sqrt(denominateur)
		# test
	#	print("a1 = {0}\na2 = {1}\nb1 = {2}\nb2 = {3}\nc1 = {4}\nc2 = {5}\ne1 = {6}\ne2 = {7}\npi = {8}\nthetaW = {9}\ndenoM = {10}".format(a1, a2, b1, b2, c1, c2, e1, e2, pi, thetaW, denominateur))
		#tajima D
		if denominateur==0:
#			return([nan, nan])
			return([0, 0])
		else:
			return([(pi - thetaW) / denominateur, thetaW]) # return [tajD, thetaW]


def achazY(n, kxy):
	# n = number of individuals
	# kx = vector of pi values for each singletons over the bin
	nPairwiseComp = n * (n-1.0)/2.0
	singletons = [ i for i in kxy if i==(n-1) ]
	nSingletons = len(kxy)
	if nSingletons < minNumbSNP: # if less than minNumbSNP of SNPs are present in the alignement --> return "nan" value
		return(nan)
#		return(0)
	else:
		an = 0.0
		bn = 0.0
		for i in range(n-1):
			an += 1.0/(i+1)
			bn += 1.0/(i+1)**2
		# f*
		f = (n - 3.0) / (an * (n-1) - n)
		# alpha
		alpha = f**2 * (an - n / (n-1.0)) + f * (an * (4.0*(n+1)/((n-1.0)**2)) - 2 * ((n+3.0)/(n-1))) - an*(8*(n+1.0)/(n*(n-1)**2)) + (n**2 + n + 60.0)/(3*n * (n-1))
		# beta
		beta = f**2 * (bn - (2*n - 1.0)/(n-1.0)**2) + f * (bn * 8.0/(n-1.0) - an * 4.0/(n*(n-1.0)) - (n**3 + 12*n**2 - 35*n + 18.0)/(n*(n-1.0)**2)) - bn * 16.0/(n*(n-1)) + an * 8.0 / (n**2 * (n-1.0)) + 2 * (n**4 + 110.0*n**2 - 255*n + 126.0) / (9.0 * n**2 * (n-1)**2)
		# mu
		mu = an - n/(n-1.0)
		# theta
		theta = (1.0 * nSingletons) / mu
		# pi singleton: sum of kxy for singletons only
		pi_singleton = 0.0
		for i in singletons:
			pi_singleton += i
		pi_singleton /= nPairwiseComp
		# compute Y* from Achaz
		numerateur = pi_singleton - f*nSingletons/an
		denominateur = alpha * theta + beta * theta**2
		if denominateur != 0:
			return((1.0 * numerateur)/(1.0 * denominateur))
		else:
			return(nan)


def LD(haplotypes, positions_bin, width, min_width):
	# haplotypes : list with all chromosomes
	# positions: list of positions of all SNPs (between 0 and 1)
	# width: size of a window (between 0 and 1, in ms unit)
	# min_width: fraction of the width used as a minimum distance between two SNP to compute LD. If width=0.1 and min_width=0.5, so, only pairs of SNP distant at
	# 	a minimum of 0.05 will be used
	min_dist = min_width*width
	nPos = len(positions_bin)
	nIndiv = len(haplotypes)

	D = []
	r_sqr = []
	for pos_a in range(nPos): # loop over positions of SNP a (to compute LD between a and b)
		pos_b = [ i for i in range(nPos) if positions_bin[i] > (positions_bin[pos_a]+min_dist) ]
		
		for b in pos_b:
			n_a1 = 0
			n_b1 = 0
			n_a1_b1 = 0
			for ind in range(nIndiv):
				a1 = haplotypes[ind][ pos_a ]
				b1 = haplotypes[ind][ b ]

				if a1 == '1':
					n_a1 += 1
					
					if b1 =='1':
						n_b1 += 1
						n_a1_b1 += 1
				else:
					if b1 == '1':
						n_b1 += 1
			
			f_a1 = n_a1/nIndiv
			f_b1 = n_b1/nIndiv
			f_a1_b1 = n_a1_b1/nIndiv
			
			D_tmp = f_a1_b1 - f_a1 * f_b1
			D.append(D_tmp)
			
			denom = f_a1 * (1-f_a1) * f_b1 * (1-f_b1)
			if denom==0:
				r_sqr_tmp=nan
			else:
				r_sqr_tmp=power(D_tmp, 2)/denom
			r_sqr.append(r_sqr_tmp)
	
	if D.count(nan)==len(D):
#		D_mean=nan
		D_mean=1
	else:
		D_mean=nanmean(D)

	if r_sqr.count(nan)==len(r_sqr):
#		r_sqr_mean=nan
		r_sqr_mean=1
	else:
		r_sqr_mean=nanmean(r_sqr)
	
	res = {'D':D_mean, 'r2':r_sqr_mean}
	return(res)	

def sort_rows(haplotypes):
	haploCount = Counter(haplotypes)
	res = ''
	for haplo_i in haploCount.most_common():
		haplo_tmp = [haplo_i[0]]
		nHaplo_tmp = haplo_i[1]
		res += '\n'.join(haplo_tmp*nHaplo_tmp) + '\n'
	return(res)

def haploStats(haplotypes, bins, positions, width, min_width):
	# haplotypes : list with all chromosomes
	# bins: dictionnary with coordinates (start/end) of each bin between 0 (first position of the chromosome) and 1 (last position of the chromosome)
	# positions: list of positions of all SNPs (between 0 and 1)
	# width: size of a window (between 0 and 1, in ms unit)
	# min_width: fraction of the width used as a minimum distance between two SNP to compute LD. If width=0.1 and min_width=0.5, so, only pairs of SNP distant at
	# 	a minimum of 0.05 will be used
	res = {}

	for bin_tmp in bins:
		pos_bin_tmp = [ i for i in range(len(positions)) if positions[i]>bins[bin_tmp]['min'] and positions[i]<bins[bin_tmp]['max'] ]
		if len(pos_bin_tmp)<=2:
			res[bin_tmp]={'binID':bin_tmp, 'nHaplo':1, 'H1':1, 'H2':nan, 'H12':1, 'H2overH1':nan, 'D':nan, 'r2':nan}

		else:
			start_tmp = min(pos_bin_tmp)
			end_tmp = max(pos_bin_tmp)
			
			positions_bin = positions[start_tmp:(end_tmp+1)]
			
			haplotype_tmp = []
			for i in haplotypes:
				haplotype_tmp.append(i[start_tmp:(end_tmp+1)])

			cnt = 0
			haploCount = Counter(haplotype_tmp)

			# nHaplo
			## number of haplotypes
			nHaplo = len(haploCount)

			# H1
			## The probability that two haplotypes, taken randomly from the population, are identical.
			if nHaplo <= 1: # if a single fixed haplotype
				H1 = 1
			else:
				H1 = sum([ power((count/nIndiv), 2) for count in haploCount.values() ])

			# H12
			## The haplotype homozygosity, only the two most common haplotypes are treated as one. Useful for looking for genetic sweeps.
			## If there is only one fixed haplotype, H12 = haplotype homozygosity = 1.
			## Cite: Garud, N. R. et al. (2015). Recent Selective Sweeps in North American Drosophila melanogaster Show Signatures of Soft Sweeps. PLOS Genetics 11 (2), pp. 1–32.

			if nHaplo <= 2:
				# if a single fixed haplotype, or only two haplotypes
				H12 = 1
			else:
				top_two_freq = [ count[1]/nIndiv for count in haploCount.most_common(2) ]
				n1 = sum([ power(count/nIndiv, 2) for count in haploCount.values() ])
				p1 = top_two_freq[0]
				p2 = top_two_freq[1]

				H12 = n1 + 2 * p1 * p2

			# H2overH1
			## The haplotype homozygosity ignoring the most common haplotype
			## Cite: Garud, N. R. et al. (2015). Recent Selective Sweeps in North American Drosophila melanogaster Show Signatures of Soft Sweeps. PLOS Genetics 11 (2), pp. 1–32.

			if nHaplo <= 1:
				H2 = nan
				H2overH1 = nan
			else:
				most_frequent = haploCount.most_common(1)[0][1] / nIndiv
				H2 = H1 - power(most_frequent, 2)
				H2overH1 = H2/H1
			
			# LD
			res_LD = LD(haplotypes=haplotype_tmp, positions_bin=positions_bin, width=width, min_width=min_width)
			res[bin_tmp]={'binID':bin_tmp, 'nHaplo':nHaplo, 'H1':H1, 'H2':H2, 'H12':H12, 'H2overH1':H2overH1, 'D':res_LD['D'], 'r2':res_LD['r2']}
	return(res)

def calc_window(rep, nRep, bins, regionSize, nIndiv, totalData, min_width):
	# function used to return list of 'positions' and 'associated pi' within different bins, over replicates in the
	# rep = the id of the surveyed replicated simulations
	# nRep = number of times the simulation had been replicated
	nPairwiseComp = (nIndiv * (nIndiv - 1.0)) / 2.0
	bins_tmp = {} # bins_tmp[bin_Id][replicate_Id][['positions'], ['kxy']]

	resLD = {} # resLD[iteration][bin]['binID', 'nHaplo', 'H1', 'H2', 'H12', 'H2overH1', 'D', 'r2']
	for i in range(nRep):
		resLD[i]=haploStats(haplotypes=totalData[i]['haplotypes'], bins=bins, positions=totalData[i]['positions'], width=width, min_width=min_width)
	
	for i in bins:
		bins_tmp[i] = {}
		for j in range(nRep):
			bins_tmp[i][j] = {} # bins_tmp[bin ID][rep ID]
			bins_tmp[i][j]['positions'] = [] # bins_tmp[bin ID][rep ID]
	#		bins_tmp[i][j]['pi'] = [] # bins_tmp[bin ID][rep ID]
			bins_tmp[i][j]['kxy'] = [] # bins_tmp[bin ID][rep ID]
	for i in range(nRep): # loop over the replicates
		for j in range(len(totalData[rep*nRep + i]['positions'])): # loop over positions
			pos_tmp = totalData[rep*nRep + i]['positions'][j]
			kxy_tmp = totalData[rep*nRep + i]['kxy'][j] # sum of pairwise differences for a SNP
			bin_value = [ k for k in bins if pos_tmp >= bins[k]['min'] and pos_tmp < bins[k]['max'] ]
			for l in bin_value:
				bins_tmp[l][i]['positions'].append(pos_tmp)
				bins_tmp[l][i]['kxy'].append(kxy_tmp)
	#		print(" ".join([ str(totalData[rep*nRep]['params'][k]) for k in totalData[rep*nRep]['params'] ])) # print test of homogeneity in parameters over replicates
	res = {}
	for i in bins: # compute statistics #1 over bins and over #2 replicates
		meanPi_tmp = [] # values of pi at bin 'i' for different replicates
		stdPi_tmp = []
		pearsonR_tmp = []
		pearsonPval_tmp = []
		tajimaD_tmp = []
		thetaW_tmp = []
		achazY_tmp = []
		nHaplo_tmp = []
		H1_tmp = []
		H2_tmp = []
		H12_tmp = []
		H2overH1_tmp = []
		D_tmp = []
		r2_tmp = []
		pos_tmp = []
		size_tmp = regionSize * (bins[i]['max'] - bins[i]['min'])
		for j in range(nRep):
			nS = len(bins_tmp[i][j]['positions']) # number of SNPs
	#		meanPi_tmp.append(mean(bins_tmp[i][j]['pi']) / size_tmp)
			meanPi_tmp.append(sum(bins_tmp[i][j]['kxy']) /(1.0 * nPairwiseComp)) # sum(Kxy over comparisons at a SNP) / number_pairwise_comparison 
			stdPi_tmp.append(stdCustom(bins_tmp[i][j]['kxy'], meanPi_tmp[j], regionSize))
			if len(bins_tmp[i][j]['positions'])<2:
#				pearson = [nan, nan]
				pearson = [1, 1]
			else:
				pearson = pearsonr(bins_tmp[i][j]['positions'], bins_tmp[i][j]['kxy'])
			pearsonR_tmp.append(pearson[0])
			pearsonPval_tmp.append(pearson[1])
			tajimaD_thetaW = tajimaD(nIndiv, meanPi_tmp[j], nS) # [tajD, thetaW]
			tajimaD_tmp.append(tajimaD_thetaW[0])
			thetaW_tmp.append(tajimaD_thetaW[1])
			achazY_tmp.append(achazY(nIndiv, bins_tmp[i][j]['kxy']))
			pos_tmp.append(mean(bins_tmp[i][j]['positions']))
			nHaplo_tmp.append(resLD[j][i]['nHaplo'])
			H1_tmp.append(resLD[j][i]['H1'])
			H2_tmp.append(resLD[j][i]['H2'])
			H12_tmp.append(resLD[j][i]['H12'])
			H2overH1_tmp.append(resLD[j][i]['H2overH1'])
			D_tmp.append(resLD[j][i]['D'])
			r2_tmp.append(resLD[j][i]['r2'])
			#def stdCustom(liste, moyenne, longueurRegion):
		res[i] = {}
		if meanPi_tmp.count(nan)==len(meanPi_tmp):
			res[i]['pi_avg'] = nan
		else:
			res[i]['pi_avg'] = round(nanmean(meanPi_tmp)/size_tmp, 5)
		
		if stdPi_tmp.count(nan)==len(stdPi_tmp):
			res[i]['pi_std'] = nan
		else:
			res[i]['pi_std'] = round(nanmean(stdPi_tmp)/size_tmp, 5)
		
		if pearsonR_tmp.count(nan)==len(pearsonR_tmp) or pearsonR_tmp.count('nan')==len(pearsonR_tmp):
#			res[i]['pearson_r'] = nan
			res[i]['pearson_r'] = 1
		else:
			res[i]['pearson_r'] = round(nanmean(pearsonR_tmp), 5)
		
		if pearsonPval_tmp.count(nan)==len(pearsonPval_tmp):
#			res[i]['pearson_pval'] = nan
			res[i]['pearson_pval'] = 0
		else:
			res[i]['pearson_pval'] = round(nanmean(pearsonPval_tmp), 5)

		if tajimaD_tmp.count(nan)==len(tajimaD_tmp):
#			res[i]['tajD'] = nan
			res[i]['tajD'] = 0
		else:
			res[i]['tajD'] = round(nanmean(tajimaD_tmp), 5)

		if thetaW_tmp.count(nan)==len(thetaW_tmp):
			res[i]['thetaW'] = nan
		else:	
			res[i]['thetaW'] = round(nanmean(thetaW_tmp)/size_tmp, 5)
		
		if achazY_tmp.count(nan)==len(achazY_tmp):
			res[i]['achazY'] = nan
		else:
			res[i]['achazY'] = round(nanmean(achazY_tmp), 5)
		
		if nHaplo_tmp.count(nan)==len(nHaplo_tmp):
			res[i]['nHaplo'] = nan
		else:
			res[i]['nHaplo'] = round(nanmean(nHaplo_tmp), 5)
		
		if H1_tmp.count(nan)==len(H1_tmp):
#			res[i]['H1'] = nan
#			if no SNP : assuming H1 = 1
			res[i]['H1'] = 1
		else:
			res[i]['H1'] = round(nanmean(H1_tmp), 5)
		
		if H2_tmp.count(nan)==len(H2_tmp):
#			res[i]['H2'] = nan
#			if no SNP : assuming H2 = 0
			res[i]['H2'] = 0
		else:
			res[i]['H2'] = round(nanmean(H2_tmp), 5)
		
		if H12_tmp.count(nan)==len(H12_tmp):
#			res[i]['H12'] = nan
#			if no SNP : assuming H12 = 1
			res[i]['H12'] = 1
		else:
			res[i]['H12'] = round(nanmean(H12_tmp), 5)
		
		if H2overH1_tmp.count(nan)==len(H2overH1_tmp):
#			res[i]['H2overH1'] = nan
#			if no SNP : assuming H2/H1 = 0
			res[i]['H2overH1'] = 0
		else:
			res[i]['H2overH1'] = round(nanmean(H2overH1_tmp), 5)
		
		if D_tmp.count(nan)==len(D_tmp):
#			res[i]['D'] = nan
#			if no SNP : assuming LD = 1
			res[i]['D'] = 1
		else:
			res[i]['D'] = round(nanmean(D_tmp), 5)
		
		if r2_tmp.count(nan)==len(r2_tmp):
#			res[i]['r2'] = nan
			res[i]['r2'] = 1
		else:
			res[i]['r2'] = round(nanmean(r2_tmp), 5)
		
#		if pos_tmp.count(nan)==len(pos_tmp):
#			res[i]['position'] = nan
#		else:
#			res[i]['position'] = round(nanmean(pos_tmp), 5)
		
		res[i]['position'] = nanmean([bins[i]['min'],bins[i]['max']])
	return(res)


# compute the average pi over replicates
# define the bin boundaries
bins = window(width, step)

# parse the ms output file
#inputFileName = "output_test.msms"
#totalData = parse_msms(inputFileName, nIndiv, nCombParam, nRep)
totalData = parse_msms_v2(inputFileName, nIndiv, nCombParam, nRep)

min_width = 0.1
#pouet = haploStats(haplotypes=totalData[0]['haplotypes'], bins=bins, positions=totalData[0]['positions'], width=width, min_width=min_width)

# compute the diversity for all SNPs, for all alignments
# commented block because: loop is now imbedded in the parse_msms function to save memory
# now: the alignments are treated during reading, not only after
#for i in totalData:
#	if totalData[i]['segsites'] != 0:
#		#tmp = comp_pi(totalData[i]['haplotypes'], totalData[i]['segsites']) # first version: return the average pi
#		#totalData[i]['pi'] = tmp['pi_allSNP'] 
#		#totalData[i]['pi_singletons'] = tmp['pi_singletons'] 
#		tmp = comp_differences(totalData[i]['haplotypes'], totalData[i]['segsites']) # second version: return the sum of differences
#		totalData[i]['kxy'] = tmp['kxy_allSNP'] 
#		totalData[i]['kxy_singletons'] = tmp['kxy_singletons'] 
#		del totalData[i]['haplotypes']	


# prepare the output files:
stats_out = ""
for i in bins:
	stats_out += 'pi_avg_bin{0}\t'.format(i)
	stats_out += 'thetaW_avg_bin{0}\t'.format(i)
	stats_out += 'pi_std_bin{0}\t'.format(i)
	stats_out += 'tajD_bin{0}\t'.format(i)
	stats_out += 'achazY_bin{0}\t'.format(i)
	stats_out += 'pearson_r_bin{0}\t'.format(i)
	stats_out += 'pearson_pval_bin{0}\t'.format(i)
	stats_out += 'nHaplo_bin{0}\t'.format(i)
	stats_out += 'H1_bin{0}\t'.format(i)
	stats_out += 'H2_bin{0}\t'.format(i)
	stats_out += 'H12_bin{0}\t'.format(i)
	stats_out += 'H2overH1_bin{0}\t'.format(i)
	stats_out += 'D_bin{0}\t'.format(i)
	stats_out += 'r2_bin{0}\t'.format(i)
stats_out = stats_out.strip() + "\n"


pos_out = ""
for i in bins:
	pos_out += "SNP_{0}\t".format(i)
pos_out = pos_out.strip() + "\n"


# treat the replicated datasets
for i in range(nCombParam): # loop over combination of parameters
#	print("simulation {0} over a total of {1}: {2}".format(i, nCombParam, time.strftime("%H:%M:%S")))
	# TODO: calling this function tajimaD(nInd, pi, nS, size)
	a = calc_window(i, nRep=nRep, bins=bins, regionSize=regionSize, nIndiv=nIndiv, totalData=totalData, min_width=min_width) # get summary statistics per bin for the replicate i
	for j in bins:
		stats_out += "{0}\t".format(a[j]['pi_avg'])
		stats_out += "{0}\t".format(a[j]['thetaW'])
		stats_out += "{0}\t".format(a[j]['pi_std'])
		stats_out += "{0}\t".format(a[j]['tajD'])
		stats_out += "{0}\t".format(a[j]['achazY'])
		stats_out += "{0}\t".format(a[j]['pearson_r'])
		stats_out += "{0}\t".format(a[j]['pearson_pval'])
		stats_out += "{0}\t".format(a[j]['nHaplo'])
		stats_out += "{0}\t".format(a[j]['H1'])
		stats_out += "{0}\t".format(a[j]['H2'])
		stats_out += "{0}\t".format(a[j]['H12'])
		stats_out += "{0}\t".format(a[j]['H2overH1'])
		stats_out += "{0}\t".format(a[j]['D'])
		stats_out += "{0}\t".format(a[j]['r2'])
		pos_out += "{0}\t".format(a[j]['position'])
	stats_out = stats_out.strip() + "\n"
	pos_out = pos_out.strip() + "\n"


outfile = open("{0}_positions.txt".format(rootOutputFileName), "w")
outfile.write(pos_out)
outfile.close()

outfile = open("{0}_sumStats.txt".format(rootOutputFileName), "w")
outfile.write(stats_out)
outfile.close()

# sorted rows
for tmp_i in totalData:
	sorted_rows = sort_rows(totalData[tmp_i]['haplotypes'])

	outfile = open('{0}_{1}_sorted_rows.txt'.format(rootOutputFileName, tmp_i), 'w')
	res = '{positions}\n{haplotypes}'.format(positions=' '.join([str(i) for i in totalData[tmp_i]['positions']]), haplotypes=sorted_rows)
	outfile.write(res)
	outfile.close()

#######
# TMP #
#######

