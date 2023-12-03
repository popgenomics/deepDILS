#!/usr/bin/python
import sys
import argparse

# Create the argument parser
parser = argparse.ArgumentParser(description='Script that reads the output of msms (INFILE) when the -oTrace argument recording frequencies over time has been specified. get_msms_trajectory will extract the frequencies and write them to the OUTFILE file.')

# add the 'infile' argument
parser.add_argument('--infile', type=str, help='Infile name')
parser.add_argument('--outfile', type=str, help='Outputfile name')

# Analyse the arguments from the command line
args = parser.parse_args()

# Access to the arguments
infile = args.infile
outfile_name = args.outfile

def extract_block(file_path):
	with open(file_path, 'r') as file:
		lines = file.readlines()

	# Initializing variables for the start and end indices of the block
	start_block = None
	end_block = None

	# Finding the start index of the block
	for i, line in enumerate(lines):
		if line.startswith('Frequency Trace:'):
			start_block = i + 1  # The next index is the start of the block

	# Finding the end index of the block
	if start_block is not None:
		for i, line in enumerate(lines[start_block:], start=start_block):
			if line.startswith('segsites:'):
				end_block = i - 1  # The previous index is the end of the block
				break  # Stop searching once the end of the block is found

	# Extracting the block of lines
	if start_block is not None and end_block is not None:
		block_lines = lines[start_block:end_block]
		if block_lines[0] == '\n':
			block_lines = block_lines[1::]
		return block_lines
	else:
		return None


def write_block(outfile_name, block):
	with open(outfile_name, 'w') as outfile:
		nPop = int((len(block[0].strip().split('\t'))-1)/2)
		header = 'time\t' + '\t'.join([ 'freq_ancestral_pop{i}\tfreq_derived_pop{i}'.format(i=i) for i in range(nPop) ]) + '\n'
		
	#		outfile.write('time\tfreq_ancestral_pop1\tfreq_derived_pop1\tfreq_ancestral_pop2\tfreq_derived_pop2\n')
		outfile.write(header)
		for line in block:
			outfile.write(line)

# Calling the function by passing the file path
block = extract_block(file_path=infile)
write_block(outfile_name=outfile_name, block=block)

