import argparse
import matplotlib.pyplot as plt

def plot_haplotype_matrix(file_path, output_file, resolution):
	# read the file
	with open(file_path, 'r') as file:
		lines = file.readlines()[1:]  # ignore the first line

	# creating the matrix from readen lines, and inversion of the order of rows
	matrix = [[int(char) for char in line.strip()] for line in lines][::-1]

	# creating the plot
	plt.figure(figsize=(resolution[0]/100, resolution[1]/100))
	plt.imshow(matrix, cmap='binary', aspect='auto', origin='lower')
	plt.axis('off')  # removing axis
	plt.xticks([])   # remove the ticks on the x-axis
	plt.yticks([])   # remove the ticks on the y-axis
	plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

	# save the picture in a PNG file with the specified resolution
	plt.savefig(output_file, bbox_inches='tight', pad_inches=0, dpi=resolution[0]/(resolution[0]/100))
	plt.close()  # close the plot after exportation

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Generate a haplotype matrix image from a msms file produced by the msms simulator')
	parser.add_argument('--input', type=str, default='simulation.msms', help='path to the input file (default: %(default)s)')
	parser.add_argument('--output', type=str, default='simulation.png', help='path to save the output PNG image (default: %(default)s)')
	parser.add_argument('-S', '--segregating_sites', type=int, default=2000, help='Number of segregating sites (default: 2000)')
	parser.add_argument('-n', '--sampled_gametes', type=int, default=40, help='Number of sampled gametes (default: 40)')

	args = parser.parse_args()
	resolution = (args.segregating_sites, args.sampled_gametes)
	plot_haplotype_matrix(args.input, args.output, resolution)
