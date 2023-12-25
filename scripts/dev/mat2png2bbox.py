import numpy as np
import matplotlib.pyplot as plt
import argparse
import glob
import os
import tarfile

def read_prior_file(base_file_name):
	"""
	Reads data from a '.prior' file corresponding to the given base file name.
	The file is expected to contain a header line and a data line.
	
	Parameters:
	- base_file_name (str): The base name of the file to read.
	
	Returns:
	- dict: A dictionary with keys from the header and corresponding values from the data line.
	"""
	prior_file = f'{base_file_name}.prior'
	with open(prior_file, 'r') as file:
		lines = file.readlines()
		header = lines[0].split()
		data = lines[1].split()
		prior_values = {header[i]: data[i] for i in range(len(header))}
		return prior_values

def scale_data(simulation_data, ranges_data):
	"""
	Normalizes each row of simulation data based on the corresponding min and max values in ranges_data.
	
	Parameters:
	- simulation_data (numpy.ndarray): Array containing the simulation data to be scaled.
	- ranges_data (numpy.ndarray): Array containing min and max values for scaling each row of simulation_data.
	
	Returns:
	- numpy.ndarray: The scaled simulation data.
	"""
	scaled_data = np.zeros_like(simulation_data)
	for i, (sim_row, range_row) in enumerate(zip(simulation_data, ranges_data)):
		min_val, max_val = range_row[0], range_row[1]
		scaled_data[i] = (sim_row - min_val) / (max_val - min_val)
	return scaled_data

def save_scaled_data(scaled_data, output_file):
	"""
	Saves scaled data to a specified output file in a text format.
	
	Parameters:
	- scaled_data (numpy.ndarray): The scaled data to be saved.
	- output_file (str): The path to the output file where the data will be saved.
	"""
	np.savetxt(output_file, scaled_data)

def save_bounding_box(object_class, x_center, y_center, x_width, y_width, output_file):
	"""
	Saves bounding box information to a specified output file.
	The information is saved in a tab-separated format.

	Parameters:
	- object_class (str): The class of the object.
	- x_center, y_center (float): Center coordinates of the bounding box.
	- x_width, y_width (float): Width and height of the bounding box.
	- output_file (str): The path to the output file.
	"""
	data_str = f"{object_class}\t{x_center}\t{y_center}\t{x_width}\t{y_width}\n"
	with open(output_file, 'w') as file:
		file.write(data_str)

def find_matching_files(base_file_name):
	"""
	Finds files matching specific patterns based on the base file name.
	Looks for files ending with '_neutral_sumStats.mat' and '_sweep_sumStats.mat'.

	Parameters:
	- base_file_name (str): The base file name to use for pattern matching.
	
	Returns:
	- tuple: A tuple containing two lists, one for each file pattern.
	"""
	neutral_file = glob.glob(f'{base_file_name}_neutral_sumStats.mat')
	sweep_file = glob.glob(f'{base_file_name}_sweep_sumStats.mat')
	return neutral_file, sweep_file

def generate_image(scaled_data, output_image):
	"""
	Generates and saves an image from the provided scaled data.
	The image is saved without axes and margins.

	Parameters:
	- scaled_data (numpy.ndarray): The data to be visualized in the image.
	- output_image (str): Path to the output image file.
	"""
	inverse_cmap = plt.cm.gray_r
	plt.figure(figsize=(scaled_data.shape[1] / 100, scaled_data.shape[0] / 100), dpi=100)
	plt.imshow(scaled_data, cmap=inverse_cmap, vmin=0, vmax=1)
	plt.axis('off')
	plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
	plt.savefig(output_image, bbox_inches='tight', pad_inches=0)

def read_sweep_positions(sweep_positions_file):
	"""
	Reads SNP (Single Nucleotide Polymorphisms) positions from a file.

	Parameters:
	- sweep_positions_file (str): The path to the file containing SNP positions.

	Returns:
	- numpy.ndarray: An array of SNP positions.
	"""
	with open(sweep_positions_file, 'r') as file:
		lines = file.readlines()
		snp_positions = np.array([float(pos) for pos in lines[1].split()])
		return snp_positions

def find_closest_snp(relative_position, snp_positions):
	"""
	Finds the index of the closest SNP (Single Nucleotide Polymorphism) to a given relative position.

	Parameters:
	- relative_position (float): The relative position to compare against SNP positions.
	- snp_positions (numpy.ndarray): An array of SNP positions.

	Returns:
	- int: The index of the closest SNP.
	"""
	return np.argmin(np.abs(snp_positions - relative_position))

def delimit_region(scaled_data, pi_neutral_mean, relative_selected_position, snp_positions):
	"""
	Delimits a region around the relative_selected_position where the values in scaled_data[0] are less than or equal to pi_neutral_mean.

	Parameters:
	scaled_data (numpy.ndarray): The scaled simulation data.
	pi_neutral_mean (float): The mean value of pi calculated from neutral simulations, used as a threshold.
	relative_selected_position (float): The relative position of interest in the genome.
	snp_positions (numpy.ndarray): Array of SNP positions in the genome.

	Returns:
	tuple: A tuple containing the center and width of the delimited region. Both are float values.
	"""

	# Find the index of the closest SNP to the relative selected position
	center_index = find_closest_snp(relative_selected_position, snp_positions)
	
	# Initialize left and right indices for the delimited region
	left_index, right_index = center_index, center_index

	# Expand the region to the left
	while left_index > 0 and scaled_data[0][left_index] <= pi_neutral_mean:
		left_index -= 1
	left_index += 1  # Adjust to include the last valid point
	
	# Expand the region to the right
	while right_index < len(scaled_data[0]) - 1 and scaled_data[0][right_index] <= pi_neutral_mean:
		right_index += 1
	right_index -= 1  # Adjust to include the last valid point

	# Calculate the center and width of the region
	region_center = (snp_positions[left_index] + snp_positions[right_index]) / 2
	region_width = snp_positions[right_index] - snp_positions[left_index]

	# Adjust the region to ensure it stays within [0, 1]
	if region_center - region_width/2 < 0:
		region_width = region_center*2
	elif region_center + region_width/2 > 1:
		region_width = (1 - region_center)*2

#	# Original function
#	x_center = round(region_center, 5)
#	x_width = round(region_width, 5)
	# Takes here as center_position and width the windows and not the SNP position
	if left_index>0:
		SNP_pos_left = float(snp_positions[left_index-1])
		left_index = (left_index-1)/len(snp_positions)
	else:
		SNP_pos_left = float(snp_positions[left_index])
		left_index = (left_index)/len(snp_positions)

	if right_index<len(snp_positions)-1:
		SNP_pos_right = float(snp_positions[right_index+1])
		right_index = (right_index+1)/len(snp_positions)
	else:
		SNP_pos_right = float(snp_positions[right_index])
		right_index = (right_index)/len(snp_positions)
	x_center = round((left_index+right_index)/2, 5)
	x_width = right_index-left_index
	y_center = 0.5
	y_width = 1
#	return 1, x_center, y_center, x_width, y_width
	return 1, x_center, y_center, x_width, y_width, SNP_pos_left, SNP_pos_right

def extract_snp_positions(base_file_name):
	"""
	Extracts SNP positions from a specified file inside a tar.gz archive.

	Parameters:
	- archive_path (str): Path to the tar.gz archive.
	- file_name (str): Name of the file inside the archive from which to extract SNP positions.

	Returns:
	- list: List of SNP positions.
	"""
	archive_path = base_file_name + '_archive_msms.tar.gz'
	file_name = base_file_name + '_sweep_0_sorted_rows.msms'
	with tarfile.open(archive_path, "r:gz") as tar:
		# Rechercher le fichier spécifique dans l'archive
		file = tar.extractfile(file_name)
		if file is not None:
			# Lire la première ligne et convertir en liste de floats
			line = file.readline().decode('utf-8')
			all_snp_positions = np.array([float(pos) for pos in line.split()])
			return all_snp_positions
		else:
			raise FileNotFoundError(f"{file_name} not found in {archive_path}")

def process_simulation_file(file_name, ranges_data, pi_neutral_mean=None, relative_selected_position=None, base_file_name=None):
	"""
	Processes a single simulation file - scales data, generates images, and handles specific operations for neutral and sweep files.
	"""
	with open(file_name, 'r') as simulation_file:
		simulation_data = np.array([[float(val) for val in line.split()] for line in simulation_file.readlines()])

	# Scale the simulation data
	scaled_data = scale_data(simulation_data, ranges_data)

	# Save the scaled data
	output_file = file_name.replace('.mat', '.rescaled')
	save_scaled_data(scaled_data, output_file)

	# Generate and save an image
	output_image = file_name.replace('.mat', '_globalPic.png')
	generate_image(scaled_data, output_image)

	if pi_neutral_mean is not None and relative_selected_position is not None and base_file_name is not None:
		# Additional operations for sweep file
		handle_sweep_file_operations(file_name, scaled_data, pi_neutral_mean, relative_selected_position, base_file_name)

def handle_sweep_file_operations(file_name, scaled_data, pi_neutral_mean, relative_selected_position, base_file_name):
	"""
	Handles operations specific to sweep files like bounding box generation.
	"""
	snp_positions_file = base_file_name + '_sweep_positions.txt'
	snp_positions = read_sweep_positions(snp_positions_file)

	object_class, x_center, y_center, x_width, y_width, snp_pos_left, snp_pos_right = delimit_region(scaled_data, pi_neutral_mean, relative_selected_position, snp_positions)
	output_file = file_name.replace('.mat', '_globalPic.txt')
	save_bounding_box(object_class, x_center, y_center, x_width, y_width, output_file)

	# Handle raw data operations
	handle_raw_data_operations(base_file_name, snp_pos_left, snp_pos_right)

def handle_raw_data_operations(base_file_name, snp_pos_left, snp_pos_right):
	"""
	Handles operations specific to raw data like defining bounding box coordinates.
	"""
	all_snp_positions = extract_snp_positions(base_file_name=base_file_name)
	nSNPs = len(all_snp_positions)
	closest_snp_index_left_rawData = find_closest_snp(snp_pos_left, all_snp_positions)
	closest_snp_index_right_rawData = find_closest_snp(snp_pos_right, all_snp_positions)
	center_raw = (closest_snp_index_left_rawData + closest_snp_index_right_rawData) / (2 * nSNPs)
	width_raw = (closest_snp_index_right_rawData - closest_snp_index_left_rawData) / nSNPs
	output_file = base_file_name + '_sweep_0_sorted_rows.txt'
	save_bounding_box(1, center_raw, 0.5, width_raw, 1, output_file)

# Main function that orchestrates the process
def main():
	parser = argparse.ArgumentParser(description="Process simulation data and generate image")
	parser.add_argument("ranges_file", help="Path to ranges file")
	parser.add_argument("base_file_name", help="Base file name without suffix")
	args = parser.parse_args()

	neutral_file, sweep_file = find_matching_files(args.base_file_name)
	if not neutral_file or not sweep_file:
		print("Matching files were not found.")
		return

	prior_values = read_prior_file(args.base_file_name)
	ranges_data = np.loadtxt(args.ranges_file)

	# Process neutral file
	process_simulation_file(neutral_file[0], ranges_data)
	
	# Calculate pi_neutral_mean for use with sweep files
	position_selected_allele = float(prior_values.get("position_selected_allele"))
	chromosome_length = float(prior_values.get("chromosome_length"))
	relative_selected_position = position_selected_allele / chromosome_length
	pi_neutral_mean = np.mean(scale_data(np.loadtxt(neutral_file[0]), ranges_data)[0])

	# Process sweep file
	process_simulation_file(sweep_file[0], ranges_data, pi_neutral_mean, relative_selected_position, args.base_file_name)

if __name__ == "__main__":
	main()

