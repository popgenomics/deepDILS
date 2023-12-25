import argparse
import os

def read_sumstats_file(file_name: str) -> tuple:
	"""Reads the contents of the sumStats file."""
	with open(file_name, 'r') as file:
		lines = file.readlines()

	headers = lines[0].strip().split('\t')
	data = [line.strip().split('\t') for line in lines[1:]]

	return headers, data

def convert_to_matrix(headers: list, data: list) -> dict:
	"""Converts the data into a matrix format."""
	stats = {}
	bin_count = len(headers) // 2  # Number of bins

	for i in range(bin_count):
		for header in headers:
			if header.endswith(f"_bin{i}"):
				key = header.replace(f"_bin{i}", "")
				if key not in stats:
					stats[key] = []
				idx = headers.index(header)
				stats[key].extend(data[j][idx] for j in range(len(data)))

	return stats

def write_to_file(file_name: str, content: str):
	"""Writes content to a file."""
	with open(file_name, 'w') as file:
		file.write(content)

def process_files(input_folder: str):
	"""Processes files in the input folder."""
	for file_name in os.listdir(input_folder):
		if file_name.endswith("_sumStats.txt"):
			file_path = os.path.join(input_folder, file_name)
			headers, data = read_sumstats_file(file_path)
			stats = convert_to_matrix(headers, data)
			output_file = file_name.replace(".txt", ".mat")
			output_path = os.path.join(input_folder, output_file)
			content = display_matrix_v2(stats)
			write_to_file(output_path, content)

def display_matrix(stats: dict) -> str:
	"""Generates the matrix representation of statistics."""
	result = '\t'
	count = 0
	for stat_i in stats:
		if count == 0:
			result += '\t'.join([str(i) for i in range(len(stats[stat_i]))]) + '\n'
		count += 1
		result += '{stat_name_i}\t{values_i}\n'.format(stat_name_i=stat_i, values_i='\t'.join(stats[stat_i]))
	return result

def display_matrix_v2(stats: dict) -> str:
	"""Generates the matrix representation of statistics."""
	"""Without row names and col names."""
	result = ''
	for stat_i in stats:
		result += '{values_i}\n'.format(values_i='\t'.join(stats[stat_i]))
	return result

def main():
	parser = argparse.ArgumentParser(description='Converts sumStats.txt files to matrix files')
	parser.add_argument('input_folder', type=str, help='Path to the folder containing sumStats.txt files')
	args = parser.parse_args()

	process_files(args.input_folder)

if __name__ == "__main__":
	main()

