import os
import argparse

def process_files(directory):
	# Initialize with None to track minimum and maximum values for each line
	min_values = None
	max_values = None

	for root, _, files in os.walk(directory):
		for file in files:
			if file.endswith('_sumStats.mat'):
				file_path = os.path.join(root, file)

				# Read the file and retrieve values as a list of lists
				with open(file_path, 'r') as f:
					lines = f.readlines()
					values = [list(map(float, line.split())) for line in lines]

					# Initialize min and max values if not yet initialized
					if min_values is None and max_values is None:
						min_values = [min(i) for i in values]
						max_values = [max(i) for i in values]
					else:
						# Update min and max values for each line
						for idx, row in enumerate(values):
							min_values[idx] = min(min_values[idx], min(row))
							max_values[idx] = max(max_values[idx], max(row))

	return list(zip(min_values, max_values))

def main(directory):
	result_table = process_files(directory)

	# Display the output table
	for row in result_table:
		print(f"{row[0]}\t{row[1]}")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Find minimum and maximum values in *_sumStats.mat files")
	parser.add_argument("directory", help="Directory path to traverse for *_sumStats.mat files")
	args = parser.parse_args()

	# Call the main function with the specified directory argument
	main(args.directory)
