import pandas as pd
import os
import argparse

def process_files(directory, output_csv):
	# Initialize a list to store data from each file
	data = []

	# Traverse through the files in the directory
	for file_name in os.listdir(directory):
		if file_name.endswith("_sumStats.txt"):
			file_path = os.path.join(directory, file_name)
			# Extract information from the file name
			file_name_parts = file_name.split('-')
			if len(file_name_parts) == 4:
				model = file_name_parts[0]
				iteration_selection = file_name_parts[3].split('_')
				dateRun = file_name_parts[1]
				timeStamp = file_name_parts[2]
				if len(iteration_selection) == 3:
					iteration = iteration_selection[0]
					selection = iteration_selection[1]
					# Read the file and store data
					with open(file_path, 'r') as file:
						# Read file data using pandas
						df = pd.read_csv(file, sep='\t')

						# Add columns: model, selection, iteration
						df['model'] = model
						df['selection'] = selection
						df['dateRun'] = dateRun 
						df['timeStamp'] = timeStamp 
						df['iteration'] = iteration

						# Rearrange columns order
						cols = ['model', 'selection', 'dateRun', 'timeStamp', 'iteration'] + [col for col in df if col not in ['model', 'selection', 'dateRun', 'timeStamp', 'iteration']]
#						cols = ['model', 'selection', 'iteration'] + [col for col in df if col not in ['model', 'selection', 'iteration']]
						df = df[cols]

						# Add the DataFrame to the list
						data.append(df)
				else:
					print(f"Issue with file: {file_name}. Skipping this file.")

	if data:  # Check if data is not empty
		# Concatenate all DataFrames into one
		global_table = pd.concat(data, ignore_index=True, sort=False)

		# Write the global table to a CSV file with tab delimiter
		global_table.to_csv(output_csv, sep='\t', index=False, na_rep='NaN')
		print(f"Combined data written to {output_csv}")
	else:
		print("No valid files found to process.")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Combine sumStats.txt files into a single CSV.')
	parser.add_argument('directory', help='Input directory containing sumStats.txt files')
	parser.add_argument('output_csv', help='Output CSV file name')

	args = parser.parse_args()
	process_files(args.directory, args.output_csv)


