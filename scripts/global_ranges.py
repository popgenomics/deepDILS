# Script that takes as input file the name of a file containing the paths to different files range_stats.txt
import sys

# Get the filename of a list of range_stats.txt files from command-line arguments
list_files_names = sys.argv[1]

# Open the file of filenames and read in each filename
infile = open(list_files_names, 'r')
list_files = []
for line in infile:
	# Strip any whitespace characters (e.g., \n) from the filename and add it to the list
	list_files.append(line.strip())
infile.close()

# Loop over each range_stats.txt file in the list
stats = {} # A dictionary to store the stats found in the files
list_stats = [] # A list to store the names of the stats found in the files
for file_tmp in list_files:
	# Open the current file and read in each line
	infile = open(file_tmp, 'r')
	for line in infile:
		line = line.strip() # Strip any whitespace characters from the line
		# Split the line into the stat name and its value
		stat = line.split(': ')[0]
		value = line.split(': ')[1]
		# If the current stat hasn't been seen before, add it to the list of stats and initialize an empty list for it in the dictionary
		if stat not in stats:
			list_stats.append(stat)
			stats[stat] = []
		# Convert the value to a float and add it to the list for the current stat in the dictionary
		stats[stat].append(float(value))
	infile.close()

# Loop over each stat name and find the minimum or maximum value, depending on whether the stat contains the word "min" or not
for i in list_stats:
	if 'min' in i:
		value_tmp = min(stats[i])
	else:
		value_tmp = max(stats[i])
	# Print out the name of the stat and its minimum/maximum value
	print('{stat}: {value}'.format(stat=i, value=value_tmp))

