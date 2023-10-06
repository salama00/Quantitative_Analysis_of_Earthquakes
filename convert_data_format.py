'''
Salama Algaz
Date: May 8, 2023
Last Updated: Septemeber 30,2023

Reads in data output (ASCII text file) from Wilber3 Iris.edu. Creates a new data
file that is easy to read for Matlab. 

Format of the ASCII text file:
Column 1: Time		(seconds) (reference to initial earthquake (Mb 7.8) at 2023-02-06T01:17:34)
Column 2: Latitude	(degrees)
Column 3: Longitude	(degrees)
Column 4: Depth		(km)

There are 458 events, so there are 458 lines.
'''

################## Packages ######################

import numpy as np
from datetime import datetime

#################### Main ########################

# Number of lines in data file.
num_lines = 458		# Number of events

# Initial event (date of earthquake that caused subsequent earthquakes/aftershocks)
date0 = datetime(2023, 2, 6, 1, 17, 34)

# Open data file to read
file = open('downloaded_data.txt', 'r');

# Create a matrix that contains time in col 1 and positions: latitude, longitude, and depth 
# in cols 2, 3, and 4 respectively for each event. Should be a 458 x 4 matrix
events = np.zeros((num_lines, 4), dtype=np.ndarray)

# Go through the data file and read in the numbers
for i in range(num_lines):
	
	# Split the line 
	line = file.readline().split('|')
	
	# Get date of event
	date = datetime.strptime(line[1], '%Y-%m-%dT%H:%M:%S')
	# Get time elapsed from initial earthquake to current event
	time_elapsed = date - date0
	
	# Time elapsed in seconds
	dt = time_elapsed.total_seconds()
	events[i][0] = dt					# Add time elapsed to matrix
	
	# Get positions
	events[i][1] = float(line[2])		# Read in latitude	(degrees)
	events[i][2] = float(line[3])		# Read in longitude	(degrees)
	events[i][3] = float(line[4])		# Read in depth		(km)

# Close file
file.close()

# Write new file for matlab:
mat_file = open('matlab_data.txt', 'w');

for e in events:
	mat_file.write(str(e[0]) + '  ' + str(e[1]) + '  ' + str(e[2]) + '  ' + str(e[3]) + '\n' )

# Close file
mat_file.close()