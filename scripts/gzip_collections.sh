#!/bin/bash

# Define source and target directories
RAW_DIR="rawdata"
PROC_DIR="procdata"

# Ensure the target directory exists
mkdir -p "$PROC_DIR"

# Iterate over directories in rawdata
for dir in "$RAW_DIR"/*; do
	if [ -d "$dir" ]; then
		# Extract the directory name
		dir_name=$(basename "$dir")
		
		# Create the tar.gz archive in procdata
		tar cf - "$dir" | pigz -c > "$PROC_DIR/${dir_name}.tar.gz"
		
		echo "Archived $dir to $PROC_DIR/${dir_name}.tar.gz"
	fi
done