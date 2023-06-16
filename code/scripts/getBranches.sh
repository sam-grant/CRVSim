# Pull the branch names and data types out of tr->Print() text file

#!/bin/bash

input_file="../txt/trkana-branches.txt" 
output_file="../txt/trkana-branches-out.txt" 

# Clear old output file
if [ -f $output_file ]; then
	rm -f $output_file && touch $output_file
fi

# Check if the file exists
if [ -f "$input_file" ]; then
  # Read the file line by line
  while read -r line; do
    # Get the line which includes the var name and dtype 
	result=$(grep "Br\s" <<< $line)
	 # Check if the result is empty
    if [ -n "$result" ]; then
		# Select everything after the first colon
		result=${result#*:}
		# Remove the trailing spaces and asterisk (expression from ChatGPT)
		result=$(echo $result | sed 's/[\*\ ]//g')
		# Rearrange string so the dtype comes type, remove the asterisk 
		result=$(echo $result | awk -F ':' '{print $2, $1}')
		# Print the result
		echo $result >> $output_file
	fi
  done < "$input_file"
else
  echo "File not found: $input_file"
fi

echo "Output file ${output_file} generated."