import re

# Specify the input and output file paths
input_file="test.txt" # "../txt/trkana-branches.txt"
output_file="out.txt" # ../txt/trkana-branches-out.txt"

# Regular expression pattern to extract branch names and data types
# pattern = r"\*\w+\s+:\w+\.\w+\s+:\s+(\w+)"
# pattern = r"\*\w+\s+:\w+\s+:\s+(\w+\.\w+)"
# pattern = r"\*\w+\s+:\s+(\w+)\s+:\s+(\w+\.\w+)"
pattern = r"^\*Br\s+\d+\s+:(\w+\.\w+)\s+:\s+(\w+)"

# Open the input and output files
with open(input_file, "r") as f_input, open(output_file, "w") as f_output:

    # Read the input file
    input_data = f_input.read()
    
    # Find all matches of the pattern in the input data
    matches = re.findall(pattern, input_data)


    
    # Write the matches to the output file
    for match in matches:
        print(match)
        f_output.write(match + "\n")
