"""
Check if the R1 and R2 files are in the same order, based on fastqc_data.txt
If not, sys.exit
If yes, print "R1 and R2 files have the same number of sequences"
"""
import zipfile
import re
import sys
import os

# Path to the ZIP archive
r1_zip_path = sys.argv[1]  # "folder.zip"
r2_zip_path = sys.argv[2]  # "folder.zip"

# File inside the ZIP to read
r1_base_name = os.path.basename(r1_zip_path)
r2_base_name = os.path.basename(r2_zip_path)
r1_file_path = r1_base_name.replace(".zip", "") + "/fastqc_data.txt"
r2_file_path = r2_base_name.replace(".zip", "") + "/fastqc_data.txt"

r1_total_sequences = 0
r2_total_sequences = 0

# Open the ZIP file
with zipfile.ZipFile(r1_zip_path, 'r') as zip_file:
    # Open the specific file
    with zip_file.open(r1_file_path) as file:
        # Read and decode the contents
        content = file.read().decode('utf-8')
        match = re.search(r"Total Sequences\t(\d+)", content)
        if match:
            r1_total_sequences = int(match.group(1))
        else:
            sys.exit("fastqc_data.txt does not contain 'Total Sequences'")

with zipfile.ZipFile(r2_zip_path, 'r') as zip_file:
    # Open the specific file
    with zip_file.open(r2_file_path) as file:
        # Read and decode the contents
        content = file.read().decode('utf-8')
        match = re.search(r"Total Sequences\t(\d+)", content)
        if match:
            r2_total_sequences = int(match.group(1))
        else:
            sys.exit("fastqc_data.txt does not contain 'Total Sequences'")

if r1_total_sequences != r2_total_sequences:
    sys.exit("R1 and R2 files have different number of sequences")
else:
    print("R1 and R2 files have the same number of sequences")
