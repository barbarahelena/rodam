import os
import csv
import re

input_file = "input_file_paths.txt"
output_file = "Samplesheet.csv"

samples = {}

# Regular expression to match sample IDs from file names
pattern = r'^(.+)\.(R[12])\.fastq\.gz$'

with open(input_file, 'r') as infile:
    for line in infile:
        line = line.strip()
        match = re.match(pattern, os.path.basename(line))
        if match:
            sample_id, read_type = match.group(1), match.group(2)
            if sample_id not in samples:
                samples[sample_id] = {"fastq_1": "", "fastq_2": ""}
            
            if read_type == "R1":
                samples[sample_id]["fastq_1"] = line
            elif read_type == "R2":
                samples[sample_id]["fastq_2"] = line

# Debug: Print the extracted sample information
for sample_id, data in samples.items():
    print("Sample ID:", sample_id)
    print("Fastq 1:", data["fastq_1"])
    print("Fastq 2:", data["fastq_2"])

with open(output_file, 'w', newline='') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(["sample", "fastq_1", "fastq_2"])
    
    for sample_id, data in samples.items():
        writer.writerow([sample_id, data["fastq_1"], data["fastq_2"]])

print("CSV file created: Samplesheet.csv")
