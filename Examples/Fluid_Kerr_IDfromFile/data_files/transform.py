import csv
import os
import glob

# Define the directory containing the CSV files
input_dir = './'

# Get a list of all CSV files in the directory
csv_files = glob.glob(os.path.join(input_dir, '*.csv'))

# Process each CSV file
for input_file in csv_files:
    output_file = input_file.replace('.csv', '_vs_r.txt')  # Create corresponding output file

    with open(input_file, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        
        # Open the new file for writing
        with open(output_file, 'w') as new_file:
            for row in csv_reader:
                for value in row:
                    new_file.write(value.strip() + '\n')  # Remove spaces and write each value on a new line

    print(f"Transformation complete: '{output_file}' has been created.")
