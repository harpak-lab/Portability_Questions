import argparse
import csv

def process(paths, keep_path, output_path):
    # Read the list of SNP IDs to keep
    with open(keep_path, 'r') as f:
        keep_ids = set(f.read().strip().split('\n'))
    
    # Define the columns we want to keep for the output
    output_columns = [
        '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'A1', 'TEST', 'OBS_CT', 
        'BETA', 'SE', 'T_STAT', 'P', 'ERRCODE'
    ]
    
    # Open the output file for writing
    with open(output_path, 'w') as output_file:
        output_writer = csv.DictWriter(output_file, fieldnames=output_columns, delimiter='\t')
        output_writer.writeheader()
        
        # Process each input file
        for path in paths:
            with open(path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                
                # Check if the file has the expected columns
                if not set(output_columns).issubset(reader.fieldnames):
                    print(f"Skipping file {path} due to missing columns")
                    continue
                
                # Filter rows based on the 'ID' column
                for line in reader:
                    if line['ID'] in keep_ids:
                        # Write only the selected columns to the output file
                        output_writer.writerow({col: line[col] for col in output_columns})

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter SNPs from .glm.linear files.')
    parser.add_argument('paths', nargs='+', help='Paths to .glm.linear files')
    parser.add_argument('-k', '--keep', required=True, help='File of SNP IDs to keep')
    parser.add_argument('-o', '--output', required=True, help='Path to the output file')
    args = parser.parse_args()
    
    process(args.paths, args.keep, args.output)
