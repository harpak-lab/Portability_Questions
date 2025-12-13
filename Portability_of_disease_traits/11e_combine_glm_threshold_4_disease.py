import argparse
import csv
import math

def process(paths, keep_path, output_path):
    # Read the list of SNP IDs to keep
    with open(keep_path, 'r') as f:
        keep_ids = set(f.read().strip().split('\n'))
    
    # Define the output columns
    fieldnames = [
        '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'A1', 'TEST', 'OBS_CT', 
        'OR', 'LOG(OR)_SE', 'Z_STAT', 'P', 'ERRCODE', 'BETA'
    ]
    
    # Open the output file and set up the CSV writer
    with open(output_path, 'w') as output_file:
        output_writer = csv.DictWriter(output_file, fieldnames=fieldnames, delimiter='\t')
        output_writer.writeheader()
        
        # Process each input file
        for path in paths:
            with open(path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for line in reader:
                    # Check if the current SNP ID is in the list to keep
                    if line['ID'] in keep_ids:
                        # Calculate BETA from OR if it's available
                        if 'OR' in line and line['OR'] != '.':
                            try:
                                line['BETA'] = str(math.log(float(line['OR'])))
                            except ValueError:
                                line['BETA'] = '.'
                        else:
                            line['BETA'] = '.'
                        
                        # Filter out extra columns that are not in fieldnames
                        filtered_line = {key: line.get(key, '.') for key in fieldnames}
                        
                        # Write the filtered line to the output file
                        output_writer.writerow(filtered_line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=(
        'Combine Plink .glm.logistic files, keeping only certain loci.'))
    parser.add_argument('paths', help='paths to .glm.logistic files', nargs='+')
    parser.add_argument('-k', '--keep', help='file of SNP IDs to keep')
    parser.add_argument('-o', '--output', help='path to the combined .glm.logistic file to be created')
    args = parser.parse_args()

    process(args.paths, args.keep, args.output)

