import argparse
import csv
import glob
import os


def process(paths, keep_path, output_path):
    with open(keep_path, 'r') as f:
        keep_ids = set(line.strip() for line in f if line.strip())

    with open(output_path, 'w', newline='') as output_file:
        output_header = [
            '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'A1', 'TEST', 'OBS_CT',
            'BETA', 'SE', 'T_STAT', 'P'
        ]
        output_writer = csv.DictWriter(output_file, fieldnames=output_header, delimiter='\t')
        output_writer.writeheader()

        for path_pattern in paths:
            for path in glob.glob(path_pattern):
                with open(path, 'r') as f:
                    reader = csv.DictReader(f, delimiter=' ')
                    reader.fieldnames = [name.strip() for name in reader.fieldnames]

                    if 'ID' not in reader.fieldnames:
                        raise ValueError(f"'ID' column not found in file {path}.")

                    for line in reader:
                        if line['ID'].strip() in keep_ids:
                            output_writer.writerow({
                                '#CHROM': line['CHROM'],
                                'POS': line['GENPOS'],
                                'ID': line['ID'],
                                'REF': line['ALLELE0'],
                                'ALT': line['ALLELE1'],
                                'A1': line['ALLELE1'],
                                'TEST': line['TEST'],
                                'OBS_CT': line['N'],
                                'BETA': line['BETA'],
                                'SE': line['SE'],
                                'T_STAT': line['CHISQ'],
                                'P': f"{10 ** (-float(line['LOG10P'])):.9g}"
                            })


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=(
        'Combine Regenie output files, keeping only certain loci.'))
    parser.add_argument('paths', help=('Paths to Regenie output files (wildcards allowed)'), nargs='+')
    parser.add_argument('-k', '--keep', required=True,
                        help='File of SNP IDs to keep')
    parser.add_argument('-o', '--output', required=True,
                        help=('Path to the combined output file to be created'))
    args = parser.parse_args()

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    process(args.paths, args.keep, args.output)

