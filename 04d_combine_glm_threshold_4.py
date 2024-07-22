import argparse
import csv


def process(paths, keep_path, output_path):
    with open(keep_path, 'r') as f:
        keep_ids = set(f.read().strip().split('\n'))
    output_file = open(output_path, 'w')
    output_writer = csv.DictWriter(output_file, fieldnames=[
        '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'A1', 'TEST', 'OBS_CT', 'BETA',
        'SE', 'T_STAT', 'P', 'ERRCODE'], delimiter='\t')
    output_writer.writeheader()
    for path in paths:
        f = open(path, 'r')
        reader = csv.DictReader(f, delimiter='\t')
        for line in reader:
            if line['ID'] in keep_ids:
                output_writer.writerow(line)
        f.close()

    output_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=(
        'Combine Plink .glm.linear files, keeping only certain loci.'))
    parser.add_argument('paths', help=('paths to .glm.linear files'), nargs='+')
    parser.add_argument('-k', '--keep', help='file of SNP IDs to keep')
    parser.add_argument('-o', '--output',
                        help=('path to the combined .glm.linear  file to be created'))
    args = parser.parse_args()

    process(args.paths, args.keep, args.output)
