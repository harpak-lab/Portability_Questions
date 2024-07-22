import argparse

import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=(
        'Filter SNPs from a GWAS output file by a p-value threshold.'))
    parser.add_argument('path', help='path to GWAS output file')
    parser.add_argument('-t', '--threshold', type=float,
                        help='Maximum p-value for included SNPs')
    parser.add_argument('-o', '--output',
                        help=('path to the SNPs file to be created. This file'
                              'has one variant identifier per line, only.'))

    args = parser.parse_args()

    (
        pd.read_csv(args.path, sep=r'\s+')
        .drop_duplicates(subset=['SNP'], keep=False)
        .loc[lambda df: df['P'] <= args.threshold, 'SNP']
        .to_csv(args.output, sep=' ', index=False, float_format='%.9g',
                header=False)
    )
