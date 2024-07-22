import argparse

import pandas as pd


def convert_dataframe(dataframe):
    return (
        dataframe
        .rename(columns={
            '#CHROM': '#CHR',
            'POS': 'BP',
            'ID': 'SNP',
            'OBS_CT': 'NMISS',
            'T_STAT': 'STAT',
        })
        .filter(items=['#CHR', 'SNP', 'BP', 'A1', 'TEST', 'NMISS', 'BETA',
                       'STAT', 'P', 'SE'])
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=(
        'Convert a Plink 2 association output (.glm.linear) format to '
        'Plink 1 association output format (.assoc)'))
    parser.add_argument('path', help=('path to a Plink 2 GWAS output file '
                                      '(.glm.linear or .glm.logistic) '
                                      'to be converted'), nargs='+')
    parser.add_argument('-o', '--output',
                        help=('path to the Plink 1 output file to be created'))
    args = parser.parse_args()

    # If multiple paths were passed, add them all to a single table
    if len(args.path) > 1:
        complete_df = pd.DataFrame()
        for path in args.path:
            df = pd.read_csv(path, sep='\t').pipe(convert_dataframe)
            complete_df = pd.concat([complete_df, df])
    # If just one path passed, just convert and save
    else:
        complete_df = (
            pd.read_csv(args.path[0], sep='\t')
            .pipe(convert_dataframe)
        )
    (
        complete_df
        .drop_duplicates(subset=['SNP'])
        .to_csv(args.output, index=False, sep='\t', float_format='%.9g')
    )
