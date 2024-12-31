import argparse
import pandas as pd


def convert_dataframe(dataframe):
    # Rename columns for compatibility
    dataframe = dataframe.rename(columns={
        'CHROM': '#CHR',
        'GENPOS': 'BP',
        'ID': 'SNP',
        'N': 'NMISS',
        'CHISQ': 'STAT',
        'ALLELE1': 'A1',
    })

    # Compute P-values from LOG10P
    if 'LOG10P' in dataframe.columns:
        dataframe['P'] = 10 ** (-dataframe['LOG10P'])
    else:
        raise KeyError("The required column 'LOG10P' is missing from the input file.")

    # Return only the relevant columns
    return dataframe.filter(items=['#CHR', 'SNP', 'BP', 'A1', 'TEST', 'NMISS', 
                                   'BETA', 'STAT', 'P', 'SE'])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=(
        'Convert a GWAS output file to a Plink 1-compatible association file'))
    parser.add_argument('path', help=('Path to GWAS output file to be converted'), nargs='+')
    parser.add_argument('-o', '--output',
                        help=('Path to the Plink 1-compatible output file'))
    args = parser.parse_args()

    # Combine multiple input files if provided
    if len(args.path) > 1:
        complete_df = pd.DataFrame()
        for path in args.path:
            df = pd.read_csv(path, delim_whitespace=True).pipe(convert_dataframe)
            complete_df = pd.concat([complete_df, df], ignore_index=True)
    else:
        complete_df = (
            pd.read_csv(args.path[0], delim_whitespace=True)
            .pipe(convert_dataframe)
        )
    (
        complete_df
        .drop_duplicates(subset=['SNP'])
        .to_csv(args.output, index=False, sep='\t', float_format='%.9g')
    )

