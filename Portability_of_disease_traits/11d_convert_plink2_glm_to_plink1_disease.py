import argparse
import pandas as pd
import numpy as np

def convert_dataframe(dataframe):
    # Determine if logistic regression output
    is_logistic = 'OR' in dataframe.columns

    # Rename columns to match PLINK 1 format
    dataframe = dataframe.rename(columns={
        '#CHROM': '#CHR',
        'POS': 'BP',
        'ID': 'SNP',
        'OBS_CT': 'NMISS',
        'T_STAT': 'STAT'
    })

    # Add conditional columns based on the type of test
    if is_logistic:
        if 'OR' in dataframe.columns:
            dataframe['BETA'] = np.log(dataframe['OR'])
        else:
            dataframe['BETA'] = np.nan  # Handle missing OR column
        # Check if 'LOG(OR)_SE' exists before assigning to 'SE'
        if 'LOG(OR)_SE' in dataframe.columns:
            dataframe['SE'] = dataframe['LOG(OR)_SE']
        else:
            dataframe['SE'] = np.nan
    else:
        # For linear regression or other models
        if 'BETA' not in dataframe.columns:
            dataframe['BETA'] = np.nan
        if 'SE' not in dataframe.columns:
            dataframe['SE'] = np.nan

    # Filter relevant columns for PLINK 1 format
    return dataframe.filter(items=['#CHR', 'SNP', 'BP', 'A1', 'TEST', 'NMISS', 'BETA', 'STAT', 'P', 'SE'])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=(
        'Convert a Plink 2 association output (.glm.linear or .glm.logistic.hybrid) '
        'format to Plink 1 association output format (.assoc)'))
    parser.add_argument('path', help=('path to a Plink 2 GWAS output file '
                                      '(.glm.linear, .glm.logistic, or .glm.logistic.hybrid) '
                                      'to be converted'), nargs='+')
    parser.add_argument('-o', '--output',
                        help=('path to the Plink 1 output file to be created'))
    args = parser.parse_args()

    # If multiple paths are provided, merge all into a single dataframe
    if len(args.path) > 1:
        complete_df = pd.DataFrame()
        for path in args.path:
            df = pd.read_csv(path, sep='\t').pipe(convert_dataframe)
            complete_df = pd.concat([complete_df, df])
    else:
        complete_df = (
            pd.read_csv(args.path[0], sep='\t')
            .pipe(convert_dataframe)
        )
    
    # Drop duplicates and save the converted file
    (
        complete_df
        .drop_duplicates(subset=['SNP'])
        .to_csv(args.output, index=False, sep='\t', float_format='%.9g')
    )

