#!/usr/bin/env python

import sys
import argparse
import os.path
import pandas as pd
import csv
import gzip
import statistics

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s",
        "--saltgenes",
        required=True,
        metavar="FILE",
        help="TSV file containing contig depths" 
    )
    parser.add_argument(
        "-t",
        "--tax",
        required=True,
        metavar="FILE",
        help="TSV file containing taxonomy per bin per sample",
    )
    parser.add_argument(
        "-o", 
        "--out", 
        required=False, 
        metavar="FILE", 
        help="Output file to save the merged table",
        default="saltgenes_tax_depths.tsv"  # Default value if --out is not provided
    )
    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)

    # Load the data into a DataFrame without headers and set custom column names
    df = pd.read_csv(args.depths, sep='\t', header=0)
    df = df.rename(columns={"bin": "Bin"})
    # print(df.columns)  # Check if the headers are read correctly

    # Strip the .fa.gz suffix from the 'Bin' column
    df["Bin"] = df["Bin"].str.replace(".fa.gz", "", regex=False)

    # Load taxonomy table
    taxonomy_df = pd.read_csv(args.tax, sep='\t')

    # Merge the tables on 'Sample' and 'Bin'
    merged_df = pd.merge(df_long, taxonomy_df, on=["Sample", "Bin"], how="inner")

    # Save the merged table to the specified output file
    merged_df.to_csv(args.out, sep='\t', index=False)

    print(f"Merged table saved to {args.out}")

if __name__ == "__main__":
    sys.exit(main())
