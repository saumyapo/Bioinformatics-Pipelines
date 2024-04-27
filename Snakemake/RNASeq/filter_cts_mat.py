#!/usr/bin/env python

# You can refer to the help manual by `python filter_cts_mat.py -h`

# argparse is a library that allows you to make user-friendly command line interfaces
import argparse

# here we are initializing the argparse object that we will modify
parser = argparse.ArgumentParser()

# we are asking argparse to require a -i or --input flag on the command line when this
# script is invoked. It will store it in the "filenames" attribute of the object. Here
# we are only asking to provide this script one file: the GTF file we are parsing
parser.add_argument("-i", "--input", help='The input file specified will be the concatenated VERSE matrix provided by snakemake',dest="input", required=True)
parser.add_argument("-o", "--output", help='The output file name and path provided by snakemake',dest="output", required=True)

# this method will run the parser and input the data into the namespace object
args = parser.parse_args()

# if you try running this on the command line and supply it a value for -i or --input
# it will show up here, stored in this object

# try just running this script and supply it a random string for the -i argument
# example: `python filter_cts_mat.py -i <your_string> -o <another_string>`

# replace the code that comes after this with the code necessary to filter the dataframe

import pandas

df = pandas.read_csv(args.input)

sample_cols = [col for col in df.columns if col != 'gene']

filtered = df.loc[(df[sample_cols] != 0).sum(axis=1) == 8]

filtered.to_csv(args.output, index=False)

