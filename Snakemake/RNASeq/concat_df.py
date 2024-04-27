#!/usr/bin/env python

# You can refer to the help manual by `python concat_df.py -h`

# argparse is a library that allows you to make user-friendly command line interfaces
import argparse

# here we are initializing the argparse object that we will modify
parser = argparse.ArgumentParser()

# we are asking argparse to require a -i or --input flag on the command line when this
# script is invoked. It will store it in the "filenames" attribute of the object
# we will be passing it via snakemake, a list of all the outputs of verse so we can
# concatenate them into a single matrix using pandas 

parser.add_argument("-i", "--input", help='A list of the VERSE output filenames provided by snakemake', dest="input", required=True, nargs='+')
parser.add_argument("-o", "--output", help='The output file name and path provided by snakemake', dest="output", required=True)

# this method will run the parser and input the data into the namespace object
args = parser.parse_args()


# if you try running this on the command line and supply it a value for -i or --input
# it will show up here, stored in this object

# try just running this script and supply it a random string for the -i argument
# example: `python concat_df.py -i <list of all files>`


import pandas
import os

concat = pandas.concat([pandas.read_csv(df, sep='\t', header=0, names = ['gene', '{}'.format(os.path.basename(df.split('.')[0]))], index_col='gene') for df in args.input], axis=1)

concat.to_csv(args.output)
