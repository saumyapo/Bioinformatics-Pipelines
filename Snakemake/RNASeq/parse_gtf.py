#!/usr/bin/env python

# You can refer to the help manual by `python parse_gtf.py -h`

# argparse is a library that allows you to make user-friendly command line interfaces
import argparse

# here we are initializing the argparse object that we will modify
parser = argparse.ArgumentParser()

# we are asking argparse to require a -i or --input flag on the command line when this
# script is invoked. It will store it in the "filenames" attribute of the object. Here
# we are only asking to provide this script one file: the GTF file we are parsing
parser.add_argument("-i", "--input", help='The input file specified will be the GTF file provided by snakemake',dest="input", required=True)
parser.add_argument("-o", "--output", help='The output file name and path provided by snakemake',dest="output", required=True)

# this method will run the parser and input the data into the namespace object
args = parser.parse_args()

# if you try running this on the command line and supply it a value for -i or --input
# it will show up here, stored in this object

# try just running this script and supply it a random string for the -i argument
# example: `python parse_gtf.py -i <your_string>`

# replace the code that comes after this with the code necessary to parse the GTF


# Regex solution using dictionary
import re
id_2_name = {}

genename = r'gene_name\s([^;]*)'
geneid = r'gene_id\s([^;]*)'

with open(args.input, 'r') as r:
    for line in r:
        if line.startswith('#'):
            continue

        gene_name = re.search(genename, line)
        gene_id = re.search(geneid, line)

        if gene_id.group().split('"')[1] in id_2_name:
            continue
        else:
            id_2_name[gene_id.group().split('"')[1]] = gene_name.group().split('"')[1]

with open(args.output, 'wt') as w:
    for k, v in id_2_name.items():
        w.write('{}\t{}\n'.format(k, v))

'''
# Solution using CSV and split

import csv

ids = []
names = []

with open(args.input, 'rt') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        if row[0][0] == '#':
            continue
        else:
            if row[2] == 'gene':
                gene_info = row[8].strip().split('; ')
                splits = [_.split('"') for _ in gene_info]
                for l in splits:
                    if l[0].strip() == 'gene_id':
                        ids.append(l[1])
                    if l[0].strip() == 'gene_name':
                        names.append(l[1])

with open(args.output, 'wt') as w:
    for geneid, genename in zip(ids, names):
        w.write('{}\t{}\n'.format(geneid, genename))
'''
