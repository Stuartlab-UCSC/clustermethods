#!/usr/bin/env python2.7

"""

Input FPKM table, output log2(TPM) plus 1

"""

import sys
import argparse
import textwrap
import copy
import numpy as np
import pandas as pd
import robust


#####  get arguments  ###################################################
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\

Input is a gene by sample table in FPKM format. Output is log2(TPM +1), 
suitable for clustering methods.

TPM = FPKM / (sum of FPKM over all genes/transcripts) * 10^6

        '''))
group = parser.add_argument_group('required arguments')
group.add_argument('inputfile', type=str, help='input tsv file')
group.add_argument('outputfile', type=str, help='output tsv file')

# optional flag

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()


df=pd.read_table(args.inputfile, sep='\t', header=0, index_col=0)

# each column is a sample
for column in df:
    dsum = df[column].sum()
    # if a value is 0, replace with 0.0 instead of trying to divide. 
    df[column] = ((df[column]/dsum)*1000000).where(df[column]>0, 0.0)
df = np.log2(df + 1)
df = df.round(2)
df.to_csv(args.outputfile, sep='\t')
