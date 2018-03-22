#!/usr/bin/env python2.7

"""
Given a table with cluster labels (as output by clustermethods.py) and a table with xy coordinates (as output by t-SNE), 
calculate leesL score fore each column in cluster labels

"""

import sys
import argparse
import textwrap
import copy
import numpy as np
import pandas as pd
import utils

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\

        '''))
group = parser.add_argument_group('required arguments')
#group.add_argument('clusterfile', type=str, help='input tsv file')
# optional flag
parser.add_argument('--clusters', type=str, required=True,  help="cluster labels file")
parser.add_argument('--xys', type=str, required=True,  help="t-SNE output xy table")

args = parser.parse_args()

#####  Read and format input  ###########################################

labeldf=pd.read_table(args.clusters, sep='\t', header=0, index_col=0)
xys=pd.read_table(args.xys, sep='\t', header=0, index_col=0)
#labeldf=pd.read_table('~/Downloads/tmp.clusters.tsv', sep='\t', header=0, index_col=0)
#xys=pd.read_table("~/Downloads/tSNE/default.tsv", sep='\t', header=0, index_col=0)

for c in labeldf.columns:
    labels = labeldf[c]
    print "{1:.03f}\t{2}\t{0}".format(c, utils.leesLScore(xys, labels), args.xys)
