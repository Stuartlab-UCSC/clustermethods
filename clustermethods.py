#!/usr/bin/env python2.7

"""
The goal of this program is to run several clustering methods on a gene by sample log2(FPKM/TPM)+1 tsv file, and to return information that allows
the researcher to decide on the 'best' one. Unfortunately there are few 'blind' measures to decide this, at least without having
a ground truth. The silhouette score is one, and a good range is supposed to be upwards of 0.6 (it goes to 1).
Another one is the cophenetic index, but it can only be calculated with the use of a distance matrix so it cannot be used
for biclustering or kmeans/kmedoids.

The program is dependent on three homemade libraries: utils, kmedoids and hierarchical. These should be in the same github repo
as the current program. Hierarchical is a slightly altered copy of sklearn.cluster.AgglomerativeClustering. See the code for details.

The module python-sklearn can be installed on Ubuntu using apt-get; pandas, scipy and numpy can be pip installed.

"""

import sys
import argparse
import textwrap
import copy
import numpy as np
import pandas as pd
import utils
import kmedoids
import time

from scipy.cluster.hierarchy import linkage, cut_tree, cophenet
from scipy.spatial.distance import squareform as sf, pdist
from scipy.stats import spearmanr
from sklearn.cluster import KMeans #, AgglomerativeClustering
from hierarchical import AgglomerativeClustering
from sklearn.metrics import silhouette_score
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.cluster.bicluster import SpectralBiclustering


def passedTime(start, end):
    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    ptime = "{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds)
    return ptime

#####  get arguments  ###################################################
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
Runs kmeans, kmedoids, biclustering and different kinds of hierarchical clustering
on input gene by sample tsv format file. Data is filtered for low expression and only the top 1000
most variable genes are used.

Program calculates the cophenetic index (for hierarchical trees) and the silhouette score. 
Note that the cophenetic index will not change for different cluster numbers: it is calculated for the whole tree

Outputs:
<args.base>.clusters.tsv	a sample by clustermethod table which can be used as 'color' upload to Tumormap
<args.base>.scores.txt		silhouette score and cophenetic index for methods
<args.base>.groups.txt		(optional) groupings of samples that nearly always cluster together no matter which method is used.

NOTE: Several methods used here are probabilistic and will NOT give the exact same output twice. 
However, usually you will find that the same method scores best on different runs.


        '''))
group = parser.add_argument_group('required arguments')
group.add_argument('inputfile', type=str, help='input tsv file')
# optional flag
parser.add_argument('--clusters',      type=int, default=5,  help="Number of expected clusters (default 5)")
parser.add_argument('--genes',         action='store_true',  help="Cluster genes instead of samples")
parser.add_argument('--maxfract',      type=float, default=0.8,  help="Maximum percentage of input samples allowed to end up in the same cluster (default 0.8)")
parser.add_argument('--subsample',     type=int, default=False,  
    help="Use this number of samples instead of the whole set. Good for running quick tests on larger datasets")
parser.add_argument('--base',          type=str, default='clusters', help="Basename for output (default clusters)")
parser.add_argument('--print_reduced', action='store_true',  
     help="Print the 1000 genes that were used for clustering in a gene by sample table, useful if you want to do any subsequent runs, or upload to Tumormap")
parser.add_argument('--print_groups',  action='store_true',  
     help="Save a list of sample groups that always occur together, no matter which clustering method was used. This is useful for defining core groups of samples. Setting this flag also adds a 'shared' column to the clusters file")
parser.add_argument('--min_groupsize', type=int, default=3,  help="Minimum number of samples in a shared group (default 3)")
parser.add_argument('--noK',           action='store_true',  help="Skip Kmeans and Kmedoids (useful if you care about the grouped samples)")

args = parser.parse_args()
if args.maxfract > 1 or args.maxfract < 0.05:
    print >>sys.stderr, "please select a maxfract between 0.05 and 1"
    sys.exit(1)

#####  Read and format input  ###########################################

start = time.time()
# gene by sample table
# Pandas dataframes interpret header (df.columns.values or list(df)) and rownames (df.index)
print >>sys.stderr, passedTime(start, time.time()), "reading", args.inputfile
df=pd.read_table(args.inputfile, sep='\t', header=0, index_col=0)

if args.subsample:
    # randomly select samples 
    if len(list(df)) > args.subsample:
        print >>sys.stderr, passedTime(start, time.time()),  "randomly selecting {} samples of {}".format(args.subsample, len(list(df)))
        ids = list(df)
        np.random.shuffle(ids)
        df = df.loc[:,ids[:args.subsample]]

# remove genes with average expression of 1 or less 
print >>sys.stderr, passedTime(start, time.time()),  "removing low expression genes (<1)"
df = df[df.mean(1)>1]
sbg = df

if not args.genes:
    # Rank genes by stdev
    print >>sys.stderr, passedTime(start, time.time()),  "selecting 1000 genes with highest stdev"
    df['stdev'] = df.std(axis=1)
    df.sort_values('stdev', inplace=True, ascending=False)
    df.drop('stdev', axis=1, inplace=True)
    # get 1000 most variable genes over all samples
    df = df.head(1000)
    
    if args.print_reduced:
        df.to_csv(args.base+'.1kgenes.tsv', sep='\t')

    sbg = df.transpose()


#####  Cluster  ######################################################

methods = []	# holds all info on each clustering method
n_clusters = args.clusters


print >>sys.stderr, passedTime(start, time.time()),  "Bicluster (probabilistic, YMMV)"
# This clusters samples and genes at the same time. It can work extremely well but is sensitive to its parameters
# In testing, n_components and n_best did well when set to args.clusters, at least in a largeish dataset (1000 samples)
biclust = SpectralBiclustering(n_clusters=n_clusters, method='log', n_components=n_clusters, n_best=n_clusters, svd_method='randomized')
biclust.fit(sbg)
labels = biclust.row_labels_
silscore = silhouette_score(sbg, labels)
biclustObject = utils.cmethod('biclust', labels, silscore, 0.0, args.maxfract)
methods.append(biclustObject)

print >>sys.stderr, passedTime(start, time.time()),  "Hierarchical clustering with a distance matrix"
# using a precomputed distance matrix (correlation and yule work equally well in testing)
# average and complete linkage both work well
s_distance = pairwise_distances(sbg, metric='correlation')
for l in ['average', 'complete']:
    ac = AgglomerativeClustering(n_clusters=n_clusters, affinity='precomputed', linkage=l)
    ac.fit_predict(s_distance)
    labels = ac.labels_
    # the cophenetic index can only be calculated on a scipy.hierarchy distance matrix, so create one.
    ac = utils.distanceMatrix(ac)
    # NOTE: the cophenetic score does not change if you alter the number of clusters
    c, coph_dists = cophenet(ac.distancematrix, sf(s_distance))
    silscore = silhouette_score(s_distance, labels, metric="precomputed")
    cname = 'hier_' + l
    hierObject = utils.cmethod(cname, labels, silscore, c, args.maxfract)
    methods.append(hierObject)

print >>sys.stderr, passedTime(start, time.time()),  "Hierarchical clustering using WARD linkage"
# Ward calculates its own euclidean distance matrix
ac = AgglomerativeClustering(n_clusters=n_clusters, linkage='ward')
ac.fit_predict(sbg)
labels = ac.labels_
ac = utils.distanceMatrix(ac)
# create euclidean distance matrix for cophenetic score
c, coph_dists = cophenet(ac.distancematrix, pdist(sbg))
silscore = silhouette_score(sbg, labels)
hierObject = utils.cmethod('hier_ward', labels, silscore, c, args.maxfract)
methods.append(hierObject)

print >>sys.stderr, passedTime(start, time.time()),  "Hierarchical clustering using spearman rank"
# unfortunately this is not an option in sklearn so use scipy
spearmat, pval = spearmanr(sbg, axis=1)
spearmat = np.round(1-spearmat, 6) # without rounding there are non-zero values on the diagonals
for l in ['average', 'complete']:
    spearman_samples =linkage(sf(spearmat), l)
    sample_groups = cut_tree(spearman_samples, n_clusters=args.clusters)
    labels = sample_groups.flatten()
    c, coph_dists = cophenet(spearman_samples, sf(spearmat))
    silscore = silhouette_score(spearmat, labels, metric="precomputed")
    cname = 'hier_spearman_' + l
    hiersaObject = utils.cmethod(cname, labels, silscore, c, args.maxfract)
    methods.append(hiersaObject)

if not args.noK:
    print >>sys.stderr, passedTime(start, time.time()),  "KMEANS  (probabilistic, YMMV)"
    # I thought max_iter was supposed to result in convergence but it doesn't seem to help much
    # run 5 times just to show range of outcome
    for i in xrange(1,6):
        kmeans = KMeans(init='k-means++', n_clusters=args.clusters, n_init=args.clusters, max_iter=100)
        labels = kmeans.fit_predict(sbg)
        silscore = silhouette_score(sbg, labels)
        cname = 'kmeans_' + str(i)
        kmeansObject = utils.cmethod(cname, labels, silscore, 0.0, args.maxfract)
        methods.append(kmeansObject)
    
    print >>sys.stderr, passedTime(start, time.time()),  "KMEDOIDS  (probabilistic, YMMV)"
    # Same issue as Kmeans. Same approach
    for i in xrange(1,6):
        medoids, clusterinfo, labels = kmedoids.kMedoids(s_distance, args.clusters)
        silscore = silhouette_score(s_distance, labels)
        cname = 'kmedoids_' + str(i)
        kmedObject = utils.cmethod(cname, labels, silscore, 0.0, args.maxfract)
        methods.append(kmedObject)


#####  Create consistent sample groups  ######################################################
# This outputs a list of samples that always occur together in a cluster, no matter which method is used
# It also adds a 'shared' column to the clusters output file. An attempt is made to give similar clusters
# similar labels, so that cluster B1 largely contains the same samples in each cluster method.

if args.print_groups:
    print >>sys.stderr, passedTime(start, time.time()),  "Finding consistent groups in all methods used"
    setlist = [i.dups for i in methods if i.ok]
    grouplist1, tally = utils.persistent_groups(copy.copy(setlist), list(sbg.index), args.min_groupsize)
    # print the groups to file
    samples = list(sbg.index)
    labels = [-1] * len(sbg.index)
    
    with open(args.base + '.groups.txt', 'w') as o:
        o.write("Group\tMembercount\n")
        for i in xrange(len(grouplist1)):
            clabel = i+1
            for x in grouplist1[i]:
                idx = samples.index(x)
                labels[idx] = clabel
            o.write("Group {}\t{}\t{}\n".format(i, len(grouplist1[i]),"\t".join(str(x) for x in grouplist1[i])))
    
    shareObject = utils.cmethod('shared', labels, 0.0, 0.0, args.maxfract)
    shareObject.ok = False
    methods.append(shareObject)
    utils.rename_labels(methods, tally)
    shareObject.orderedlabels = ['' if x=='0' else x for x in shareObject.orderedlabels]

with open(args.base + '.scores.txt', 'w') as s:
    s.write("silhouette score\tcophenetic score\tdensity score\tmethod\n")
    for i in methods:
        if i.ok and not (i.name == 'shared'):
            if i.cscore == 0.0:
                s.write("{0:.2f}\tNA\t{1}\n".format(i.silscore, i.name))
            else:
                s.write("{0:.2f}\t{1:0.2f}\t{2}\n".format(i.silscore, i.cscore, i.name))


#####  Output cluster file  ###################################################
#samples on rows, clusters in columns - prepend a letter to the labels so TumorMap will show them properly

outdf = pd.DataFrame([i.orderedlabels for i in methods if i.ok])
outdf.columns = df.columns
outdf = outdf.transpose()
outdf.columns = [i.name for i in methods if i.ok]
outdf.to_csv(args.base+'.clusters.tsv', sep='\t')


