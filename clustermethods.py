#!/usr/bin/env python2.7

"""
THIS IS A WORK IN PROGRESS

The goal of this program is to run several clustering methods on a gene by sample FPKM tsv file, and to return information that allows
the researcher to decide on the 'best' one. Unfortunately there are few 'blind' measures to decide this, at least without having
a ground truth. The silhouette score is one, and a good range is supposed to be upwards of 0.6 (it goes to 1).

The current output consists of these silhouette scores and a list of samples that group together in each algorithm.

The program is dependent on three homemade libraries: robust, kmedoids and hdbscan. These should be in the same github repo
as the current program.

The module python-sklearn can be installed on Ubuntu using apt-get; hdbscan, scipy and numpy can be pip installed.

"""

import sys
import argparse
import textwrap
import copy
import numpy as np
import pandas as pd
import robust
import kmedoids
import hdbscan
import time

from scipy.cluster.hierarchy import linkage, cut_tree
from scipy.spatial.distance import squareform as sf
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.cluster.bicluster import SpectralBiclustering
import sklearn.preprocessing as pp
from sklearn.preprocessing import QuantileTransformer, StandardScaler

def passedTime(start, end):
    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    ptime = "{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds)
    return ptime

#####  get arguments  ###################################################
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
Runs kmeans, kmedoids, hdbscan, biclustering and hierarchical clustering with three different distance matrices
on input gene by sample tsv format file. Data is normalized and filtered for low expression inside the code.

        '''))
group = parser.add_argument_group('required arguments')
group.add_argument('inputfile', type=str, help='input tsv file')
# optional flag
parser.add_argument('--clusters', type=int, default=5,  help="Number of expected clusters (default 5)")
parser.add_argument('--no_subsample', action='store_true',  help="Use the full set (Program uses 100 samples by default)")

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

#####  Read and format input  ###########################################

start = time.time()
# gene by sample table
# Pandas dataframes interpret header (df.columns.values or list(df)) and rownames (df.index)
# this is harder to do with a direct numpy ndarray
print >>sys.stderr, passedTime(start, time.time()), "reading", args.inputfile
df=pd.read_table(args.inputfile, sep='\t', header=0, index_col=0)

print >>sys.stderr, passedTime(start, time.time()),  "filtering and normalizing samples..."
# randomly select 100 samples (comment this out if you want to do a full scale run)
if not args.no_subsample:
    ids = list(df)
    np.random.shuffle(ids)
    df = df.loc[:,ids[:100]]

# remove genes with average expression of 1 fpkm or less (this removes ~45000 genes, keeping 15,000)
df = df[df.mean(1)>1]
print >>sys.stderr, df.shape

# (quantile) normalize - it actually makes no difference for finding the most variable genes
#qt = StandardScaler(copy=False)  # this is another method
qt = QuantileTransformer(random_state=0, copy=False)
#qt.fit_transform(df.transpose())

# Verena's method works too:
#for col in list(df):
#    df[col] = (df[col] - min(df[col]))/(max(df[col])-min(df[col]))

# Rank genes by stdev
colorder = np.argsort(np.std(df,1))
# Sort df by stdev, highest first
df.reindex(index=list(colorder.sort_values(ascending=False).index))
# get 1000 most variable genes over all samples
df = df.head(1000)

### TODO: MUST FILTER OUT LOW SCORING SAMPLES

# Some methods want a gene_by_sample dataframe, others a sample_by_gene
gbs = df
sbg = df.transpose()


#####  Cluster  ######################################################
methods = []

# keep track of cluster assignments and silhouette scores
labels = dict()
silscore = dict()

print >>sys.stderr, passedTime(start, time.time()),  "Clustering..."

print >>sys.stderr, passedTime(start, time.time()),  "HDBSCAN"
## HDBSCAN
# this method allows samples to not be assigned to a cluster at all
# I can't really make it work yet with anything > min_cluster_size=2, even though that
# does output more than 2 clusters.
clusterer = hdbscan.HDBSCAN(min_cluster_size=2)
labels = clusterer.fit_predict(sbg)

# remove unclustered samples from array
keepme = [i for i, j in enumerate(labels.tolist()) if j != -1]
reduce=sbg.ix[keepme]
# how many clusters?
print "HBSCAN cluster count",len(set(labels))
silscore = silhouette_score(reduce, labels[keepme])
hdbObject = robust.cmethod('hdbscan', labels, silscore, args.clusters)
methods.append(hdbObject)

print >>sys.stderr, passedTime(start, time.time()),  "Bicluster"
## BICLUSTERING
# This clusters samples and genes at the same time, supposedly allowing less influence
# of genes with really similar expresson patterns. I haven't played with the parameters yet
n_clusters = (args.clusters, 5)
biclust = SpectralBiclustering(n_clusters=n_clusters, method='log', random_state=0)
biclust.fit(sbg)
labels = biclust.row_labels_
silscore = silhouette_score(sbg, labels)
biclustObject = robust.cmethod('biclust', labels, silscore, args.clusters)
methods.append(biclustObject)

print >>sys.stderr, passedTime(start, time.time()),  "Hier Pearson Complete"
## HIERARCHICAL CLUSTERING, PEARSON, COMPLETE LINKAGE
# depends on a distance matrix, which can be created in a variety of ways. Have not tried others yet
s_distance = pairwise_distances(sbg, metric='correlation')
hier_clust_samples =linkage(sf(s_distance), 'complete')
sample_groups = cut_tree(hier_clust_samples, n_clusters=args.clusters)
labels = sample_groups.flatten()
silscore = silhouette_score(s_distance, labels, metric="precomputed")
hierpcObject = robust.cmethod('hier_pearson_complete', labels, silscore, args.clusters)
methods.append(hierpcObject)

print >>sys.stderr, passedTime(start, time.time()),  "Hier Pearson Average"
## HIERARCHICAL CLUSTERING, PEARSON, AVERAGE LINKAGE
avg_samples =linkage(sf(s_distance), 'average')
sample_groups = cut_tree(avg_samples, n_clusters=args.clusters)
labels = sample_groups.flatten()
silscore = silhouette_score(s_distance, labels, metric="precomputed")
hierpaObject = robust.cmethod('hier_pearson_average', labels, silscore, args.clusters)
methods.append(hierpaObject)

print >>sys.stderr, passedTime(start, time.time()),  "Hier Spearman Average"
## HIERARCHICAL CLUSTERING, SPEARMAN, AVERAGE LINKAGE
# Rank based clustering seems to work a bit better according to Jaskowiak paper and Max Hauessler.
spearmat, pval = spearmanr(sbg, axis=1)
spearmat = np.round(1-spearmat, 6) # without rounding there are non-zero values on the diagonals
spearman_samples =linkage(sf(spearmat), 'average')
sample_groups = cut_tree(spearman_samples, n_clusters=args.clusters)
labels = sample_groups.flatten()
silscore = silhouette_score(spearmat, labels, metric="precomputed")
hiersaObject = robust.cmethod('hier_spearman_average', labels, silscore, args.clusters)
methods.append(hiersaObject)

print >>sys.stderr, passedTime(start, time.time()),  "KMEANS"
## KMEANS
# This is probabilistic and doesn't give the same outcome twice.
# I thought max_iter was supposed to result in convergence but it doesn't seem to help much
# run 10 times, keep the one with best silhouette score
bestsil = 0
for i in xrange(100):
    kmeans = KMeans(init='k-means++', n_clusters=args.clusters, n_init=args.clusters, max_iter=100)
    labels = kmeans.fit_predict(sbg)
    silscore = silhouette_score(sbg, labels)
    if silscore > bestsil:
        kmeansObject = robust.cmethod('kmeans', labels, silscore, args.clusters)
        bestsil = silscore
methods.append(kmeansObject)

print >>sys.stderr, passedTime(start, time.time()),  "KMEDOIDS"
## KMEDOIDS
# Same issue as Kmeans. It is supposed to be better for RNASeq than Kmeans but I haven't seen 
# much to convince me
# run 10 times, keep the one with best silhouette score
bestsil = 0
for i in xrange(100):
    medoids, clusterinfo, labels = kmedoids.kMedoids(s_distance, args.clusters)
    silscore = silhouette_score(s_distance, labels)
    if silscore > bestsil:
        kmedObject = robust.cmethod('kmedoids', labels, silscore, args.clusters)
        bestsil = silscore
methods.append(kmedObject)


#####  Validate results  ######################################################

for i in methods:
    if i.ok:
        print "{0:.2f} silhouette score for {1}".format(i.silscore, i.name)


minsize = 3
print >>sys.stderr, passedTime(start, time.time()),  "Finding consistent groups"
setlist = [i.dups for i in methods if i.ok]

# print consistent groups of samples
grouplist1 = robust.persistent_groups(copy.copy(setlist), list(sbg.index), minsize)
for i in xrange(min(args.clusters, len(grouplist1))):
    print "Group {}, size {}: {}".format(i, len(grouplist1[i])," ".join(str(x) for x in grouplist1[i]))


####  Hold out #####################################################################
# if we hold out one method in the group calculation, do we get larger groups?
# not sure how to do this other than brute force it

#membercount = dict()
bestgroup = grouplist1
maxsize = 0
namelist = [i.name for i in methods if i.ok]
for i in xrange(len(namelist)):
    newset = setlist[::]
    del(newset[i])
    count = 0
#    membercount[namelist[i]] = 0
    grouplist = robust.persistent_groups(newset[::], list(sbg.index), minsize)
    # count size of first groups
    for x in xrange(min(len(grouplist), args.clusters)):
         count += len(grouplist[x])
    # allow a bit of wiggle
    if count > maxsize +5:
        print "best group currently the one without {}, groupsize is then {}".format(namelist[i], count)
        bestgroup = grouplist
        maxsize = count
    print i, namelist[i], count

for i in xrange(min(args.clusters, len(grouplist1))):
    print "Group {}, size {}: {}".format(i, len(bestgroup[i])," ".join(str(x) for x in bestgroup[i]))





