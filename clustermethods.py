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

from scipy.cluster.hierarchy import linkage, cut_tree
from scipy.spatial.distance import squareform as sf
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.cluster.bicluster import SpectralBiclustering
from sklearn.preprocessing import QuantileTransformer, StandardScaler

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

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

#####  Read and format input  ###########################################

# gene by sample table
# Pandas dataframes interpret header (df.columns.values or list(df)) and rownames (df.index)
# this is harder to do with a direct numpy ndarray
print >>sys.stderr, "reading", args.inputfile
df=pd.read_table(args.inputfile, sep='\t', header=0, index_col=0)

print >>sys.stderr, "filtering and normalizing samples..."
# randomly select 100 samples (comment this out if you want to do a full scale run)
ids = list(df)
np.random.shuffle(ids)
df = df.loc[:,ids[:100]]

# remove genes with average expression of 1 fpkm or less (this removes ~45000 genes, keeping 15,000)
df = df[df.mean(1)>1]

# (quantile) normalize - it actually makes no difference for finding the most variable genes
#qt = StandardScaler(copy=False)  # this is another method
qt = QuantileTransformer(random_state=0, copy=False)
qt.fit_transform(df.transpose())

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

# keep track of cluster assignments and silhouette scores
labels = dict()
silscore = dict()

print >>sys.stderr, "Clustering..."

## HDBSCAN
# this method allows samples to not be assigned to a cluster at all
# I can't really make it work yet with anything > min_cluster_size=2, even though that
# does output more than 2 clusters.
clusterer = hdbscan.HDBSCAN(min_cluster_size=2)
labels['hdbscan'] = clusterer.fit_predict(sbg)
# remove unclustered samples from array
keepme = [i for i, j in enumerate(labels['hdbscan'].tolist()) if j != -1]
reduce=sbg.ix[keepme]
silscore['hdbscan'] = silhouette_score(reduce, labels['hdbscan'][keepme])

## BICLUSTERING
# This clusters samples and genes at the same time, supposedly allowing less influence
# of genes with really similar expresson patterns. I haven't played with the parameters yet
n_clusters = (args.clusters, 5)
biclust = SpectralBiclustering(n_clusters=n_clusters, method='log', random_state=0)
biclust.fit(sbg)
labels['biclust']=biclust.row_labels_
silscore['biclust'] = silhouette_score(sbg, labels['biclust'])

## HIERARCHICAL CLUSTERING, PEARSON, COMPLETE LINKAGE
# depends on a distance matrix, which can be created in a variety of ways. Have not tried others yet
s_distance = pairwise_distances(sbg, metric='correlation')
hier_clust_samples =linkage(sf(s_distance), 'complete')
sample_groups = cut_tree(hier_clust_samples, n_clusters=args.clusters)
labels['hier_pearson_complete'] = sample_groups.flatten()
silscore['hier_pearson_complete'] = silhouette_score(s_distance, labels['hier_pearson_complete'], metric="precomputed")

## HIERARCHICAL CLUSTERING, PEARSON, AVERAGE LINKAGE
avg_samples =linkage(sf(s_distance), 'average')
avg_sample_groups = cut_tree(avg_samples, n_clusters=args.clusters)
labels['hier_pearson_avg'] = avg_sample_groups.flatten()
silscore['hier_pearson_avg'] = silhouette_score(s_distance, labels['hier_pearson_avg'], metric="precomputed")

## HIERARCHICAL CLUSTERING, SPEARMAN, AVERAGE LINKAGE
# Rank based clustering seems to work a bit better according to Jaskowiak paper and Max Hauessler.
spearmat, pval = spearmanr(sbg, axis=1)
spearmat = np.round(1-spearmat, 6) # without rounding there are non-zero values on the diagonals
spearman_samples =linkage(sf(spearmat), 'average')
spearman_sample_groups = cut_tree(spearman_samples, n_clusters=args.clusters)
labels['hier_spearman_avg'] = spearman_sample_groups.flatten()
silscore['hier_spearman_avg'] = silhouette_score(spearmat, labels['hier_spearman_avg'], metric="precomputed")

## KMEANS
# This is probabilistic and doesn't give the same outcome twice.
# I thought max_iter was supposed to result in convergence but it doesn't seem to help much
kmeans = KMeans(init='k-means++', n_clusters=args.clusters, n_init=args.clusters, max_iter=100)
labels['kmeans'] = kmeans.fit_predict(sbg)
silscore['kmeans'] = silhouette_score(sbg, labels['kmeans'])

## KMEDOIDS
# Same issue as Kmeans. It is supposed to be better for RNASeq than Kmeans but I haven't seen 
# much to convince me
medoids, clusterinfo, labels['kmed'] = kmedoids.kMedoids(s_distance, args.clusters)
silscore['kmed'] = silhouette_score(s_distance, labels['kmed'])


#####  Remove skewed results  ######################################################

# If the vast majority of samples (currently 80%) end up in the same cluster, something's wrong.
rmlist = []
for id, l in labels.iteritems():
    retval = id
    for c in set(l):
        ct = list(l).count(c)
        if (float(ct)/len(l)) > 0.8:
            print "{} cluster {} occurs {} times, skipping this method".format(id, str(c), ct)
            rmlist.append(id)
            break

for id in rmlist:
    del(silscore[id])
    del(labels[id])

#####  Validate results  ######################################################

for method, score in silscore.iteritems():
    print "{} silhouette score for {}".format(score, method)

# Methods often agree on core groups of samples. This part of the code is not robust, it currently
# depends on commenting out lines to see the effects of removing one or more of the methods.
# it will also die horribly if the skew filter removed the labels dictionary entry for the method.
setlist = []
#setlist.append(robust.list_duplicates(labels['kmeans']))
#setlist.append(robust.list_duplicates(labels['kmed']))
setlist.append(robust.list_duplicates(labels['hier_spearman_avg']))
setlist.append(robust.list_duplicates(labels['hier_pearson_avg']))
setlist.append(robust.list_duplicates(labels['hier_pearson_complete']))
#setlist.append(robust.list_duplicates(labels['hdbscan']))
setlist.append(robust.list_duplicates(labels['biclust']))

# print consistent groups of samples
robust.persistent_groups(copy.copy(setlist), list(sbg.index))



