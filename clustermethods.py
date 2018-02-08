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
from sklearn.cluster import KMeans, AgglomerativeClustering
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
Runs kmeans, kmedoids, biclustering and different kinds of hierarchical clustering
on input gene by sample tsv format file. Data is filtered for low expression and only the top 1000
most variable genes are used.

Outputs:
<args.base>.clusters.tsv	a sample by clustermethod table which can be used as 'color' upload to Tumormap
<args.base>.shared.txt		groupings of samples that nearly always cluster together no matter which method is used.

        '''))
group = parser.add_argument_group('required arguments')
group.add_argument('inputfile', type=str, help='input tsv file')
# optional flag
parser.add_argument('--clusters', type=int, default=5,  help="Number of expected clusters (default 5)")
parser.add_argument('--maxfract', type=float, default=0.8,  help="Maximum percentage of input samples allowed to end up in the same cluster (default 0.8)")
parser.add_argument('--no_subsample', action='store_true',  help="Use the full set (Program uses 1000 samples by default)")
parser.add_argument('--base', type=str, default='clusters', help="Basename for output (default clusters)")
parser.add_argument('--min_groupsize', type=int, default=3,  help="Minimum number of samples in a shared group (default 3)")

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
if args.maxfract > 1 or args.maxfract < 0.05:
    print >>sys.stderr, "please select a maxfract below between 0.05 and 1"
    sys.exit(1)

#####  Read and format input  ###########################################

start = time.time()
# gene by sample table
# Pandas dataframes interpret header (df.columns.values or list(df)) and rownames (df.index)
print >>sys.stderr, passedTime(start, time.time()), "reading", args.inputfile
df=pd.read_table(args.inputfile, sep='\t', header=0, index_col=0)

# randomly select 1000 samples (comment this out if you want to do a full scale run)
if not args.no_subsample:
    print >>sys.stderr, passedTime(start, time.time()),  "randomly selecting 1000 samples"
    ids = list(df)
    np.random.shuffle(ids)
    df = df.loc[:,ids[:1000]]

# remove genes with average expression of 1 or less 

print >>sys.stderr, passedTime(start, time.time()),  "removing low expression genes (<1)"
df = df[df.mean(1)>1]

# Rank genes by stdev
print >>sys.stderr, passedTime(start, time.time()),  "selecting 1000 genes with highest stdev"
colorder = np.argsort(np.std(df,1))
# Sort df by stdev, highest first
df.reindex(index=list(colorder.sort_values(ascending=False).index))
# get 1000 most variable genes over all samples
df = df.head(1000)

#df.to_csv('red.tsv', sep='\t')
#sys.exit()

### TODO: MUST FILTER OUT LOW SCORING SAMPLES

# Some methods want a gene_by_sample dataframe, others a sample_by_gene
gbs = df
sbg = df.transpose()


#####  Cluster  ######################################################

methods = []	# holds all info on each clustering method
n_clusters = args.clusters

print >>sys.stderr, passedTime(start, time.time()),  "Clustering..."


print >>sys.stderr, passedTime(start, time.time()),  "Bicluster (probabilistic, YMMV)"
# This clusters samples and genes at the same time. It can work extremely well but is sensitive to its parameters
# In testing, n_components and n_best did well when set to args.clusters, at least in a largeish dataset (1000 samples)
biclust = SpectralBiclustering(n_clusters=n_clusters, method='log', n_components=n_clusters, n_best=n_clusters, svd_method='randomized')
biclust.fit(sbg)
labels = biclust.row_labels_
silscore = silhouette_score(sbg, labels)
biclustObject = robust.cmethod('biclust', labels, silscore, args.maxfract)
methods.append(biclustObject)


print >>sys.stderr, passedTime(start, time.time()),  "Hierarchical clustering"
# using a precomputed distance matrix (correlation and yule work equally well in testing)
# average and complete linkage both work well
s_distance = pairwise_distances(sbg, metric='correlation')
for l in ['average', 'complete']:
    ac = AgglomerativeClustering(n_clusters=n_clusters, affinity='precomputed', linkage=l)
    ac.fit_predict(s_distance)
    labels = ac.labels_
    silscore = silhouette_score(s_distance, labels, metric="precomputed")
    cname = 'hier_' + l
    hierObject = robust.cmethod(cname, labels, silscore, args.maxfract)
    methods.append(hierObject)

# using WARD clustering
ac = AgglomerativeClustering(n_clusters=n_clusters, linkage='ward')
ac.fit_predict(sbg)
labels = ac.labels_
silscore = silhouette_score(sbg, labels)
hierObject = robust.cmethod('hier_ward', labels, silscore, args.maxfract)
methods.append(hierObject)

# using spearman rank based matrix - unfortunately not an option in sklearn so use scipy
spearmat, pval = spearmanr(sbg, axis=1)
spearmat = np.round(1-spearmat, 6) # without rounding there are non-zero values on the diagonals
for l in ['average', 'complete']:
    spearman_samples =linkage(sf(spearmat), l)
    sample_groups = cut_tree(spearman_samples, n_clusters=args.clusters)
    labels = sample_groups.flatten()
    silscore = silhouette_score(spearmat, labels, metric="precomputed")
    cname = 'hier_spearman_' + l
    hiersaObject = robust.cmethod(cname, labels, silscore, args.maxfract)
    methods.append(hiersaObject)

print >>sys.stderr, passedTime(start, time.time()),  "KMEANS  (probabilistic, varies wildly)"
# I thought max_iter was supposed to result in convergence but it doesn't seem to help much
# run 5 times just to show range of outcome
for i in xrange(1,6):
    kmeans = KMeans(init='k-means++', n_clusters=args.clusters, n_init=args.clusters, max_iter=100)
    labels = kmeans.fit_predict(sbg)
    silscore = silhouette_score(sbg, labels)
    cname = 'kmeans_' + str(i)
    kmeansObject = robust.cmethod(cname, labels, silscore, args.maxfract)
    methods.append(kmeansObject)

print >>sys.stderr, passedTime(start, time.time()),  "KMEDOIDS  (probabilistic, varies wildly)"
# Same issue as Kmeans. Same approach
for i in xrange(1,6):
    medoids, clusterinfo, labels = kmedoids.kMedoids(s_distance, args.clusters)
    silscore = silhouette_score(s_distance, labels)
    cname = 'kmedoids_' + str(i)
    kmedObject = robust.cmethod(cname, labels, silscore, args.maxfract)
    methods.append(kmedObject)


#####  Validate results  ######################################################

for i in methods:
    if i.ok:
        print "{0:.2f} silhouette score for {1}".format(i.silscore, i.name)



print >>sys.stderr, passedTime(start, time.time()),  "Finding consistent groups in all methods used"
setlist = [i.dups for i in methods if i.ok]

# print consistent groups of samples
grouplist1 = robust.persistent_groups(copy.copy(setlist), list(sbg.index), args.min_groupsize)
for i in xrange(min(args.clusters, len(grouplist1))):
    print "Group {}, size {}: {}".format(i, len(grouplist1[i])," ".join(str(x) for x in grouplist1[i]))


####  Hold out #####################################################################
# To avoid having one method mess up the shared groups, do an all vs all comparison with one holdout each round
# Then keep the least fragmented set

bestgroup = grouplist1
maxsize = 0
namelist = [i.name for i in methods if i.ok]
for i in xrange(len(namelist)):
    newset = setlist[::]
    del(newset[i])
    count = 0
    grouplist = robust.persistent_groups(newset[::], list(sbg.index), args.min_groupsize)
    # count size of first groups
    for x in xrange(min(len(grouplist), args.clusters)):
         count += len(grouplist[x])
    # allow a bit of wiggle
    if count > maxsize +5:
        print >>sys.stderr, "best group currently the one without {}, groupsize is then {}".format(namelist[i], count)
        bestgroup = grouplist
        maxsize = count
    print >>sys.stderr, i, namelist[i], count

# Treat the shared groups like any other cluster method. Any sample not in a group gets a -1
# print the groups to file
samples = list(sbg.index)
labels = [-1] * len(sbg.index)

with open(args.base + '.groups.txt', 'w') as o:
    o.write("Group\tMembercount\n")
    for i in xrange(len(bestgroup)):
        clabel = i+1
        for x in bestgroup[i]:
            idx = samples.index(x)
            labels[idx] = clabel
        #tumtype = [x.split('_')[0] for x in bestgroup[i]]
        o.write("Group {}\t{}\t{}\n".format(i, len(bestgroup[i]),"\t".join(str(x) for x in bestgroup[i])))

shareObject = robust.cmethod('shared', labels, "0.0", args.maxfract)
methods.append(shareObject)

#####  Output cluster file  ###################################################
#samples on rows, clusters in columns - prepend a letter to the labels so TumorMap will show them properly
outdf = pd.DataFrame([i.strlabels for i in methods if i.ok])
outdf.columns = df.columns
outdf = outdf.transpose()
outdf.columns = [i.name for i in methods if i.ok]
outdf.to_csv(args.base+'.clusters.tsv', sep='\t')

#for a in ['euclidean', 'l1', 'l2', 'manhattan', 'cosine']:
#    for l in ['average', 'complete']:
#        ac = AgglomerativeClustering(n_clusters=n_clusters, affinity=a, linkage=l)
#        ac.fit_predict(sbg)
#        labels = ac.labels_
#        ajr = adjusted_rand_score(truth, labels)
#        rs = robust.randScore(truth, labels)
#        silscore = silhouette_score(s_distance, labels, metric="precomputed")
#        print "{0:.2f} AJR\t{4:.02f} AR\t{1:.02f} silscore\tAG {2} {3}".format(ajr, silscore, d, l, rs)