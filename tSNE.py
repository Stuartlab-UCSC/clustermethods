#!/usr/bin/python

# TEMPORARY CODE 
# play with t-SNE settings and plotting output

import sys, os, re, getopt, copy
import numpy as np
import pandas as pd

from sklearn.manifold import TSNE

# gene by sample table
# this interprets header (df.columns.values or list(df)) and rownames (df.index)
# this is harder to do with a direct numpy ndarray
df=pd.read_table('/home/jeltje/Downloads/tpm_set1.tsv', sep='\t', header=0, index_col=0)
sbg=pd.DataFrame.transpose(df)


# create range of tSNE outputs

#for pp in [15, 50]:
#    for ea in [4,12,20,30]: 
#        for lr in [50, 200, 500, 900]:
pp=15
ea=12
ct=0
lr=900
for lr in [200,500,900]:
    ct = 0
    for i in xrange(5):
        ct+=1
        print lr, ct
        tsne = TSNE(perplexity=pp, early_exaggeration=ea, learning_rate=lr)
        tsne.fit(sbg)
        fname = "~/Downloads/tSNE/tsne_lr{0}_kl{1}.{2}.tsv".format(lr, int(round(tsne.kl_divergence_*1000)), ct)
        print fname
        xys = pd.DataFrame(tsne.embedding_, index=sbg.index, columns=['x', 'y'])
        xys.to_csv(fname, sep="\t", index_label='#ID')





#ct = 0
#    print "tSNE", ct
#    ax = plotlist[ct]
#    projection = TSNE(perplexity=50, learning_rate=lr, random_state=2).fit_transform(sbg)
#    ax.set_title("Learning rate=%d" % lr)
#    ax.scatter(*projection.T, s=20, color=colors, linewidth=0, alpha=0.7)
#    ct += 1
#plt.show()

#outname="./gene_dendro.png"
#plt.savefig(outname)


#### Colored scatterplot(s)

#mkeys = mcolors.cnames.keys()
#shuffle(mkeys)
#colors = [mkeys[i] for i in labels]

#(fig, subplots) = plt.subplots(2, 2, figsize=(15, 8))
#plotlist= [subplots[0,0], subplots[0,1], subplots[1,0], subplots[1,1]]

#ct = 0
#for perplexity in [15, 50]:
#for ea in [4,12,20,30]:   # Makes almost no difference, leave at default 12
#for lr in [50, 200, 500, 900]:
#    print "tSNE", ct
#    ax = plotlist[ct]
#    projection = TSNE(perplexity=50, learning_rate=lr, random_state=2).fit_transform(sbg)
#    ax.set_title("Learning rate=%d" % lr)
#    ax.scatter(*projection.T, s=20, color=colors, linewidth=0, alpha=0.7)
#    ct += 1
#plt.show()

#outname="./gene_dendro.png"
#plt.savefig(outname)

