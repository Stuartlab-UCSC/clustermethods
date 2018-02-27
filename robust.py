import sys
import numpy as np
import pandas as pd
import scipy.spatial.distance as dist
from scipy.misc import comb

"""
Helper functions for clustermethods.py
"""

from collections import defaultdict
from sklearn.metrics import adjusted_rand_score
import copy

def persistent_groups(setlist, samples, minsize):
    """Find samples that always group together in groups of minsize. For instance 
       in [0, 1, 2, 3, 1] vs ["a", "bra", "ca", "da", "bra"], positions 2 and 5 occur together (as 2 in 
       the first list and "bra" in the second). The lists can contain anything, however -1 values are ignored.
    """
    tally = common_idx(setlist.pop(0), setlist.pop(0))
    while(setlist):
        tally = common_idx(tally, setlist.pop(0))
    #Match sample IDs to their indices; each set in setlist is a group
    retlist = []
    for myset in tally:
        slist = sorted([samples[i] for i in myset])
        if len(slist) < minsize:
            continue
	retlist.append(slist)
    # sort so that the largest groups are listed first
    retlist.sort(key=len, reverse=True)
    return retlist, tally



def list_duplicates(mylist):
    """Identify sets of positions with the same value"""
    tally = defaultdict(set)
    for i in xrange(len(mylist)):
        tally[mylist[i]].add(i)
    list_of_sets = []
    for key, locs in tally.items():
        if (len(locs)>1 and key >=0):
            list_of_sets.append(locs)
    return list_of_sets


def common_idx(X, Y):
    """Take lists of numbers, return sets with common numbers"""
    tally = []
    for xset in X:
        for yset in Y:
            have_common = xset.intersection(yset)
            if len(have_common) >1:
                tally.append(have_common)
    return tally

def adj_rand_score(truth, labels):
    """Orphan function to calculate adjusted rand score"""
    return adjusted_rand_score(truth, labels)

def adj_rand_score_old(labeldict):
    """Orphan function to calculate adjusted rand score"""
    dictcopy = copy.copy(labeldict)
    for a_label in labeldict:
        del(dictcopy[a_label])
        for b_label in dictcopy:
            print "{0:0.2f}\t{1} vs {2}".format(adjusted_rand_score(labeldict[a_label], dictcopy[b_label]), a_label, b_label)

class cmethod():
    def __init__(self, name, labels, silscore, cscore, fract):
        self.name = name
        self.cscore = cscore
        self.silscore = silscore
	self.labels = labels
        self.strlabels = ['A'+str(i) if i>=0 else '0' for i in labels]
        self.dups = list_duplicates(self.labels)
        self.skewCheck(fract)
    def skewCheck(self, fract):
        self.ok = True
        # If more than <fract> fraction of samples end up in the same cluster, something's wrong.
        for c in set(self.labels):
            ct = list(self.labels).count(c)
            if (float(ct)/len(self.labels)) > fract:
                print >>sys.stderr, "{} cluster {} occurs {} times, skipping this method".format(self.name, str(c), ct)
                self.ok = False
                break

def rename_labels(methodlist, tally):
    """Attempt to give shared groups the same labels in all methods"""
    # persistent_tally is a list of sets of indices to the (str)labels
    # order tally by length of list
    tally.sort(key=len, reverse=True)
    ct = 0
    for m in methodlist:
        m.replaced = list(set(m.strlabels))
        m.orderedlabels = m.strlabels[:]
    for group in tally:
        ct += 1
        newlabel = 'B'+str(ct)
        # loop over all cmethods
        for m in methodlist:
            labelmatches = [m.strlabels[g] for g in group]
            # get the most occurring label match (this is a problem with ties!)
            label = max(set(labelmatches), key=labelmatches.count)
            # only replace once
            if label in m.replaced:
                if label != '0':
                    m.orderedlabels = [newlabel if x==label else x for x in m.orderedlabels]
                m.replaced.remove(label)

def randScore(truth, labels):
    """Unadjusted Rand Index from  https://stats.stackexchange.com/questions/89030/rand-index-calculation answer 3"""
    truth = listBin(truth)
    labels = listBin(labels)
    tp_plus_fp = comb(np.bincount(truth), 2).sum()
    tp_plus_fn = comb(np.bincount(labels), 2).sum()
    A = np.c_[(truth, labels)]
    tp = sum(comb(np.bincount(A[A[:, 0] == i, 1]), 2).sum()
             for i in set(truth))
    fp = tp_plus_fp - tp
    fn = tp_plus_fn - tp
    tn = comb(len(A), 2) - tp - fp - fn
    return (tp + tn) / (tp + fp + fn + tn)


def listBin(mylist):
    """Turn a list of strings into a list of numbers:  ["A", "A", "B", "A"] becomes [0, 0, 1, 0]"""
    uniqs = list(set(mylist))
    return [uniqs.index(i) for i in mylist]

def distanceMatrix(clusterObject):
    """Add a scipy hierarchy style distance matrix to the input clusterobject"""
    # what is the lowest number in the input children (usually 0 or 1)?
    minct = min(min(x, y) for x, y in clusterObject.children_)
    childcount = {key:1 for key in xrange(minct, len(clusterObject.children_)+1)}
    counters = []
    for s, e in clusterObject.children_:
        ct = childcount[s] + childcount[e]
        counters.append(ct)
        childcount[len(childcount)] = ct
    clusterObject.distancematrix = np.column_stack([clusterObject.children_, clusterObject.distance, counters])
    return clusterObject

# the functions below were all copied from code in TumorMap's https://github.com/ucscHexmap/compute/calc/
def ztransDF(df):
    '''
    :param df: pandas structure, series or data frame
    :return: z-score normalized version of df.
    '''
    return ((df - df.mean(axis=0)) / df.std(axis=0, ddof=0))

def spatialWeightMatrix(xys):
    '''
    :param xys: x-y positions for nodes on the map
    :return: col X col inverse euclidean distance matrix,
    row normalized to sum to 1.
    '''
    invDist = inverseEucDistance(xys)

    # Self comparisons to 0
    notZd = pd.DataFrame(invDist,index=xys.index,columns=xys.index)

    return (notZd / notZd.sum(axis=1))

def inverseEucDistance(xys):

    distmat = dist.squareform(dist.pdist(xys,'euclidean'))
    return 1 / (1 + distmat)

def catSSS(catLee):
    #below we are twice over counting (matrix is symetric) but also twice over dividing, so that's not a mistake

    #         average of on_diagonal                           average of off diagonal
#    return (catLee.trace() / catLee.shape[0]) - ((catLee.sum() - catLee.trace())/ (catLee.shape[1]**2 - catLee.shape[1]))
    return [(catLee.trace() / catLee.shape[0]), ((catLee.sum() - catLee.trace())/ (catLee.shape[1]**2 - catLee.shape[1]))]


# code below copied from calc/leesL.py
def leesL(Z, V):
    VTV = np.dot(V.transpose(), V)
    ZTVTVZ = np.dot(np.dot(Z.transpose(),  VTV), Z)
    return ZTVTVZ / VTV.sum().sum()

def leesLScore(xys, labels, split=False):
    """
    Returns the distance score between a set of coordinates and features (here cluster labels)
    """
    datMat = pd.get_dummies(labels).apply(pd.to_numeric)
    datMat = datMat.loc[xys.index]
    # in case we lost any categories for the map we are looking at
    datMat = datMat[datMat.columns[datMat.sum() != 0]] 
    datMat = ztransDF(datMat)
    wm = spatialWeightMatrix(xys)
    [posScore, negScore]= catSSS(
              leesL(
                 datMat, wm
              )
            )
    if split:
        return [posScore, negScore]
    else:
        return posScore+negScore

