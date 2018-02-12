import sys
import numpy as np
from scipy.misc import comb

"""
This set of functions turns a list of repeating numbers (such as cluster assignments) into sets of  
positions that always co-occur.
for instance in [0, 1, 2, 3, 1] vs ["a", "bra", "ca", "da", "bra"], positions 2 and 5 occur together (as 2 in 
the first list and "bra" in the second). The lists can contain anything, however -1 values are ignored.
"""

from collections import defaultdict
from sklearn.metrics import adjusted_rand_score
import copy

def persistent_groups(setlist, samples, minsize):
    """Find samples that always group together in groups of minsize"""
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
    return retlist


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

def adj_rand_score(labeldict):
    """Orphan function to calculate adjusted rand score"""
    dictcopy = copy.copy(labeldict)
    for a_label in labeldict:
        del(dictcopy[a_label])
        for b_label in dictcopy:
            print "{0:0.2f}\t{1} vs {2}".format(adjusted_rand_score(labeldict[a_label], dictcopy[b_label]), a_label, b_label)

class cmethod():
    def __init__(self, name, labels, silscore, fract):
        self.name = name
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

