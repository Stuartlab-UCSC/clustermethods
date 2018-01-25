"""
This set of functions turns a list of repeating numbers (such as cluster assignments) into sets of  
positions that always co-occur.
for instance in [0, 1, 2, 3, 1] vs ["a", "bra", "ca", "da", "bra"], positions 2 and 5 occur together (as 2 in 
the first list and "bra" in the second). The lists can contain anything, however -1 values are ignored.
"""

from collections import defaultdict
from sklearn.metrics import adjusted_rand_score
import copy

def persistent_groups(setlist, samples):
    """Find samples that always group together"""
    tally = common_idx(setlist.pop(0), setlist.pop(0))
    while(setlist):
        tally = common_idx(tally, setlist.pop(0))
    #Match sample IDs to their indices; each set in setlist is a group
    for myset in tally:
        dada = sorted([samples[i] for i in myset])
        print "Group:",len(dada)," ".join(str(x) for x in dada)

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


