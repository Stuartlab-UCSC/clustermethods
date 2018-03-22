# Clustering_methods
Code to run different clustering methods on a dataset and compare the results.

THIS IS A WORK IN PROGRESS

The goal of this program is to run several clustering methods on a gene by sample FPKM tsv file, and to return information that allows
the researcher to decide on the 'best' one. Unfortunately there are few 'blind' measures to decide this, at least without having
a ground truth. The silhouette score is one, and a good range is supposed to be upwards of 0.6 (it goes to 1).

The current output consists of these silhouette scores and a list of samples that group together in each algorithm.

The program is dependent on three homemade libraries: utils, kmedoids and hdbscan. These should be in the same github repo
as the current program.

The module python-sklearn can be installed on Ubuntu using apt-get; hdbscan, scipy and numpy can be pip installed.
