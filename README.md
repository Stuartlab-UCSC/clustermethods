# Clustering methods
The goal of this program is to run several clustering methods on a gene by sample expression file in `tsv` format, and to return a table with cluster labels
for each sample and each method used. A silhousette score is also calculated and output in a `.scores.txt` file

## Dependencies
```
sklearn
pandas
scipy
numpy
```
The module `python-sklearn` can be installed on Ubuntu using `apt-get`; `pandas`, `scipy` and `numpy` can be pip installed.

## Before you start
### 1. Determine number of clusters

Most clustering methods require a preset number of clusters. Unfortunately there's no way to automatically determine the optimum number of clusters so you will have to do some preliminary analysis.

The best way is to visualize your data using a dimensionality reduction method such as [t-SNE](https://github.com/Stuartlab-UCSC/clustermethods/blob/master/tSNE.py) or [TumorMap](http://tumormap.ucsc.edu/). Simply input your [batch corrected](https://sysbiowiki.soe.ucsc.edu/node/323) gene by sample table and look at the resulting figure.

If your data comes out as a single blob on either of these methods, you likely won't be able to get multiple clusters from the clustering algorithms. Of course you can request three clusters and the programs will give you something, but it's likely meaningless and will differ between methods.

### 2. Log2(TPM/FPKM+1)

If you haven't already done so, convert your FPKM input to log2(TPM+1). The fpkmToLog2TPMplus1.py in this repository does exactly what it says.

## Running the program
### 3. Clustering methods

The program clustermethods.py in this repository takes your gene by sample table and your required number of clusters. First, it removes the genes with low expression and then it selects the top 1000 most variable genes. Using fewer genes [tends to lead to better outcomes](https://www.ncbi.nlm.nih.gov/pubmed/28778489). It then runs the following clustering methods:

-   [Spectral biclustering](https://www.ncbi.nlm.nih.gov/pubmed/12671006). This method clusters genes and samples at the same time.
-   Hierarchical clustering with a correlation based distance matrix and average linkage
-   Hierarchical clustering with a correlation based distance matrix and complete linkage
-   Hierarchical clustering with ward linkage (creates a euclidean distance matrix)
-   Hierarchical clustering with spearman rank and average linkage
-   Hierarchical clustering with spearman rank and complete linkage
-   Kmeans
-   Kmedoids

If you want to cluster genes instead of samples, set the `--genes` flag. The program will now run on the filtered gene set, without selecting the 1000 most variable genes.

###  4. The workflow

 - After reading in the input table, the program selects a subset of samples if you have set `--subsample`. 
 - It then removes any genes with an average expression of 1. Remember that this is log2 of (TPM plus 1). If you want to use for instance raw FPKM, you need to change this number in the code.
 - It calculates the standard deviation of expression of all genes, then orders them from most to least variable.
 - The top 1000 genes are selected and used as input.
 - Input to the various cluster methods is either this table, or a pairwise distance correlation matrix that is calculated during the run.

### 5. Scoring

For each method, the program caculates the [silhouette score](https://en.wikipedia.org/wiki/Silhouette_(clustering)) and for the hierarchical methods it also calculates the [cophenetic index](https://en.wikipedia.org/wiki/Cophenetic_correlation) (which can only be calculated on dendrograms). These scores say something about how well the method fits the data, but note that the cophenetic index will not change if you alter the cluster number. This is because it gets calculated on the full tree, not on the number of clusters you request (which is a step that happens afterwards).

### 6. Output
```
<args.base>.clusters.tsv 
<args.base>.scores.txt 
<args.base>.1kgenes.tsv (optional)
<args.base>.groups.txt (optional)
```
To help you visualize the results, the program outputs a [Tumormap](http://tumormap.ucsc.edu/)-ready input feature map (`clusters.tsv`) with labels for each clustering method. The `scores.txt` file contains silhouette and (where possible) cophenetic index scores.
You can also opt to output the reduced gene set that was used for clustering and use this as Tumormap's layout input or as input to another iteration to the clustermethods.py program. In both programs this will increase upload and calculation speed.

You can opt to get a list of grouped samples. These are samples that end up in the same cluster in all the methods used in the program.

### 7. Making changes

The clustermethod.py code and the associated libraries are well commented and you should feel free to check out a copy and play with it.  
Most of the clustering methods used are from [Python's scikit learn package](http://scikit-learn.org/stable/), which has extensive documentation and example code. You should be able to add additional methods from this package, or from [scipy clustering](https://docs.scipy.org/doc/scipy/reference/cluster.html).


