# PROMISING (PRioritization Of Modules In Sets of Implicated Non-intersecting Genes)

The PROMISING method finds promising leads in the candidate genesets created by genome-wide association studies. It proritizes the candidates based on their cohesiveness in an underlying protein-protein functional network under the assumption that the *actual* genes will be functionally similar. For example, the causal genes may operate in the same complex or engage in the same pathway.

This is a work in progress. The method embedded in this code will be published soon.

## Building

The promising binary can be built with the standard configure/make process.

```
./configure
make
```

At the moment, this has only been tested on a Linux machine with gcc. 

### Dependencies

A BLAS library is requrired to compute the Regularized Laplacian of the given protein network. A fast version is recommended for networks with more than a thousand nodes.



## Usage

```
promising  [-c] [-m <string>] [-o <string>] [-p <int>] -g <string>
           -s <string> [--] [--version] [-h]


Where: 
   -g <string>,  --groups <string>
     (required)  Groups file

   -s <string>,  --similarities <string>
     (required)  Similarity matrix
	 
   -m <string>,  --method <string>
     Scoring method, SUM, MAX, MAX-CLIQUE, MAX-3SETS, MAX-4SETS

   -o <string>,  --outfile <string>
     Output summary file

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

```

### Required arguments

The program has two required arguments, one defining the gene sets, the other defining a matrix of similarities between all genes. All other arguments/flags are optional.

### Gene sets (`-g`)

GMT files are used to define gene sets. In a GMT file, each line represents one gene set. All items within a line are tab delimited. The first column is the set's name, the second a description, and the next columns are the genes in the set. Note that the gene sets must be disjoint.

Example:
```
Set 1    description    A    B    C
Set 2    description    D    E
Set 3    description    F    G    H    I
```


### Similarity matrix (`-s`)

The similarity matrix defines the similarities among all the genes contained within the gene sets. Matrices are described in a simple, human readable format. The first line has the tab-delimited column/row names (colum and row order must be the same). Then the matrix is given line-by-line. Columns are delmited with a tab. The continuous values may be in scientific notation. Note that matrices are assumed to be symmetric (hence half of the matrix is redundant information-- we intend to use an alternative matrix format in the future). 

Example:
```
A    B    C    D
0.0  0.0  0.0  1.0
0.0  0.0  1.0  0.0
0.0  1.0  0.0  0.0
1.0  0.0  0.0  0.0
```


### Optional arguments

The following arguments may be supplied.

#### Output summary file (`-o`)

Output the results in a tab-delimited format with three columns (locus, gene, score). Example:

```
locus    gene            score
Locus5   FANCE           0.372794
Locus1   FANCC           0.321366
Locus7   FANCD2          0.311326
Locus13  FANCA           0.296143
Locus2   FANCF           0.295335
Locus10  BRCA2           0.267046
Locus10  RP11-298P3.4    0.245363
Locus3   PALB2           0.220358
Locus9   ERCC4           0.198534
Locus11  FANCI           0.191027
Locus12  MYOG            0.189491
```


#### Method (`-m`)

The method used to calculate candidate scores. Possible methods as `SUM`, `MAX`, `MAX-CLIQUE`, `MAX-3SETS`, and `MAX-4SETS`. The default is `MAX-4SETS`.

*SUM*

A candidate's score is the sum of its similarities to all candidates in opposing genesets. The summation is scaled by the size of the given gene set, because sets are not required to contain the same number of candidates.

*MAX*

A candidate's score is the sum of the most similar candidate in each opposing set.

*MAX-CLIQUE*

A candidate's score is calculated by first finding its best partner on each opposing set, then taking the sum of all pairwise similarities among the partners and the candidate.

*MAX-3SETS*

Considering every combination of gene sets of size three that include a particular candidate, the candidate is scored by taking the sum of the *MAX* method (see above) for the candidate in each of the subsets.


*MAX-4SETS*

Considering every combination of gene sets of size four that include a particular candidate, the candidate is scored by taking the sum of the *MAX* method (see above) for the candidate in each of the subsets. Same as *MAX-3SETS*, but considers 4 gene sets at a time.




## Example

An example network has been included in the "data" directory. A few example genesets are found in "data/examples". Trying out the example Fanconi Anemia geneset, 

```
promising -s data/string_reglaplacian_notextmining_network.tsv -g data/examples/fanconi_anemia.gmt -o summary.txt
```

will produce the following output to stdout.

```
Pulling out necessary subnetwork from file.
Scoring genes.

Top 10 genes
Gene	Score
----------------
FANCI	0.1427
FANCE	0.132256
FANCF	0.131841
FANCD2	0.124038
FANCC	0.119533
FANCA	0.113659
UBE2T	0.0783118
BRCA2	0.0736614
PALB2	0.064135
RP11-298P3.4	0.0455477
```

These top hits just happen to be the actual Fanconi Anemia causal genes. Huzzah!

It will also write a slighly longer summary to the file "summary.txt".

Genes can also be scored on an empirical p-value (see upcoming publication). To do so, use the -p flag with the number of iterations. For example,

```
promising -s data/string_reglaplacian_notextmining_network.tsv -g data/examples/fanconi_anemia.gmt -o summary.txt -p 10000
```

will run the method 10,000 times with random sets of genes the same size as what was input. It will then calculate empirical p-values. It will also report adjusted p-values using the Bonferroni method to control for the FWER.



# Calculating Regularized Laplacian kernel on network

The above priotization method assumes the Regularized Laplacian graph kernel has been applied to protein-protein networks. We provide a tool that will apply the kernel on network files.

## Usage

```
USAGE: 

   reglaplacian  [-k <string>] [-a <float>] -o <string> -n <string>
                 [--] [--version] [-h]

Where: 

   -k <string>,  --kernel <string>
     Graph kernel (default is regularize laplacian)

   -a <float>,  --alpha <float>
     Alpha parameter

   -o <string>,  --outfile <string>
     (required)  Output matrix file

   -n <string>,  --network <string>
     (required)  Network file
```


### Input file


Networks are defined with three entries per line: protein1, protein2, and the score. The items in each line are tab delimited.

Example:
```
A    D    1.0
B    C    1.0
```


## License

Copyright © 2016 kca

Distributed under the GNU public license version 3.

