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

## File formats

GMT files are used to define gene sets. In a GMT file, each line represents one gene set. All items within a line are tab delimited. The first column is the set's name, the second a description, and the next columns are the genes in the set.

Example:
```
Set 1    description    A    B    C
Set 2    description    D    E
Set 3    description    F    G    H    I
```

Networks are defined with three entries per line: protein1, protein2, and the score. The items in each line are tab delimited.

Example:
```
A    D    1.0
B    C    1.0
```

Matrices are described in a simple, human readable format. The first line has the tab-delimited column names. The the matrix is given line-by-line, each continuous value delmited with a tab. The continuous values may be in scientific notation.

Example:
```
A    B    C    D
0.0  0.0  0.0  1.0
0.0  0.0  1.0  0.0
0.0  1.0  0.0  0.0
1.0  0.0  0.0  0.0
```

## Usage
```
USAGE: 

   promising  [-d <string>] [-m <string>] [-s <int>] [-o
              <string>] [-p <int>] -g <string> [--]
              [--version] [-h]

Where: 

   -o <string>,  --outfile <string>
     Output summary file

   -p <int>,  --pval <int>
     Pvalue iterations

   -g <string>,  --groups <string>
     (required)  Groups file

   -m <string>,  --matrix <string>
     (required)  Matrix file
	 
   -d <string>,  --degree_groups <string>
     Degree groups for node permutations

```

An example network has been included in the "data" directory. A few example genesets are found in "data/examples".

### Examples

Trying out the example Fanconi Anemia geneset, 

```
promising -n data/string_reglaplacian_notextmining_network.tsv -g data/examples/fanconi_anemia.gmt -o summary.txt
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
promising -n data/string_reglaplacian_notextmining_network.tsv -g data/examples/fanconi_anemia.gmt -o summary.txt -p 10000
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

## License

Copyright © 2016 kca

Distributed under the GNU public license version 3.

