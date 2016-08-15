# PROMISING (PRioritization On Modules In Sets of Implied Non-intersecting Genes)

The PROMISING method finds promising leads in the candidate genesets created by genome-wide association studies. It proritizes the candidates based on their cohesiveness in an underlying protein-protein functional network under the assumption that the *actual* genes will be functionally similar. For example, the causal genes may operate in the same complex or engage in the same pathway.

This is a work in progress. The method embedded in this code will be published soon.

## Building

The promising binary can be built with the standard configure/make process.

./configure
make


At the moment, this has only been tested on a Linux machine with gcc. 

## Usage
```
Brief USAGE: 
   ./bin/promising  [-o <string>] [-p <int>] -g <string> -n <string> [--]
                    [--version] [-h]

```

```
Required arguments:
-n NETWORK_FILE: Network file for prioritization.
-g GENESETS: Genesets to priortize.

Optional arguments:
-p PVAL_ITERATIONS: Calculate empirical p-values with given number of iterations.
-o OUTPUT_FILENAME: Output a full list of genes with scores to the given file.
```

An example network has been included in the "data" directory. A few example genesets are found in "data/examples".

### Examples

Trying out the example Fanconi Anemia geneset, 

```
./bin/promising -n data/string_reglaplacian_notextmining_network.tsv -g data/examples/fanconi_anemia.txt -o summary.txt
```

will produce the following output to stdout.

```
Output
Top 10 genes
Gene	Score
----------------
FANCI	0.546999
FANCD2	0.526318
FANCE	0.521776
FANCF	0.502351
FANCA	0.481958
FANCC	0.48177
PALB2	0.425826
BRCA2	0.419089
UBE2T	0.406832
RP11-298P3.4	0.32322
```

These top hits just happen to be the actual Fanconi Anemia causal genes. Huzzah!

It will also write a slighly longer summary to the file "summary.txt".

Genes can also be scored on an empirical p-value (see upcoming publication). To do so, use the -p flag with the number of iterations. For example,

```
./bin/promising -n data/string_reglaplacian_notextmining_network.tsv -g data/examples/fanconi_anemia.txt -o summary.txt -p 10000
```

will run the method 10,000 times with random sets of genes the same size as what was input. It will then calculate empirical p-values. It will also report adjusted p-values using the Bonferroni method to control for the FWER.



## License

Copyright © 2016 kca

Distributed under the GNU public license version 3.

