# Program: Aclust
(Alignment and clustering of protein sequences)

ACLUST is a C-language program that generates a phylogenetic tree from protein sequence Fasta input.

A binary tree is based on a NNJ (Nearest Neighbor Joining) algorithm in the space of all-vs-all pairwise sequence distances computed using
a modified version of ScoreDist (Sonnhammer & Hollich 2005).

Pairwise alignments are computed from local alignments (Smith & Waterman 1981) or interpolated from a user-supplied multiple sequence aligment.

Additional trees may be built in orthogonal coordinate space after embedding the distance matrix using metric matrix
distance geometry (Crippen & Havel 1988).

### Program input

One or more FASTA files. To input Fasta records in a multiple alignment use the -maln flag. Otherwise aclust
will run all-vs-all pairwise alignments.

### Output

The program produces some messages on stdout and stderr. In addition (with a common prefix):

prefix.aln.js	pairwise alignments JSON-like format (default ON, can be disabled)

prefix.aln.txt	pairwise alignments Text (non-standard) format (default OFF, can be enabled)

prefix.dmx.txt	pairwise scoredistances (default ON)

prefix.dree.txt	tree in Newick format, computed using NNJ in distance space

optional:

prefix.tree0.txt	tree in Newick format, computed from NNJ in the space of embedded orthogonal coordinates

prefix.tree.txt		tree in Newick format, computed using recursive NNJ/embed

### Usage:

Download the code: <pre>git clone git@github.com:GarryGippert/Aclust.git</pre>

<pre>cd Aclust</pre>

Compile aclust program (requires gcc or similar C-compiler): <pre>cd src; make</pre>

Run the program: <pre>aclust -s ../dat/BLOSUM62.dat my.fa</pre>

If you move or run the program from another location, you must refer to the correct
path to the substitution score matrix file Aclust/dat/BLOSUM62.dat using the -s parameter.

<pre>aclust -s somepath/Aclust/dat/BLOSUM62.dat my.fasta</pre>

Additional help is available using the -h flag <pre>bin/aclust -h</pre>

### References

Scoredist : [Sonnhammer & Hollich, 2005](https://pubmed.ncbi.nlm.nih.gov/15857510/).

Local alignment : [Smith & Waterman, 1981](https://pubmed.ncbi.nlm.nih.gov/7265238).

Distance matrix embedding : [Crippen & Havel, 1988](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.540110212) and references therein. <i>(reference to MDS - Multi-Dimensional Scaling - may be appropriate here.)</i>

Blosum62 : [Henikoff & Henikoff, 1992](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC50453/).

### Author

Garry Paul Gippert, Bioengineering, Danish Technical University, Lyngby, Denmark. MANUSCRIPT IN PREPARATION

### Licensing

GNU public license 3.0 - [License Text](https://www.gnu.org/licenses/gpl-3.0.html#license-text)
