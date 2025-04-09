# Program: Aclust
(Alignment and clustering of protein sequences)

ACLUST is a C-language program that generates a phylogenetic tree from protein sequence Fasta input.

### Program input

One or more protein sequence files in FASTA format. If input records represent a multiple alignment use the -maln flag.
Otherwise aclust will run all-vs-all pairwise alignments.

### Program output

One or more trees in Newick format. Alignments in JSON-adjacent and/or Text formats. Distance matrix in text format.
Scattered diagnostics are printed variously to stdout and stderr.

### Output files

All output files share a common prefix (-p this)

Default output files:

* prefix.aln.js		pairwise alignments in JSON-like format
* prefix.dmx.txt	pairwise scoredistances
* prefix.dree.txt	distance NNJ tree in Newick format

Optional output files:

* prefix.aln.txt	pairwise alignments Text (non-standard) format (default OFF, can be enabled)
* prefix.tree0.txt	tree in Newick format, computed from NNJ in the space of embedded orthogonal coordinates
* prefix.tree.txt		tree in Newick format, computed using recursive NNJ/embed

### Inner workings

Pairwise alignments are interpolated from a user-supplied multiple sequence alignment, or computed using a modified local alignment algorithm (based on Smith & Waterman 1981).
The unpublished modification to SW includes a gap-crossover allowance (for gaps that start in one sequence and end in the other).

By convention BLOSUM62 substitution scores are used.

Pairwise distances are based on a modified ScoreDist function (Sonnhammer & Hollich 2005). The unpublished modification normalizes the expectation score to
sequence length rather than alignment length. By convention the shorter sequence is used.

The first tree is computed using a Nearest Neighbor Joining (NNJ) algorithm in the space of all-vs-all pairwise sequence distances. Distances are recomputed using either
branch-averaging leaf-leaf distance averaing may be used as determined by a command line parameter -dave.

A second tree is (optionally) is based on NNJ in the space of 20-dimensional orthogonal coordinates obtained by embedding the distance matrix using distance geomery (Crippen & Havel 1988).

A third tree is (optionally and more slowly) obtained by a recursive application of NNJ and re-embedding to each subbranch of the initial embed tree.

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
