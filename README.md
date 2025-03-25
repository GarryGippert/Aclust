# Program: Aclust
(Alignment and clustering of protein sequences)

Contents: A C-language program **aclust.c** that generates phylogenetic trees from Fasta protein sequence input.

Garry Paul Gippert, Bioengineering, Danish Technical University, Lyngby, Denmark. MANUSCRIPT IN PREPARATION

Briefly:

- A matrix of pairwise distances is computed from sequence alignments. Pairwise alignments are written in txt and JSON format.

- A Nearest Neighor Joining (NNJ) tree is computed from the distance matrix and written in Newick format.

- The distance matrix is embedded into orthogonal coordinates and NNJ used to construct a second tree, written Newick format.

- The initial embed tree is refined by recursive sub-branch reembedding and NNJ.

### Usage:
Compile the program: <pre>cd src; make; make install</pre>

Run the program: <pre>bin/aclust -s dat/BLOSUM62.dat my.fa</pre>

Help with command line parameters: <pre>bin/aclust -h</pre>

### Help with command line parameters
<b><pre>bin/aclust -h</pre></b>
<pre>
aclust [parameters_and_flags] my.fasta [another.fasta ...] 

Control aligment:
	-s <string>		location of file containing BLOSUM62 substitution score matrix
				default '../dat/BLOSUM62.dat'

	-go <double>		gap open penalty
				default 12.0

	-ge <double>		gap extend penalty
				default 1.0

	-gx <integer>		gap maximum crossover length
				default 100, set to 0 to deactivate crossover allowance

	-m			flag MSA = infer pairwise alignments from input MSA rather than SW alignment
				default off

	-nonself		deactivates self alignments  (obscure usecase)

Control tree generation:
	-e <char>		choose D for only distance tree, S for additional first embed tree, F for additional full recursive embed tree
				default 'F'

	-d <integer>		embed dimension = number of orthogonal coordinates to derive from distance matrix before NNJ
				default 20

Control output:
	-p <string>		prefix for output files
				default = name of first fasta input file

	-j			deactivates writing of JSON alignment file
				default on

	-v			activates more verbose output (usually on stderr)
				default off, can be invoked more than once
</pre>

#### Input files
Protein sequences are read from one or more Fasta files.
If the -m flag is given, input is assumed to be already in a multiple alignment.

#### Output files
All output files share a common prefix (-p option).

<pfx>.aln.txt	Computed or derived pairwise alignments in non-standard TXT format.

<pfx>.aln.js	Computed or derived pairwise alignments in JSON format. Individual lines should be parsed separately.

<pfx>.dmx.txt	Distance matrix file in simple 3-column format. Can in principle be used to derive SSN.

<pfx>.dree.txt	Newick format tree from NNJ on distance matrix

<pfx>.tree0.txt	Newick format tree from NNJ on initial embedding of distances into orthogonal coordinate space.

<pfx>.tree.txt	Newick format tree from recursive reembedding/NNJ starting from initial tree.

### More details

MANUSCRIPT IN PREPARATION

1. For each pair of input sequences I and J, a distance DsubI,J is computed using a renormalized version of ScoreDist (Sonnhammer & Hollich, 2005). By default, alignments are computed using a local alignment algorithm (Smith & Waterman, 1981), affine gap penalties, and Blosum62 (Henikoff & Henikoff, 1992) amino-acid substitution scores. Input an MSA Fasta file (-m flag) and inferred pairwise alignment scores and derived distances are computed using the same method. <i>Scoredist normalization is part of the subject material in the coming paper.</i> The distance matrix is written to (shared_prefix implied) **_dmx.txt**.

2. The distance matrix is embedded into orthogonal coordinates using metric matrix distance geometry (Crippen & Havel, 1988). Each sequence is represented by a point in this space <i>(x,y,z,....) or (w0,w1,w2,...,wM-1), for dimension M usually 20</i>.
   
(2a). A distance-space nearest-neighbor joining tree is constructed from the initial distance matrix and written to **_dree.txt**. Branch lengths to parent are computed using weighted branch counts.
'''

3. A nearest-neighbor joining algorithm in orthogonal coordinate space is used to compute a tree. In this algorithm, sequence points are leaf nodes. Thereafter, the two nearest leaf or internal nodes are replaced by one node at their weighted average position, and the process repeated until a single (the root) node remains.  The initial embed tree is written to **_tree0.txt**.

4. Tree refinement. Starting with the initial root node, trees for left- and right sub-branches are recomputed independently, considering only the subset of points they hold. This process is repeated recursively until the entire tree has been traversed, a computationally expensive refinement. the tree topology by gradually shedding deleterious effects introduced by non-metric distances arising predominantly from (often many) low-homology comparisons within the input sequence pool.

#### Caveats
To understand and explore:

The Metric Matrix Distance Geometry algorithm at the heart of this work was brought into this field out of curiosity. With this algorithm the question 'What does **sequence space** look like?' turned into 'What about embedding a distance matrix derived from pairwise alignments and put it up in 3D?' and from there to 'How to generate Phylogenetic trees automatically?'.

Key drawbacks:
- Compute resources needed to generate sequence alignments.
- Branch distortion and neighbor 'trapping'. <i>Large matrix values dominate the first eigenvectors, preventing correct 'unfurling' of higher-dimensional details - crumpling local tree topology - until one reaches a sufficient level of matrix exhaustion. To make matters worse, large distances are both the least accurate, and the most abundant (most randomly-selected protein pairs are unrelated, but still might return an 'alignment'). Tree refinement - a computationally expensive procedure - ameliorates branch crumpling, but does not overcome left-right branch 'freeze-out' that may place potential neighbors on opposite sides of a division early in the calculation. The latter behaviour has been observed when the sequence pool contains both full-length and short fractional-length proteins.)</i>
- Sequence length limits. 10k amino-acid characters.
- Sequence count limits. 10k sequence entries. This is overkill; due to algorithm complexity the CPU time becomes prohibitive long before this. <i>(Hint: gain experience by starting with a few 10s to few 100s of sequences at first, before gradually moving to larger sequence pools.)</i>

### Program Author
Garry Paul Gippert, Bioengineering, Danish Technical University, Lyngby, Sealand, Denmark. Please contact me with suggestions or questions. Aclust was developed for and with a bunch of the world's great colleagues, starting way back in 2007.

### Provenance and Licensing
Aclust was conceived by Garry P. Gippert while employed at Novozymes
A/S. A derived source code has been made available here under GNU
General Public License v3.0. Date of this GitHub repository Nov 2023.

Additionally where applicable:

This license covers all content within the provided data package.
Other content within the package that e.g. is obtained from external
sources, collaborators or public sources is not covered but should
respect all license restrictions that may be associated.

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software") to utilize, copy, modify, merge, distribute and to
permit persons to whom the Software is furnished to do so, subject
to the following conditions:

THE CONTENT CANNOT BE INTEGRATED OR SOLD FOR COMMERCIAL USE IN ANY
FORM. ALL MATERIAL WITHIN THIS DATA PACKAGE IS FOR NON-COMMCERIAL
USE ONLY.

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

### References

Scoredist : [Sonnhammer & Hollich, 2005](https://pubmed.ncbi.nlm.nih.gov/15857510/).

Local alignment : [Smith & Waterman, 1981](https://pubmed.ncbi.nlm.nih.gov/7265238).

Distance matrix embedding : [Crippen & Havel, 1988](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.540110212) and references therein. <i>(reference to MDS - Multi-Dimensional Scaling - may be appropriate here.)</i>

Blosum62 : [Henikoff & Henikoff, 1992](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC50453/).

GNU public license 3.0 - [License Text](https://www.gnu.org/licenses/gpl-3.0.html#license-text)
