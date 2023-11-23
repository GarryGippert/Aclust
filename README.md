# Program: Aclust
(Alignment and clustering of protein sequences)

Contents: A single C-language source program **aclust.c** that generates phylogenetic trees from Fasta protein sequence input.

Garry Paul Gippert, MANUSCRIPT IN PREPARATION

Briefly:

- A matrix of pairwise distances is computed from sequence alignments.

- A binary tree is computed directly from the distance matrix using NNJ (Nearest Neighor Joining) and count-weighted distance averaging when computing branch lengths to parent.

- A second binary tree is computed, also using NNJ, but in the space of orthogonal coordinates obtained by eigenvalue decomposition of the metric matrix derived from the distance matrix (as if it were Euclidean, which in most cases is far from).

- Recursive sub-branch reembedding refines topology and branch lengths, but (drawback) freezes early tree left-right subdivisions.

### Usage:
Compile the program: <pre>cd src; make; make install</pre>

Run the program: <pre>bin/aclust -s dat/BLOSUM62.dat my.fa</pre>

Help with command line parameters: <pre>bin/aclust -h</pre>

### Help with command line parameters
<b><pre>bin/aclust -h</pre></b>
<pre>
aclust [parameters_and_flags] my.fasta [another.fasta ...] 

Required parameters:
	-s 'path'		file location of substitution score matrix (could be 'dat/BLOSUM62.txt')
Optional parameters:
	-p 'string'		prefix for all output files (default=name of first input fasta file)
	-d integer		embed dimension (default 20)
Optional flags:
	-m 			activates to interpret input Fasta as MSA
Less important flags:
	-j			deactivates writing of JSON alignment file
	-nonself		deactivates self alignments
	-v			activates more verbose output
</pre>

#### Input and output files
Input Fasta
- Sequences may be pre-aligned (MSA). Otherwise local (SW) alignments are computed.

Output Newick trees, alignments, and distance matrix (sharing a common prefix)
- _aln.txt- Alignments (text)
- _aln.js - Alignments (JSON-parsable text)
- _dmx.txt - Distance Matrix (text)
- **_dree.txt** - Distance Matrix tree (newick)
- **_tree0.txt** - Initial Embed tree (newick)
- **_tree.txt** - Refined Embed tree (newick)

### Bit longer introduction

MANUSCRIPT IN PREPARATION

1. For each pair of input sequences I and J, a distance DsubI,J is computed using a renormalized version of ScoreDist (Sohnhammer & Hollich, 2005). By default, alignments are computed using a local alignment algorithm (Smith & Waterman, 1981), affine gap penalties, and Blosum62 (Henikoff & Henikoff, 1992) amino-acid substitution scores. Input an MSA Fasta file (-m flag) and inferred pairwise alignment scores and derived distances are computed using the same method. <i>Scoredist normalization is part of the subject material in the coming paper.</i> The distance matrix is written to (shared_prefix implied) **_dmx.txt**.

2. The distance matrix is embedded into orthogonal coordinates using metric matrix distance geometry (Crippen & Havel, 1988). Each sequence is represented by a point in this space <i>(x,y,z,....) or (w0,w1,w2,...,wM-1), for dimension M usually 20</i>.
   
(2a). A distance-space nearest-neighbor joining tree is constructed from the initial distance matrix and written to **_dree.txt**. Branch lengths to parent are computed using weighted branch counts.
'''

3. A nearest-neighbor joining algorithm in orthogonal coordinate space is used to compute a tree. In this algorithm, sequence points are leaf nodes. Thereafter, the two nearest leaf or internal nodes are replaced by one node at their weighted average position, and the process repeated until a single (the root) node remains.  The initial embed tree is written to **_tree0.txt**.

4. Tree refinement. Starting with the initial root node, trees for left- and right sub-branches are recomputed independently, considering only the subset of points they hold. This process is repeated recursively until the entire tree has been traversed, a computationally expensive refinement. the tree topology by gradually shedding deleterious effects introduced by non-metric distances arising predominantly from (often many) low-homology comparisons within the input sequence pool.

### Program Author
Garry Paul Gippert, Bioengineering, Danish Technical University, Lyngby, Sealand, Denmark. Please contact me with suggestions or questions. Aclust was developed for and with a bunch of the world's great colleagues, starting way back in 2007.

### Provenance and Licensing
Key software and conceptual elements used in Aclust were developed by Garry Paul Gippert while employed at Novozymes A/S, Denmark. The software was kindly relicensed by Novozymes A/S back to Garry Paul Gippert, in Jan 2022, with the licencing text given below. The current Aclust repository contains a cherry-picked and shrink-wrapped selection from that material, therefore the same licensing conditions are likely to apply here.
<pre>**LICENSE.txt** 
Copyright 2022 Novozymes A/S

This license covers all content within the provide data 
package, that originates from Novozymes A/S and delivered as 
'Garry-Paul-Gippert_data-extract_1.tar.gz' on 31.01.2022. 
Other content within the package that e.g. is obtained from 
external sources, collaborators or public sources is not 
covered but should respect all license restrictions that may 
be associated. 

Permission is hereby granted, free of charge, to any person 
obtaining a copy of this software and associated documentation 
files (the "Software") to utilize, copy, modify, merge, 
distribute and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions:

THE CONTENT CANNOT BE INTEGRATED OR SOLD FOR COMMERCIAL USE IN 
ANY FORM. ALL MATERIAL WITHIN THIS DATA PACKAGE ORIGINATING 
FROM NOVOZYMES A/S IS FOR NON-COMMCERIAL USE ONLY.

The above copyright notice and this permission notice shall be 
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
DEALINGS IN THE SOFTWARE.
</pre>

The spelling mistake is not mine. To the extent it does not conflict with the above conditions Aclust is made available here under GNU General Public License v3.0. Date of this GitHub repository Nov 2023.

### References

Scoredist : [Sonnhammer & Hollich, 2005](https://pubmed.ncbi.nlm.nih.gov/15857510/).

Local alignment : [Smith & Waterman, 1981](https://pubmed.ncbi.nlm.nih.gov/7265238).

Distance matrix embedding : [Crippen & Havel, 1988](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.540110212) and references therein. <i>(reference to MDS - Multi-Dimensional Scaling - may be appropriate here.)</i>

Blosum62 : [Henikoff & Henikoff, 1992](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC50453/).

GNU public license 3.0 - [License Text](https://www.gnu.org/licenses/gpl-3.0.html#license-text)
