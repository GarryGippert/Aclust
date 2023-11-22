# Program: Aclust
(Alignment and clustering of protein sequences)

The program generate a phylogenetic tree from Fasta protein sequence input.

Input Fasta file
- Sequences may be pre-aligned (MSA). Otherwise local (SW) alignments are computed.

Output Newick file (**_tree.txt**) *output files do share a common prefix
- Derived or computed pairwise alignments are written to local files (**_aln.txt** and/or **_aln.js**) enriched with scoredist values.
- Derived distance matrix is written to local file (**_dmx.txt**).
.
### Brief introduction

MANUSCRIPT IN PREPARATION

1. For each pair of input sequences, a distance is computed using a modified version of ScoreDist (Sohnhammer & Hollich, 2005). Input sequences may represent a multiple alignment, OR, if the input sequences are unaligned, pairwise protein sequence alignments are computed using a local alignment algorithm (Smith & Waterman, 1981) with affine gap penalties. Alignment scores are computed using Blosum62 (Henikoff & Henikoff, 1992) amino-acid substitution scores in either case.

2. The distance matrix is embedded into orthogonal coordinates using metric matrix distance geometry (Crippen & Havel, 1988). Each sequence is represented by a point in this space <i>(x,y,z,w0,w1,w2,...)</i>.

3. A nearest-neighbor joining algorithm in orthogonal coordinate space is used to compute a tree. In this algorithm, sequence points are leaf nodes. Thereafter, the two nearest leaf or internal nodes are replaced by one node at their weighted average position, and the process repeated until a single (the root) node remains.

4. Starting with the initial root node, trees for left- and right sub-branches are recomputed independently, considering only the subset of points they hold. This process - repeated recursively until the entire tree has been traversed - is computationally expensive, but improves the tree topology by gradually shedding deleterious effects introduced by non-metric distances arising predominantly from (often many) low-homology comparisons within the input sequence pool.

### Program Author
Garry Paul Gippert, Bioengineering, Danish Technical University, Lyngby, Sealand, Denmark. Please contact me with suggestions or questions.

### Provenance and Licensing
Key software and conceptual elements used in Aclust were developed by Garry Paul Gippert while employed at Novozymes A/S, Denmark. The software was kindly relicensed by Novozymes A/S back to Garry Paul Gippert, in Jan 2022, under conditions that it not be commercialized. The text of that agreement is reproduced here:
<pre>LICENSE.txt 
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

To the extent it does not conflict with the above conditions, Aclust is made available under GNU General Public License v3.0. Date of this GitHub repository Nov 2023.

### References

Scoredist : [Sonnhammer & Hollich, 2005](https://pubmed.ncbi.nlm.nih.gov/15857510/).

Local alignment : [Smith & Waterman, 1981](https://pubmed.ncbi.nlm.nih.gov/7265238).

Distance matrix embedding : [Crippen & Havel, 1988](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.540110212) and references therein. <i>(reference to MDS - Multi-Dimensional Scaling - may be appropriate here.)</i>

Blosum62 : [Henikoff & Henikoff, 1992](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC50453/).

GNU public license 3.0 - [License Text](https://www.gnu.org/licenses/gpl-3.0.html#license-text)
