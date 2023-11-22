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

1. For each pair of input sequences, a distance is computed using a modified version of ScoreDist (Sohnhammer & Hollich, 2005). Input sequences may represent a multiple alignment, OR, if the input sequences are unaligned, pairwise protein sequence alignments are computed using a local alignment algorithm (Smith & Waterman, 1981) with affine gap penalties. The Blosum62 amino-acid substitution score matrix is used to determine alignment match scores.

2. The distance matrix is embedded into orthogonal coordinates using metric matrix distance geometry (for example Crippen & Havel, 1988).

3. A nearest-neighbor joining (bifurcating) tree is computed in the orthogonal coordinate space. This is an iterative process in which two (nearest) nodes are replaced by one node at their weighted average position, etc., until a single root node is reached.

4. Starting with the root node of the tree, points within each left and right sub-branch are independently re-embedded (step 2) and re-joined (step 3). The procedure is repeated recursively, gradually removing deleterious effects caused by large-but-inaccurate distances that arise between sequences sharing low homology.

### Program Author
Garry Paul Gippert, DTU Bioengineering, Lyngby, Sealand, Denmark. GarryG 'at' dtu 'dot' dk

### Provenance
Key software and conceptual elements used in Aclust were developed by Garry Paul Gippert while employed at Novozymes A/S, Denmark. The software was kindly relicensed from Novozymes A/S back to Garry Paul Gippert, in 2022, under conditions that it not be commercialized. Aclust is made available under GNU General Public License v3.0.
