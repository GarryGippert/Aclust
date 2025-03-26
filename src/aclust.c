/* ACLUST

Aclust generates a phylogenetic tree from sequence pairwise distances,
computed either by pre-aligned (multiple aligned) FASTA input, or by
generating all-vs-all pairwise local alignments using SW.

Please contact the author Garry Paul Gippert, GarryG@dtu.dk, DTU Bioengineering,
Danish Technical University, with questions or for more information.

Distance-matrix approach.

NNJ in distance space, and embedded coordinate space.

1. For each pair of input sequences, a distance is computed using a modified version of ScoreDist
(Sonnhammer & Hollich, 2005). Input sequences may represent a multiple alignment, OR, if the input
sequences are unaligned, pairwise protein sequence alignments are computed using a local alignment
algorithm (Smith & Waterman, 1981) with affine gap penalties. The Blosum62 amino-acid substitution
score matrix is used to determine alignment match scores.

2. The distance matrix is embedded into orthogonal coordinates using metric matrix distance geometry
(for example Crippen & Havel, 1988).

3. A nearest-neighbor joining (bifurcating) tree is computed in the orthogonal coordinate space.
This is an iterative process in which two (nearest) nodes are replaced by one node at their weighted
average position, etc., until a single root node is reached.

4. Starting with the root node of the tree, points within each left and right sub-branch are
independently re-embedded (step 2) and re-joined (step 3). The procedure is repeated recursively,
gradually removing deleterious effects caused by large-but-inaccurate distances that arise between
sequences sharing low homology.

Input file(s) contain protein sequences in Fasta format.

Output files produced: (all with shared prefix)
        prefix_aln.txt  Text (somewhat human readable) alignments, useful in trouble-shooting.
        prefix_aln.js   JSON-parsable details of each pairwise alignment, useful to parse later.
        prefix_dmx.txt  Distance matrix in plain text format <labelI> <labelJ> <distanceIJ>.  not needed.
	prefix_dree.txt	Distance matrix NNJ tree, Newick format.
        prefix_tree0.txt Metric matrix embed tree, Newick-format.
        prefix_tree.txt Recursively refined embed tree, Newick-format.

ACLUST was developed and written by Garry Paul Gippert, and packaged as a single C source file in 2023.

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

ACLUST is made available online at GitHub GarryGippert:Aclust
(https://github.com/GarryGippert/Aclust) with a GNU General Public
License v3.0, and may be used, modified and distributed according
to the principles therein, as long as the above license text from
Novozymes is adhered to, and is preserved and distributed within all
copies of this code and derivative works. Signed,
Garry Paul Gippert Nov 23, 2023
*/
/* ACLUST Additional notes:

The gap affine penalty used here allows a gap 'cross over' with no additional gap-open penalty.  A
cross-over gap starts in one alignment string and ends in the other alignment string, with no
intervening match state.

CROSS-OVER GAP:

In the following, O is the gap opening penalty, and e is the gap extension penalty. In the present
work, O = 12, e = 1. The crossover allowance reduces the gap penalty by 
	(2 * O + 4 * e)   -   (O + 5 * e)    =    O - e
compared to twice opening a gap. (Fictitious example for illustration.)

	usual	     OeeOee		normal cost of two opened gaps
	cross	     Oeeeee		with cross-over allowance costs O-e less
	aln1:	ACDEF---KLMNPQRSTV
	aln2:	ACDEFGHI---NPQRSTV

Re-embedding of isolated sub-branches has the effect of gradually reducing deleterious effects of
long-range, inaccurate distances.  Probably both the choice of distance function, and algorithmic
choices such as only taking 20 eigenvalues/vectors, contribute to distortions of local topology when
including long-range distances.  CAVEAT: Nodes have been observed to be 'trapped' in the wrong initial
branch, for example when comparing full-length sequences and sequence fragments.

PADDED AND GAPPED ALIGNMENTS:

The GAP symbol '-' is used to indicate non-matched positions within the local alignment. The PAD symbol
'+' is used to indicate regions of sequence that fall outside the local alignment, and allows the full
input sequences to be reproduced from the output alignment string.

	>aln1
	ACGHIKNPQRVWY
	>aln2
	DEFGHIKLMNPQRST

For example if we align the two sequences above, we get the folling local alignment having 10 aligned
positions (starting with G, ending with R), 8 matched positions including a gap of length 2, and a
total pad length of 20 which accounts for subsequences found outside the local alignment.

	Pair aln1 13 x aln2 15
	AC+++GHIK--NPQRVWY++
	++DEFGHIKLMNPQR+++ST
	Ascore 33 oi 2 o2 3 Plen 20 Alen 10 Mlen 8 Ilen 8 Glen 2 Olen 1 Clen 8 Nlen 8
	Mscore 46.000000 M1 46 M2 46 MR -8 SD0 18.6776 SD1 39.6954 SD2 50.6719 SD 39.6954

ACLUST was developed and written by Garry Paul Gippert, and packaged as a single C source file in 2023.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXSEQUENCELEN 10000
#define MAXFILENAME 1000
#define MAXLINELEN 1000
#define MAXWORDLEN 100
#define DEFAULT_SCORE_MATRIX "../dat/BLOSUM62.dat"

#define ALIGN_GAP_CHAR '-'
#define ALIGN_GAP_VAL -99.9
#define ALIGN_GAP_NDX -1
#define ALIGN_PAD_CHAR '+'
#define ALIGN_PAD_VAL -99.9
#define ALIGN_PAD_NDX -2
#define ALIGN_GAP	0x0001	/* enable GAP */
#define ALIGN_PAD	0x0002	/* enable PAD */
#define ALIGN_CROSS	0x0004	/* enable gap cross-over */
#define ALIGN_BODY	0x0008	/* enable gap body-gaps (not recommended) */


#define	EPSILON		1.0e-8
#define NEARZERO(a)	(fabs(a) < 10.0*EPSILON ? 0.0 : (a))
#define SIGN(a)		( (a) < 0.0 ?   (-1.0) :    (1.0) )

int p_v = 0;			/* verbose flag, set to 1 for additional diagnostic output */
double p_go = 10.0;		/* gap open = first gap penalty */
double p_ge = 0.5;		/* gap extend = next gap penalty */
int p_gx = 100;			/* maximum gap crossover length (set to 0 to deactivate) */
int p_maln = 0;			/* read multiple alignment */
int p_jaln = 1;			/* write alignments as JSON */
int p_taln = 0;			/* write alignment as text */
int p_baln = 0;			/* write alignment as binary */
int p_nonself = 0;		/* do not align with self (show only off-diagonal elements */
int p_metadata = 0;		/* print tree node metadata */

char *scorematrixfile = NULL;
char *f_dmxfilename = NULL;
char *oprefix = NULL;

FILE *jsnfp = NULL;		/* file pointer for writing alignment JSON line-by-line */
FILE *alnfp = NULL;		/* file pointer for writing alignment free text */

void j_opn(FILE * fp)
/* open json line */
{
#ifdef JSONLONG
	fprintf(fp, "[");
#else
	fprintf(fp, "{");
#endif
}

void j_cmt(FILE * fp, char *comment)
/* output jsonified comment */
{
	if (comment)
		fprintf(fp, ", \"comment\": \"%s\"", comment);
}

void j_url(FILE * fp, char *url)
/* output jsonified url */
{
	if (url)
		fprintf(fp, ", \"url\": \"%s\"", url);
}

#define NO 0
#define YES 1

void j_int(FILE * fp, int comma, char *key, int value, char *comment, char *url)
/* output jsonified integer */
{
#ifdef JSONLONG
	fprintf(fp, "{\"key\": \"%s\", \"value\": %d", key, value);
	j_cmt(fp, comment);
	j_url(fp, url);
	fprintf(fp, "}");
#else
	fprintf(fp, "\"%s\": %d", key, value);
#endif
	if (comma == YES)
		fprintf(fp, ", ");
}

void j_str(FILE * fp, int comma, char *key, char *value, char *comment, char *url)
/* output jsonified string */
{
#ifdef JSONLONG
	fprintf(fp, "{\"key\": \"%s\", \"value\": \"%s\"", key, value);
	j_cmt(fp, comment);
	j_url(fp, url);
	fprintf(fp, "}");
#else
	fprintf(fp, "\"%s\": \"%s\"", key, value);
#endif
	if (comma == YES)
		fprintf(fp, ", ");
}

void j_dbl(FILE * fp, int comma, char *key, double value, char *comment, char *url)
/* output jsonified string */
{
#ifdef JSONLONG
	fprintf(fp, "{\"key\": \"%s\", \"value\": %g", key, value);
	j_cmt(fp, comment);
	j_url(fp, url);
	fprintf(fp, "}");
#else
	fprintf(fp, "\"%s\": %g", key, value);
#endif
	if (comma == YES)
		fprintf(fp, ", ");
}

void j_cls(FILE * fp)
/* close JSON line */
{
#ifdef JSONLONG
	fprintf(fp, "]\n");
#else
	fprintf(fp, "}\n");
#endif
}

double *double_vector(int n)
/* return double * vector of length n*/
{
	double *v;
	int i;
	if ((v = (double *)malloc(n * sizeof(double))) == NULL)
		fprintf(stderr, "failed to allocate double * vector\n"), exit(1);
	for (i = 0; i < n; i++)
		v[i] = 0;
	return (v);
}

double *double_vector_copy(int n, double *u)
/* return copy of double * vector */
{
	double *v = double_vector(n);
	memcpy(v, u, n * sizeof(double));
	return (v);
}

extern void double_vector_drand48(int n, double *v, double lower, double upper)
/* assign double *vector random values in range lower to upper */
{
	int i;
	for (i = 0; i < n; i++)
		v[i] = lower + (upper - lower) * drand48();
}

extern void double_vector_normal(int n, double *v)
/* normalize double * vector to unit length */
{
	double s = 0.0;
	int i;
	for (i = 0; i < n; i++)
		s += v[i] * v[i];
	s = sqrt(s);
	for (i = 0; i < n; i++)
		v[i] /= s;
}

void double_vector_free(int n, double *v)
/* free double * vector */
{
	if (v)
		free((char *)v);
	v = NULL;
}

double **double_matrix(int ni, int nj)
/* return double ** matrix having dimensions ni * nj
// but leave second dimension unallocated if nj <= 0 */
{
	double **m;
	int i, j;
	if ((m = (double **)malloc(ni * sizeof(double *))) == NULL)
		fprintf(stderr, "Unable to allocate double ** matrix\n"), exit(1);
	for (i = 0; i < ni; i++)
		m[i] = NULL;
	if (nj > 0)
		for (i = 0; i < ni; i++) {
			m[i] = double_vector(nj);
			for (j = 0; j < nj; j++)
				m[i][j] = 0.0;
		}
	return (m);
}

void double_matrix_free(int ni, int nj, double **m)
/* free double ** matrix */
{
	int i;
	for (i = 0; i < ni; i++)
		double_vector_free(nj, m[i]);
	free((char *)m);
	m = NULL;
}

void double_matrix_print(int ni, int nj, double **m, char **lab)
/* print double ** matrix, with optional labels */
{
	int i, j;
	int print_half = 1;
	for (i = 0; i < ni; i++) {
		if (lab != NULL)
			printf("%20s :", lab[i]);
		for (j = 0; j <= (print_half ? i : nj - 1); j++)
			printf(" %5.1f", (float)m[i][j]);
		printf("\n");
	}
}

char *char_vector(int n)
{
	char *str;
	if ((str = (char *)malloc((n + 1) * sizeof(char))) == NULL)
		fprintf(stderr, "Unable to allocate char * vector\n"), exit(1);
	strcpy(str, "");
	str[n] = 0;
	return (str);
}

char *char_string(char *inp)
{
	int n = strlen(inp);
	char *str = char_vector(n);
	strcpy(str, inp);
	str[n] = 0;
	return (str);
}

char **char_matrix(int ni, int nj)
/* return char ** matrix having dimensions ni * nj
// but leave second dimension unallocated if nj <= 0 */
{
	char **m;
	int i, j;
	if ((m = (char **)malloc(ni * sizeof(char *))) == NULL)
		fprintf(stderr, "Unable to allocate char ** matrix\n"), exit(1);
	for (i = 0; i < ni; i++)
		m[i] = NULL;
	if (nj > 0)
		for (i = 0; i < ni; i++) {
			m[i] = char_vector(nj);
			for (j = 0; j < nj; j++)
				m[i][j] = ' ';
		}
	return (m);
}

void char_matrix_free(int ni, int nj, char **m)
{
	int i;
	for (i = 0; i < ni; i++)
		if (m[i] != NULL)
			free(m[i]);
	m = NULL;
}

int *int_vector(int n)
{
	int *v;
	int i;
	if ((v = (int *)malloc(n * sizeof(int))) == NULL)
		fprintf(stderr, "Unable to allocate int * vector\n"), exit(1);
	for (i = 0; i < n; i++)
		v[i] = 0;
	return (v);
}

int *int_vector_ramp(int n)
{
	int *v = int_vector(n);
	int i;
	for (i = 0; i < n; i++)
		v[i] = i;
	return (v);
}

void int_vector_free(int n, int *v)
{
	free((char *)v);
	v = NULL;
}

int **int_matrix(int ni, int nj)
/* allocate integer matrix having dimensions ni * nj
// but leave second dimension unallocated if nj <= 0 */
{
	int **m;
	int i, j;
	if ((m = (int **)malloc(ni * sizeof(int *))) == NULL)
		fprintf(stderr, "Unable to allocate int ** matrix\n"), exit(1);
	for (i = 0; i < ni; i++)
		m[i] = NULL;
	if (nj > 0)
		for (i = 0; i < ni; i++) {
			m[i] = int_vector(nj);
			for (j = 0; j < nj; j++)
				m[i][j] = 0.0;
		}
	return (m);
}

void int_matrix_free(int ni, int nj, int **m)
{
	int i;
	for (i = 0; i < ni; i++)
		int_vector_free(nj, m[i]);
	free((char *)m);
	m = NULL;
}

/* BLOSUM */

double **blosum_mtx = NULL;
int nb = 0;
#define MAXN 30
char alphabet[MAXN];

double blosum_zscore(double mscore, int mlen)
/* Compute Blosum62 Zscore from average random Blosum score = -mlen, and stddev score = 2 n^0.5,
	for gapless fragments, unclear to use ascore or mscore */
{
        if (mlen <= 0) 
                return (-99.9); 
        return ((mscore + (double) mlen) / (2.0 * sqrt((double) mlen)));
}

static int print_blosum_pscore_parameters = 0;

double blosum_pscore(double zscore, int qlen)
/* Pscore = -log10(Pvalue). from extreme value distrubution.
   Probability of Zscore P(z)
        P(z) = exp[(A-z)/B - exp[(A-z)/B]]/B
        D(z) = exp[- exp[(A-z)/B]]
   where
        A = 3.30  + 0.387  * ln(qlen);
        B = 0.393 + 0.0585 * ln(qlen);
   The probability that zscore >= z is given by 1 - D(z).
   We put this in -log10 units.
   LARGER PSCORE IS BETTER.  GPG 040614
*/
{
        if (qlen < 0)
                return 0.0;

	/* Extreme value distribution constants to three sig figs */
#define EVD_A_CONST 3.30
#define EVD_A_SLOPE 0.387
#define EVD_B_CONST 0.393
#define EVD_B_SLOPE 0.0585
        double A = EVD_A_CONST + EVD_A_SLOPE * log((double) qlen);
        double B = EVD_B_CONST + EVD_B_SLOPE * log((double) qlen);
        double P = exp((A - zscore) / B - exp((A - zscore) / B)) / B;
        double D = exp(-exp((A - zscore) / B));

        /* limit of pscore resolution */
#define PSCORE_LIMIT 1.0e-15
        double f = 1.0 - D;
        f = (f > PSCORE_LIMIT ? f : PSCORE_LIMIT);
        double S = -log10(f);
#ifdef DEBUG
        printf("# blosum_pscore Z %f Q %d A %g B %g P %g D %g => S %g\n",
               zscore, qlen, A, B, P, D, S);
#endif

        /* print parameters once */
        if (!print_blosum_pscore_parameters) {
                fprintf(stderr, "Blosum PSCORE parameters:\n");
                fprintf(stderr, " Q = %d\n", qlen);
                fprintf(stderr, " A = EVD_A_CONST %g + EVD_A_SLOPE %g * log(Q) = %g\n", EVD_A_CONST, EVD_A_SLOPE, A);
                fprintf(stderr, " B = EVD_B_CONST %g + EVD_B_SLOPE %g * log(Q) = %g\n", EVD_B_CONST, EVD_B_SLOPE, B);
                fprintf(stderr, " Z = %g\n", zscore);
                fprintf(stderr, " P = exp((A - Z)/B - exp((A - Z)/B))/B = %g\n", P);
                fprintf(stderr, " D = exp(-exp((A - Z)/B)) = %g\n", D);
                fprintf(stderr, " S = -log10(1.0 - D) = %g (blosum pscore)\n", S);
                print_blosum_pscore_parameters++;
        }
        return (S);
}

static double g_lambda = 0.267000;
static double g_kappa = 0.041000;

double natscore(double score)
{
        return (score * g_lambda - log(g_kappa));
}               
                        
double bitscore(double score)
{
        return (natscore(score) / log(2.0));
} 

int apos(char c)
/* return index of character c in scorematrix alphabet */
{
	char *l = strchr(alphabet, c);
	if (!l)
		fprintf(stderr, "Could not find char '%c' in alphabet '%s'\n", c, alphabet), exit(1);
	return l - alphabet;
}

#define UNDEFINED_SCOREMATRIX_VALUE 999.9
double scorematrix_element(char a, char b)
/* return scorematrix matrix element for characters a, b */
{
	/* ignore pad positions (outside of the local alignment) */
	if (a == ALIGN_PAD_CHAR || b == ALIGN_PAD_CHAR) {
		if (p_v)
			fprintf(stderr, "Element A %c PAD or B %c PAD %c\n", a, b, ALIGN_PAD_CHAR);
		return (UNDEFINED_SCOREMATRIX_VALUE);
	}
	if (a == ALIGN_GAP_CHAR || b == ALIGN_GAP_CHAR) {
		if (p_v)
			fprintf(stderr, "Element A %c GAP or B %c GAP %c\n", a, b, ALIGN_GAP_CHAR);
		return (UNDEFINED_SCOREMATRIX_VALUE);
	}
	int aa = apos(a), bb = apos(b);
	if (aa < 0 && bb < 0) {
		fprintf(stderr, "Element A %c pos %d, or B %c pos %d, not in '%s'\n", a, aa, b, bb, alphabet);
		return (UNDEFINED_SCOREMATRIX_VALUE);
	}
	return (blosum_mtx[aa][bb]);
}

int *residue_type_index_vector(char *fseq)
{
	int n = strlen(fseq), i;
	int *v = int_vector(n);
	for (i = 0; i < n; i++) {
		int a = apos(fseq[i]);
		v[i] = a;
	}
	return (v);
}

void print_scorematrix()
{
	int i, j;
	printf("# recovered scorematrix\n");
	printf(" ");
	for (j = 0; j < nb; j++)
		printf("  %c", alphabet[j]);
	printf("\n");
	for (i = 0; i < nb; i++) {
		printf("%c ", alphabet[i]);
		for (j = 0; j < nb; j++)
			printf("%2g ", blosum_mtx[i][j]);
		printf("\n");
	}
	/* Note: Command line 'diff' to compare the read and expected values of the matrix. */
}

void read_scorematrix(char *filename)
/* READ BLOSUM62 amino-acid substitution score matrix, and alphabet from a file
// in the very specific format provided with this source code.
// allocates and assigned global variables blosum_mtx, nb and alphabet */
{
	char line[MAXLINELEN], text[MAXLINELEN];
	int lineno = 0;
	FILE *fp;
	char word[10], c;
	int i, j, cnt = 0;
	double v;
	strcpy(alphabet, "");
	nb = 0;
	if ((fp = fopen(filename, "r")) == NULL)
		fprintf(stderr, "Blosum62 file %s not readable\n", filename), exit(1);
	fprintf(stderr, "Blosum62 file %s readable\n", filename);
	while (fgets(line, sizeof(line), fp) != NULL) {
		lineno++;
		if (sscanf(line, "%[^\n]", text) != 1)
			fprintf(stderr, "Cannot scan to newline\nlineno %2d:%s\n", lineno, line), exit(1);
		/* skip comments */
		if (text[0] == '#')
			continue;
		/* read matrix alphabet */
		if (text[0] == ' ') {
			for (i = 0; i < strlen(text); i++) {
				if (text[i] != ' ') {
					alphabet[nb++] = text[i];
					if (nb >= MAXN)
						fprintf(stderr, "Alphabet count %d exceeds MAXN %d\n", nb, MAXN), exit(1);
					alphabet[nb] = 0;
				}
			}
			i = 0;
			blosum_mtx = double_matrix(nb, nb);
			fprintf(stderr, "N %d alphabet %s\n", nb, alphabet);
			continue;
		}
		/* read matrix elements, */
		/* first character of each line should match next in alphabet */
		if (sscanf(text, "%c %[^\n]", &c, line) != 2)
			fprintf(stderr, "Cannot sscanf %%c:%s\n", text), exit(1);
		if (c != alphabet[i])
			fprintf(stderr, "Mismatch char '%c' != alphabet[%d] '%c'\n", c, i, alphabet[i]), exit(1);
		strcpy(text, line);
		j = 0;
		while (j < nb && sscanf(text, "%lf %[^\n]", &v, line) > 0) {
			if (p_v)
				fprintf(stderr, "i %2d '%c' j %2d '%c' value %f\n", i, alphabet[i], j, alphabet[j], v);
			strcpy(text, line);
			blosum_mtx[i][j] = v;
			j++;
			cnt++;
		}
		i++;
	}
	fprintf(stderr, "Read %d matrix elements, %d alphabet characters %s from %s\n",
		cnt, nb, alphabet, filename);
	fclose(fp);
	if (p_v)
		print_scorematrix();
}

/* SECTION FASTA */

/* global variables for input sources and fasta entries */
int g_nent = 0, g_nsrc = 0;
#define MAXENTRIES 10000
char *facc[MAXENTRIES];
char *fseq[MAXENTRIES];
int fsrc[MAXENTRIES];		/* source fasta file names, worst case one sequence per file */
/* int  *frti[MAXENTRIES]; residue index deactivated */

void print_fasta()
{
	int i;
	for (i = 0; i < g_nent; i++)
		printf(">%s\n%s\n", facc[i], fseq[i]);
}

/* used in parsing text */
char line[MAXLINELEN], text[MAXLINELEN], acc[MAXLINELEN], seq[MAXSEQUENCELEN];

void read_fasta(char *filename)
/* read simple sequence fasta file side effect: updates char *facc and *fseq and fsrc */
{
	int lineno = 0;
	FILE *fp;
	if ((fp = fopen(filename, "r")) == NULL)
		fprintf(stderr, "Fasta file %s not readable\n", filename), exit(1);
	fprintf(stderr, "Fasta file %s readable\n", filename);
	strcpy(acc, "");
	strcpy(seq, "");
	while (fgets(line, sizeof(line), fp) != NULL) {
		lineno++;
		if (sscanf(line, "%[^\n]", text) != 1) {
			fprintf(stderr, "Cannot sscanf %%[^\n]:%s\n", line);
			continue;
		}
		/* skip comments */
		if (text[0] == '#')
			continue;
		/* read accession (ignore rest of Fasta header) */
		if (text[0] == '>') {
			if (g_nent >= MAXENTRIES)
				fprintf(stderr, "Fasta count g_nent %d > MAXENTRIES %d\n", g_nent, MAXENTRIES), exit(1);
			if (strlen(acc) > 0) {
				facc[g_nent] = char_string(acc);
				fseq[g_nent] = char_string(seq);
				fsrc[g_nent] = g_nsrc;
				/* frti[g_nent] = residue_type_index_vector(fseq[g_nent]); */
				strcpy(acc, "");
				strcpy(seq, "");
				g_nent++;
			}
			if (sscanf(text, ">%s", acc) != 1)
				fprintf(stderr, "Cannot sscanf >%%s:%s\n", text), exit(1);
			continue;
		}
		/* otherwise accumulate sequence */
		strcat(seq, text);
	}
	if (g_nent >= MAXENTRIES)
		fprintf(stderr, "Count of fasta sequence records %d exceeds MAXENTRIES %d\n", g_nent, MAXENTRIES), exit(1);
	if (strlen(acc) > 0) {
		facc[g_nent] = char_string(acc);
		fseq[g_nent] = char_string(seq);
		fsrc[g_nent] = g_nsrc;
		/* frti[g_nent] = residue_type_index_vector(fseq[g_nent]); */
		strcpy(acc, "");
		strcpy(seq, "");
		g_nent++;
	}
	fprintf(stderr, "Read %d Fasta entries from %s\n", g_nent, filename);
	fclose(fp);
	if (p_v)
		print_fasta();
}

/* ALIGNMENT */

/* use some global variables to save time allocating/deallocating memory */
double **global_score_matrix = NULL, **global_match_matrix = NULL;
int global_seqlen = 0;

void pair_score_matrix(int f1, int f2)
/* provide Blosum62 substitution score matrix for pair of fasta elements fi, fj */
{
	char *s1 = fseq[f1], *s2 = fseq[f2];
	int i, j, n1 = strlen(s1), n2 = strlen(s2);
	if (n1 > global_seqlen || n2 > global_seqlen){
		if (global_score_matrix) {
			double_matrix_free(global_seqlen + 1, global_seqlen + 1, global_score_matrix);
			double_matrix_free(global_seqlen, global_seqlen, global_match_matrix);
		}
		global_seqlen = (n1 > n2 ? n1 : n2);
		global_score_matrix = double_matrix(global_seqlen + 1, global_seqlen + 1);
		global_match_matrix = double_matrix(global_seqlen, global_seqlen);
	}
	/* body of score matrix contains sequence I vs sequence J substitutions */
	for (i = 0; i < n1; i++)
		for (j = 0; j < n2; j++)
			global_score_matrix[i][j] = scorematrix_element(s1[i], s2[j]);
	/* edges of score matrix contain sequence I and sequence J self-scores */
	for (i = 0; i < n1; i++)
		global_score_matrix[i][n2] = scorematrix_element(s1[i], s1[i]);
	for (j = 0; j < n2; j++)
		global_score_matrix[n1][j] = scorematrix_element(s2[j], s2[j]);
}

double **global_T = NULL, **global_U = NULL, **global_V = NULL;
int global_N = 0;

double align_score(char *s1, char *s2, int n1, int n2, double **S, double **M, double fg, double ng, int *o1, int *o2, int align_flag)
/* Generate optimal local (Smit & Waterman 1981) alignment path using affine gap penalties
 * including an unpublished crossover gap allowance (Gippert 2001). Return alignment score
 * and (indirectly) sequence offsets. */
{
	int i, j;
	double ascore = 0.0;
	*o1 = -1;
	*o2 = -1;

	if (n1 > global_N || n2 > global_N) {
		if (global_T) {
			double_matrix_free(global_N + 1, global_N + 1, global_T);
			double_matrix_free(global_N + 1, global_N + 1, global_U);
			double_matrix_free(global_N + 1, global_N + 1, global_V);
		}
		global_N = (n1 > n2 ? n1 : n2);
		global_T = double_matrix(global_N + 1, global_N + 1);
		global_U = double_matrix(global_N + 1, global_N + 1);
		global_V = double_matrix(global_N + 1, global_N + 1);
	} else {
		for(i = 0; i <= global_N; i++)
			for (j = 0; j <= global_N; j++)
				global_T[i][j] = global_U[i][j] = global_V[i][j] = 0.0;
	}

	/* initialize transition matrices */
	double **T = global_T, **U = global_U, **V = global_V;

	/* additional tmp variables and row pointers */
	double *Ti, *Tp, *Si, *Vi, *Vp, *Ui, *Mi, *Up, t1, t2, t3, Tij, Uij, Vij;

	/* initialize far edge */
	for (j = n2 - 1; j >= 0; j--) {
		Tij = T[n1][j + 1];
		T[n1][j] = Tij;
		U[n1][j] = Tij;
		V[n1][j] = Tij - fg;
	}

	/* backwards from last row */
	for (i = n1 - 1; i >= 0; i--) {
		Si = S[i];
		Mi = M[i];
		Ti = T[i];
		Ui = U[i];
		Vi = V[i];
		Tp = T[i + 1];
		Up = U[i + 1];
		Vp = V[i + 1];
		Ti[n2] = Tp[n2];
		Ui[n2] = Ti[n2] - fg;
		Vi[n2] = Ti[n2];

		/* backwards from last column */
		for (j = n2 - 1; j >= 0; j--) {
			Tij = Tp[j + 1] + Si[j];
			Mi[j] = Tij;

			/* insertion in sequence 2 */
			t1 = Vp[j] - ng;
			t2 = Tp[j] - fg;
			Vij = (t1 > t2 ? t1 : t2);
			if (align_flag & ALIGN_CROSS) {
				t3 = Up[j] - ng;
				Vij = (t3 > Vij ? t3 : Vij);
			}
			Tij = (Vij > Tij ? Vij : Tij);

			/* insertion in sequence 1 */
			t1 = Ui[j + 1] - ng;
			t2 = Ti[j + 1] - fg;
			Uij = (t1 > t2 ? t1 : t2);
			if (align_flag & ALIGN_CROSS) {
				t3 = Vi[j + 1] - ng;
				Uij = (t3 > Uij ? t3 : Uij);
			}
			Tij = (Uij > Tij ? Uij : Tij);

			if (Tij > ascore) {
				ascore = Tij;
				if (p_v > 1)
					fprintf(stderr, "ASCORE i %d%c j %d%c ascore %g\n", i, s1[i], j, s2[j], ascore);
				(*o1) = i;
				(*o2) = j;
			}
			Tij = (Tij < 0.0 ? 0.0 : Tij);
			Ti[j] = Tij;
			Vi[j] = Vij;
			Ui[j] = Uij;
		}
	}
	return (ascore);
}

/* ALN structure */
typedef struct aln {
	struct aln *next;
	char *name1, *name2, *seq1, *seq2, *aln1, *aln2;	/* accession, sequence string and alignment string of each pair */
	int len1, len2, start1, start2, end1, end2;		/* sequence lengths and start/end coordinates (1-based) */
	int plen, alen, mlen, ilen, glen, olen, clen, nlen;	/* counts related to alignment */
	double gapcost, ascore, mscore, aprime, mprime, ab, mb, mscore1, mscore2, mscorer, sd0, sd1, sd2, sd;	/* properties of the alignment */
	double zscore, pscore;					/* related to CE alignments */
	double score, evalue, bitscore;				/* related to BLAST alignemnts */
} ALN;

ALN *aln_alloc()
/* return an allocated but empty ALN object */
{
	ALN *A = NULL;
	if ((A = (ALN *) malloc(sizeof(ALN))) == NULL)
		fprintf(stderr, "Could not allocate ALN object\n"), exit(1);
	A->next = NULL;
	A->name1 = A->name2 = A->seq1 = A->seq2 = A->aln1 = A->aln2 = NULL;
	A->len1 = A->len1 = A->start1 = A->start2 = A->end1 = A->end2 = -1;
	A->plen = A->alen = A->mlen = A->ilen = A->glen = A->olen = A->clen = A->nlen = 0;
	A->ascore = A->mscore = A->aprime = A->mprime = A->ab = A->mb = A->mscore1 = A->mscore2 = A->mscorer = A->mscore = A->sd0 = A->sd1 = A->sd2 = A->sd = 0.0;
	A->zscore = A->pscore = A->score = A->evalue = A->bitscore = -99.9;
	return (A);
}

void aln_free(ALN *A)
/* recursively call memory free of alignment object or list */
{
	ALN *N = A->next;
	if (A->name1) free(A->name1), A->name1 = NULL;
	if (A->name2) free(A->name2), A->name2 = NULL;
	if (A->seq1) free(A->seq1), A->seq1 = NULL;
	if (A->seq2) free(A->seq2), A->seq2 = NULL;
	if (A->aln1) free(A->aln1), A->aln1 = NULL;
	if (A->aln2) free(A->aln2), A->aln2 = NULL;
	free((char *)A);
	if (N)
		aln_free(N);
}

ALN *aln_obj(char *name1, char *name2, char *seq1, char *seq2, char *aln1, char *aln2)
/* return a populated ALN object */
{
	ALN *A = aln_alloc();
	if (name1) A->name1 = char_string(name1);
	if (name2) A->name2 = char_string(name2);
	if (seq1) A->seq1 = char_string(seq1);
	if (seq2) A->seq2 = char_string(seq2);
	if (aln1) A->aln1 = char_string(aln1);
	if (aln2) A->aln2 = char_string(aln2);
	/* do that with statistics here! */
	return(A);
}

void aln_write_json(ALN *A)
/* Write simple JSON, note, sets global variable jsnfp, which must be pre-initialized to NULL */
{
	if (! jsnfp) {
		char *filename = char_vector(strlen(oprefix) + strlen(".aln.js") + 1);
		sprintf(filename, "%s%s", oprefix, ".aln.js");
		fprintf(stderr, "jsnfile %s\n", filename);
		if ((jsnfp = fopen(filename, "w")) == NULL)
			fprintf(stderr, "Unable to open JSON file %s for writing\n", filename), exit(1);
		free(filename);
	}
	/* JSON format parsable line-by-line */
	j_opn(jsnfp);
	/* input */
	j_str(jsnfp, YES, "name1", A->name1, NULL, NULL);
	j_int(jsnfp, YES, "len1", A->len1, NULL, NULL);
	j_str(jsnfp, YES, "name2", A->name2, NULL, NULL);
	j_int(jsnfp, YES, "len2", A->len2, NULL, NULL);
	/* alignment strings */
	j_str(jsnfp, YES, "aln1", A->aln1, NULL, NULL);
	j_int(jsnfp, YES, "start1", A->start1, NULL, NULL);
	j_int(jsnfp, YES, "end1", A->end1, NULL, NULL);
	j_str(jsnfp, YES, "aln2", A->aln2, NULL, NULL);
	j_int(jsnfp, YES, "start2", A->start2, NULL, NULL);
	j_int(jsnfp, YES, "end2", A->end2, NULL, NULL);
	/* alignment counts and scores */
	j_int(jsnfp, YES, "plen", A->plen, NULL, NULL);
	j_int(jsnfp, YES, "alen", A->alen, NULL, NULL);
	j_int(jsnfp, YES, "mlen", A->mlen, NULL, NULL);
	j_int(jsnfp, YES, "ilen", A->ilen, NULL, NULL);
	j_int(jsnfp, YES, "glen", A->glen, NULL, NULL);
	j_int(jsnfp, YES, "olen", A->olen, NULL, NULL);
	j_int(jsnfp, YES, "clen", A->clen, NULL, NULL);
	j_int(jsnfp, YES, "nlen", A->nlen, NULL, NULL);
	j_dbl(jsnfp, YES, "ascore", A->ascore, NULL, NULL);
	j_dbl(jsnfp, YES, "mscore", A->mscore, NULL, NULL);
	j_dbl(jsnfp, YES, "mscore1", A->mscore1, NULL, NULL);
	j_dbl(jsnfp, YES, "mscore2", A->mscore2, NULL, NULL);
	j_dbl(jsnfp, YES, "mscorer", A->mscorer, NULL, NULL);
	/* score distances */
	j_dbl(jsnfp, YES, "sd0", A->sd0, NULL, NULL);
	j_dbl(jsnfp, YES, "sd1", A->sd1, NULL, NULL);
	j_dbl(jsnfp, YES, "sd2", A->sd2, NULL, NULL);
	j_dbl(jsnfp, NO, "sd", A->sd, NULL, NULL); /* last data element receives a NO to solve a json-related issue */
	j_cls(jsnfp);
	fflush(jsnfp);
}

void aln_write_stderr(ALN *A)
/* write text alignment to stderr */
{
	/* free-form text */
	fprintf(stderr, "Align %s %d x %s %d", A->name1, A->len1, A->name2, A->len2);
	fprintf(stderr, " O1 %d O2 %d", A->start1, A->start2);
	fprintf(stderr, " Plen %d Alen %d Mlen %d Ilen %d Glen %d Olen %d Clen %d Nlen %d",
		A->plen, A->alen, A->mlen, A->ilen, A->glen, A->olen, A->clen, A->nlen);
	fprintf(stderr, " Ascore %f Mscore %f", A->ascore, A->mscore);
	fprintf(stderr, " M1 %g M2 %g MR %g", A->mscore1, A->mscore2, A->mscorer);
	fprintf(stderr, " SD0 %g SD1 %g SD2 %g SD %g\n", A->sd0, A->sd1, A->sd2, A->sd);
	fprintf(stderr, "%s\n%s\n", A->aln1, A->aln2);
	fprintf(stderr, "\n");	/* extra newline for human readability */
}

void aln_write_text(ALN *A)
/* write plain structured text for alignment, requires that alnfp is initially NULL */
/* special case of supplied fp causes a write and exit */
{
	if (! alnfp) {
		char *filename = char_vector(strlen(oprefix) + strlen(".aln.txt") + 1);
		sprintf(filename, "%s%s", oprefix, ".aln.txt");
		if ((alnfp = fopen(filename, "w")) == NULL)
			fprintf(stderr, "Unable to open aln file %s for writing\n", filename), exit(1);
		free(filename);
	}
	/* free-form text */
	fprintf(alnfp, "Align %s %d x %s %d", A->name1, A->len1, A->name2, A->len2);
	fprintf(alnfp, " O1 %d O2 %d", A->start1, A->start2);
	fprintf(alnfp, " Plen %d Alen %d Mlen %d Ilen %d Glen %d Olen %d Clen %d Nlen %d",
		A->plen, A->alen, A->mlen, A->ilen, A->glen, A->olen, A->clen, A->nlen);
	fprintf(alnfp, " Ascore %f Mscore %f", A->ascore, A->mscore);
	fprintf(alnfp, " M1 %g M2 %g MR %g", A->mscore1, A->mscore2, A->mscorer);
	fprintf(alnfp, " SD0 %g SD1 %g SD2 %g SD %g\n", A->sd0, A->sd1, A->sd2, A->sd);
	fprintf(alnfp, "%s\n%s\n", A->aln1, A->aln2);
	fprintf(alnfp, "\n");	/* extra newline for human readability */
	fflush(alnfp);
}

void aln_write_binary(ALN *A)
/* stub for writing aln binary */
/* stub for writing aln binary, requires that blnfp is initially NULL */
{
	return;
}

#define MAXSCOREDIST 9999.9
double compute_scoredistance(double ma, double mr, double m1, double m2, double scale)
/* unpublished Gippert, G.P, ca 2009, sequence-length normalization of ScoreDist from Sonnhammer & Hollich 2005 */
{
	double num = ma - (mr * scale);
	double den = ((m1 + m2) / 2.0 - mr) * scale;
	double e = num / den, sd;
	if (e <= 0.0) {
		fprintf(stderr, "scoredist: e %g <= 0.0, set dist to MAXSCOREDIST %g\n", e, MAXSCOREDIST);
		sd = MAXSCOREDIST;
	}
	else {
		sd = -log(e) * 100;
	}
	return (sd);
}

/* Find alignment from cumultative matching scores */
ALN *align_ali(char *seq1, char *seq2, int len1, int len2, int o1, int o2, double **S, double **M, int align_flag)
/* return match score and indirectly *alen and alignment strings *a1 and *a2 */
{
	char *aln1, *aln2, *t1, *t2;
	int plen, alen, mlen, ilen, glen, olen, clen, nlen, i, j, k, l, nk, nl;
	double max, tmp, *scovec, ascore, gscore, mscore, mscore1, mscore2, mscorer, zs, ps;
	ALN *a = NULL;

	if (o1 != -1 && o2 != -1) {
		aln1 = aln2 = NULL;
		plen = alen = mlen = ilen = glen = olen = clen = nlen = 0;
		ascore = gscore = mscore = mscore1 = mscore2 = mscorer = ps = 0.0;
		zs = -99.9;
		t1     = char_vector(len1 + len2 + 2);
		t2     = char_vector(len1 + len2 + 2);

		if (align_flag & ALIGN_PAD) {
			for (i = 0; i < o1; i++) {
				t1[plen] = seq1[i];
				t2[plen] = ALIGN_PAD_CHAR;
				plen++;
			}
			for (j = 0; j < o2; j++) {
				t1[plen] = ALIGN_PAD_CHAR;
				t2[plen] = seq2[j];
				plen++;
			}
		}

		/* inject first aligned position into score, counts and alignment strings */
		k = o1;
		l = o2;
		mscore  += S[k][l];
		mscore1 += S[k][len2];
		mscore2 += S[len1][l];
		mscorer -= 1.0;
		t1[plen] = seq1[k];
		t2[plen] = seq2[l];
		max = M[k][l];
		plen++;
		alen++;
		mlen++;
		/* count identical, conserved, and non-negative matches */
		if (seq1[k] == seq2[l])
			ilen++;
		if (S[k][l] > 0.0)
			clen++;
		if (S[k][l] >= 0.0)
			nlen++;

		while ((k != len1 - 1) && (l != len2 - 1) && max > 0) {
			/* Find subsequent match */
			nk = k + 1;
			nl = l + 1;
			max = M[nk][nl];
			/* gap in column */
			for (i = k + 2; i < len1; i++) {
				tmp = M[i][l + 1] - p_go - 1.0 * (i - (k + 2)) * p_ge;
				if (tmp > max) {
					max = tmp;
					nk = i;
					nl = l + 1;
				}
			}
			/* gap in row */
			for (j = l + 2; j < len2; j++) {
				tmp = M[k + 1][j] - p_go - 1.0 * (j - (l + 2)) * p_ge;
				if (tmp > max) {
					max = tmp;
					nk = k + 1;
					nl = j;
				}
			}
			/* cross-over gap in both column and row ? */
#ifdef RAW_CODE
			/* this version is the most 'correct', although the crossover gap is not
			 * described in any other source code for SW local alignment, afaik. */
			if (align_flag & ALIGN_CROSS) {
				for (i = k + 1; i < len1; i++)
				for (j = l + 1; j < len2; j++) {
					int g = i - (k + 1) + j - (l + 1);
					tmp = M[i][j] - (g ? p_go + p_ge * (g - 1) : 0.0);
					if (tmp > max) {
						max = tmp;
						nk = i;
						nl = j;
					}
				}
			}
#else
			/* this version limits the search scope for a suitable crossover match to p_gx */
			if (p_gx && align_flag & ALIGN_CROSS) {
				int g, h, x = len1 - (k + 1) + len2 - (l + 1);
				x = (x > p_gx ? p_gx : x);
				for (g = 0; g < x; g++)
				for (h = 0; h <= g; h++) {
					i = k + 1 + h, j = l + 1 + g - h;
					if (i < len1 && j < len2) {
						tmp = M[i][j] - (g ? p_go + p_ge * (g - 1) : 0.0);
						if (tmp > max) {
							max = tmp;
							nk = i;
							nl = j;
						}
					}
				}
			}
#endif
			if (max > 0) {
				/* Insert gaps */
				for (i = k + 1; i < nk; i++) {
					t1[plen] = seq1[i];
					t2[plen] = '-';
					plen++;
					alen++;
				}
				for (j = l + 1; j < nl; j++) {
					t1[plen] = '-';
					t2[plen] = seq2[j];
					plen++;
					alen++;
				}
				/* and update gap/open count and total gapscore (penalty) */
				int g = nk - (k + 1) + nl - (l + 1);
				glen += g;
				olen += (g > 0 ? 1 : 0);
				gscore += (g > 0 ? p_go + (g - 1) * p_ge : 0.0);

				/* inject next match position into score, counts and alignment strings */
				k = nk;
				l = nl;
				mscore  += S[k][l];
				mscore1 += S[k][len2];
				mscore2 += S[len1][l];
				mscorer -= 1.0;
				t1[plen] = seq1[k];
				t2[plen] = seq2[l];
				max = M[k][l];
				plen++;
				alen++;
				mlen++;
				/* count identical, conserved, and non-negative matches */
				if (seq1[k] == seq2[l])
					ilen++;
				if (S[k][l] > 0.0)
					clen++;
				if (S[k][l] >= 0.0)
					nlen++;
			}
		}
		/* apply C-terminal padding. Position (nk, nl) refer to the index number of the last match */
		if (align_flag & ALIGN_PAD) {
			for (i = nk + 1; i < len1; i++) {
				t1[plen] = seq1[i];
				t2[plen] = ALIGN_PAD_CHAR;
				plen++;
			}
			for (j = nl + 1; j < len2; j++) {
				t1[plen] = ALIGN_PAD_CHAR;
				t2[plen] = seq2[j];
				plen++;
			}
		}

		/* terminate alignment strings */
		t1[plen] = 0;
		t2[plen] = 0;

		if (mlen > 0) {
			zs = blosum_zscore(mscore, mlen);
			/* the 'query' length is arbitrarily defined as the length
		   	of the first sequence */
			if (zs > 0.0)
				ps = blosum_pscore(zs, len1);
		}

		/* ALIGN contains copy of input sequences */
		if ((a = aln_alloc()) == NULL)
			fprintf(stderr, "cannot allocate ALN\n"), exit(1);
		a->len1 = len1;
		a->seq1 = char_string(seq1);
		a->len2 = len2;
		a->seq2 = char_string(seq2);
		a->start1 = o1 + 1, a->end1 = nk;	/* 1-based sequence start and end coordinates */
		a->start2 = o2 + 1, a->end2 = nl;	/* 1-based sequence start and end coordinates */
		a->aln1 = char_string(t1), free(t1);
		a->aln2 = char_string(t2), free(t2);
		a->plen = plen, a->alen = alen, a->mlen = mlen, a->glen = glen, a->ilen = ilen, a->olen = olen, a->clen = clen, a->nlen = nlen;
		a->gapcost = gscore, a->ascore = mscore - gscore, a->mscore = mscore, a->mscore1 = mscore1, a->mscore2 = mscore2, a->mscorer = mscorer, a->aprime = natscore(ascore), a->mprime = natscore(mscore), a->ab = bitscore(ascore), a->mb = bitscore(mscore), a->zscore = zs, a->pscore = ps;

		/* compute score distances */

		// normalized to shortest sequence length
		int minlen = (len1 < len2 ? len1 : len2);
		double scale =  (double)minlen/(double)mlen;
		a->sd = compute_scoredistance(mscore, mscorer, mscore1, mscore2, scale);

		// normalized to sequence1 length
		scale = (double)len1/(double)mlen;
		a->sd1 = compute_scoredistance(mscore, mscorer, mscore1, mscore2, scale);

		// normalize to sequence2 length
		scale = (double)len2/(double)mlen;
		a->sd2 = compute_scoredistance(mscore, mscorer, mscore1, mscore2, scale);

		// normalize to alignment length (original Sohnhammer Scoredist)
		scale = (double)alen/(double)mlen;
		a->sd0 = compute_scoredistance(mscore, mscorer, mscore1, mscore2, scale);
	}
	return(a);
}

int p_strict = 0;		/* set to 1, and no deviation is allowed in the recomputation of alignment score */

void align_stats(ALN *A, int expected_plen, double expected_ascore)
/* compute ALN statistics */
{
	if (A->aln1 == NULL) fprintf(stderr, "align_stats A->aln1 == NULL\n"), exit(1);
	if (A->aln2 == NULL) fprintf(stderr, "align_stats A->aln2 == NULL\n"), exit(1);

	// plen = count total (padded + gapped) length of alignment
	if ((A->plen = strlen(A->aln1)) != strlen(A->aln2))
		fprintf(stderr, "strlen(A->aln1) %ld != strlen(A->aln2) %ld\n%s\n%s\n", strlen(A->aln1), strlen(A->aln2), A->aln1, A->aln2), exit(1);

	// I forgot the use case for testing plen == expected_plen ...
	if (expected_plen > 0 && A->plen != expected_plen)
		fprintf(stderr, "A->plen %d != expected_plen %d\n", A->plen, expected_plen), exit(1);

	// alen = count aligned positions
	// mlen = count matched positions
	// ilen = count identical positions
	// glen = count gap positions
	// olen = count gap openings
	// clen = count conservative matches(score > 0)
	// nlen = count non - negative matches(score >= 0)
	// ascore = sum alignment match score minus total gap cost
	// mscore = sum alignment match score
	// mscore1 = sum sequence1 match score over aligned region
	// mscore2 = sum sequence2 match score over aligned region
	// mscorer = sum alignment random match score * /
	//o1, o2 = 1 - based sequence offsets to start of local alignment
	int i, isg = 0, o1 = -1, o2 = -1, n1 = 0, n2 = 0;

	/* compute gaps and match scores from the alignment strings */
	double sum_score = 0.0;
	for (i = 0; i < A->plen; i++) {
		char a = A->aln1[i];
		char b = A->aln2[i];
		/* track residue position n1 and n2 (1-based) */
		if (a != ALIGN_PAD_CHAR && a != ALIGN_GAP_CHAR)
			n1++;
		if (b != ALIGN_PAD_CHAR && b != ALIGN_GAP_CHAR)
			n2++;
		double s = scorematrix_element(a, b);
		if (p_v) {
			printf("p_go %g x olen %d + p_ge %g x (glen %d - olen %d)\n", p_go, A->olen, p_ge, A->glen, A->olen);
			if (fabs(s) < 100)
				sum_score += s;
			double gap_score = p_go * (float)(A->olen) + p_ge * (float)(A->glen - A->olen);
			printf("aln[ %d ] : '%d%c' vs '%d%c' :  score %g, sumscore %g, gapscore %g, total_score %g\n",
		       		i, n1, a, n2, b, s, sum_score, gap_score, (sum_score - gap_score));
		}
		/* ignore pad positions (outside of the local alignment) */
		if (a == ALIGN_PAD_CHAR || b == ALIGN_PAD_CHAR)
			continue;
		/* aligned positions */
		/* set offsets */
		if (o1 < 0)
			o1 = n1, o2 = n2;
		A->alen++;
		if (a == ALIGN_GAP_CHAR || b == ALIGN_GAP_CHAR) {
			A->glen++;
			if (!isg)
				A->olen++;
			isg = 1;
			continue;
		}
		/* matched positions */
		A->mlen++;
		if (a == b)
			A->ilen++;
		if (s > 0)
			A->clen++;
		if (s >= 0)
			A->nlen++;
		A->mscore += s;
		A->mscore1 += scorematrix_element(a, a);
		A->mscore2 += scorematrix_element(b, b);
		A->mscorer -= 1.0;	/* random match expectation */
		isg = 0;
	}
	if (A->seq1 && strlen(A->seq1) != n1) fprintf(stderr, "Sequence len(seq1) %ld != recomputed sequence length n1 %d\n", strlen(A->seq1), n1), exit(1);
	if (A->seq2 && strlen(A->seq2) != n2) fprintf(stderr, "Sequence len(seq2) %ld != recomputed sequence length n2 %d\n", strlen(A->seq2), n2), exit(1);
	A->start1 = o1;
	A->start2 = o2;
	A->len1 = n1;
	A->len2 = n2;

	A->gapcost = A->olen * p_go + (A->glen - A->olen) * p_ge;
	A->ascore = A->mscore - A->gapcost;

	/* compute score distances */
	/* original score distance of Sonnhammer & Hollich normalized to alignment length */
	double scale = (double)(A->alen) / (double)(A->mlen);
	A->sd0 = compute_scoredistance(A->mscore, A->mscorer, A->mscore1, A->mscore2, scale);

	/* normalize to first sequence length */
	scale = (double)n1 / (double)(A->mlen);
	A->sd1 = compute_scoredistance(A->mscore, A->mscorer, A->mscore1, A->mscore2, scale);

	/* normalize to second sequence length */
	scale = (double)n2 / (double)(A->mlen);
	A->sd2 = compute_scoredistance(A->mscore, A->mscorer, A->mscore1, A->mscore2, scale);

	/* preferred score distance to minimum length */
	int minlen = (n1 < n2 ? n1 : n2);
	scale = (double)minlen / (double)(A->mlen);
	A->sd = compute_scoredistance(A->mscore, A->mscorer, A->mscore1, A->mscore2, scale);

	/* TODO - identify the use case for this */
	if (expected_ascore > 0.0 && A->ascore != expected_ascore) {
		fprintf(stderr, "Pair %s %s computed A->ascore %g != expected %g\n", A->name1, A->name2, A->ascore, expected_ascore);
		aln_write_stderr(A);
		if (1 || p_strict)
			exit(1);
	}

	/* TODO: reintroduce the ALN data structure which contains the entire pairwise alignment
		then call to write it to JSON
		and call to write it to TALN (text alignment)
		and call to write it to BALN (binary alignment)

	   It might be better to decouple (again) alignment generation from tree generation.
		for example after writing the binary aligment to read it again as a way
		of determining the distance matrix, rather than populating the distance matrix
		one alignment at a time.

	   This would need to be tied into reading a distance matrix.

	   This would also need to be reconciled with reading a multiple alignment. Does it really make sense
	   to save the multiple alignment as an entire set of pairwise alignments? Spacewise not.

	   Eventually we could skip the text output, as it is a legacy format not part of any current workflow.
	   However, in case people have ever seen this code and use it, maybe the ascii part could be kept.

	   GarryG March 17, 2025
	*/
}

ALN *pair_malign(int fi, int fj)
/* return inferred alignment of fasta elements fi, fj */
{
	/* pre-aligned fasta sequences */
	char *a1 = fseq[fi], *a2 = fseq[fj];
	ALN *A = aln_obj(facc[fi], facc[fj], NULL, NULL, fseq[fi], fseq[fj]);
	align_stats(A, -1, -1.0);

	return(A);
}

ALN *pair_align(int fi, int fj)
/* return local alignment of  fasta elements fi, fj */
{
	int align_flag = ALIGN_GAP | ALIGN_PAD | ALIGN_CROSS;
	char *s1 = fseq[fi], *s2 = fseq[fj];
	int o1, o2, n1 = strlen(s1), n2 = strlen(s2);

	pair_score_matrix(fi, fj);	/* allocate global score and match matrix */
	double **sx = global_score_matrix, **mx = global_match_matrix;

	/* compute optimal alignment and alignment score */
	double ascore = align_score(s1, s2, n1, n2, sx, mx, p_go, p_ge, &o1, &o2, align_flag);
	ALN *A = align_ali(s1, s2, n1, n2, o1, o2, sx, mx, align_flag);
	A->name1 = char_string(facc[fi]);
	A->name2 = char_string(facc[fj]);

	return(A);
}

double **global_dmx = NULL;
double **global_imx = NULL;

double **align_fasta()
/* return all-vs-all pairwise distance matrix from multiple alignment or computed pairwise alignments
// return a symmetric matrix containing align scoredistance normalized
// to the shortest sequence length of the pair */
{
	int i, inext, j;
	global_dmx = double_matrix(g_nent, g_nent);	/* distance matrix */
	global_imx = double_matrix(g_nent, g_nent);	/* identity matrix */
	ALN *A;
	for (i = 0; i < g_nent; i++) {
		for (j = (p_nonself ? i + 1 : i); j < g_nent; j++) {
			A =  (p_maln ? pair_malign(i, j) : pair_align(i, j));

			/* write alignments */
			if (p_baln)
				aln_write_binary(A);
			if (p_jaln)
				aln_write_json(A);
			if (p_taln)
				aln_write_text(A);

			/* populate distance matrix with minlen score distance */
			global_dmx[i][j] = global_dmx[j][i] = A->sd;
			global_imx[i][j] = global_imx[j][i] = 100.0 * (double)(A->ilen)/(double)(A->mlen);

			aln_free(A);
		}
	}
	return (global_dmx);
}

/* EMBED */

int p_dim = 20;			/* embed dimension: default 20 */
char p_e = 'F';			/* 'D' = distance only, 'S' = single only, 'F' = full recursive embed */
int p_ilim = 1000;		/* embed iteration limit: default 1000 */
double p_clim = 1.0e-12;	/* embed polynomial convergence limit, default 1e-12 */

double eigvec(int n, double **mx, double *v, double *t)
/* Determine most significant eigenvector of matrix m by successive approximation
// (Crippen & Havel) */
{
	int i, j, count = 0;
	double norm;
	double_vector_drand48(n, t, -1.0, 1.0);
	double_vector_normal(n, t);
	double ratio = 100.0, value = 0.0, prev;

	do {
		prev = value;

		/* rotate trial vector by metric matrix */

		for (i = 0; i < n; i++) {
			v[i] = 0.0;
			for (j = 0; j < n; j++)
				v[i] += mx[i][j] * t[j];
		}

		/* compute dot product */

		value = 0.0;
		for (i = 0; i < n; i++)
			value += v[i] * t[i];

		/* check unity condition */

		if (value == 0.0)
			break;

		ratio = fabs((value - prev) / value);

		if (p_v)
			printf("# EIGVAL iter(%d) prev(%e) value(%e) conv(%e)\n", count, prev, value, ratio);

		/* normalize */

		norm = 0.0;
		for (i = 0; i < n; i++)
			norm += v[i] * v[i];
		norm = sqrt(norm);

		for (i = 0; i < n; i++)
			t[i] = v[i] / norm;

		count++;

	} while (count < p_ilim && ratio > p_clim);

	for (i = 0; i < n; i++)
		v[i] = t[i];

#ifdef DEBUG
	printf("# %12.4f eigenvalue found %3d / %3d iter., %g / %g conv.\n", value, count, p_ilim, ratio, p_clim);
#endif

	return (value);
}

void matrix_deflate(int n, double **mx, double *v, double e)
/* eliminate eigenspace of last eigenvalue from matrix, which reduces
 * its rank by one (from Crippen & Havel 1988) */
{
	int i, j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			mx[i][j] -= e * v[i] * v[j];
			mx[i][j] = NEARZERO(mx[i][j]);
		}
	}
}

double **metric_matrix(int n, double **d)
/* produce metric matrix from distance matrix d */
{
	double **mx = double_matrix(n, n), *comsqr = double_vector(n), radsqr = 0.0, dissqr;
	int i, j, nneg = 0;
	for (i = 0; i < n - 1; i++) {
		for (j = i + 1; j < n; j++) {
			dissqr = d[i][j] * d[i][j];
			radsqr += dissqr;
			comsqr[i] += dissqr;
			comsqr[j] += dissqr;
		}
	}
	radsqr /= (double)n;
	for (i = 0; i < n; i++) {
		comsqr[i] -= radsqr;
		comsqr[i] /= (double)n;
		if (comsqr[i] < 0.0)
			fprintf(stderr, "Info: metric_matrix, dim %d, squared center of mass %lf < 0.0, count %d\n",
				i, comsqr[i], ++nneg);
	}
	for (i = 0; i < n; i++)
		mx[i][i] = comsqr[i];
	for (i = 0; i < n - 1; i++) {
		for (j = i + 1; j < n; j++) {
			dissqr = d[i][j] * d[i][j];
			mx[j][i] = mx[i][j] = 0.5 * (comsqr[i] + comsqr[j] - dissqr);
		}
	}
	double_vector_free(n, comsqr);
	return (mx);
}

void double_vector_scale(int n, double *v, double scale)
/* multiple vector by scalar value */
{
	int i;
	for (i = 0; i < n; i++)
		v[i] *= scale;
}

double **embed_dmx(int n, double **d)
/* N-dimensional embedding of distance matrix using metric matrix distance geometry.
// (Crippen & Havel, 1988, p 311). Uses matrix exhaustion and matrix deflation. */
{
	int i, j;
	double **mx = metric_matrix(n, d);	/* metric matrix (dot products of COM vectors) */
	double **v = double_matrix(p_dim, n);	/* eigenvectors (dim x n) */
	double *e = double_vector(p_dim);	/* eigenvalues */
	double *t = double_vector(n);	/* tmp vector, reused many times in eigvec */
	double **coord;		/* coordinates (n x dim) */

	for (j = 0; j < p_dim; j++) {
		e[j] = eigvec(n, mx, v[j], t);
		matrix_deflate(n, mx, v[j], e[j]);
#ifdef DEBUG
		printf("Eigenvector j %d %10g %10g\n", j, e[j], -SIGN(e[j]) * sqrt(fabs(e[j])));
		for (i = 0; i < n; i++)
			printf(" i %d %10g %10g\n", i, v[j][i], v[j][i] * (-SIGN(e[j]) * sqrt(fabs(e[j]))));
#endif
		double_vector_scale(n, v[j], -SIGN(e[j]) * sqrt(fabs(e[j])));
	}
	double_matrix_free(n, n, mx);
	double_vector_free(p_dim, e);
	double_vector_free(n, t);
	/* return coordinates */
	coord = double_matrix(n, p_dim);
	for (i = 0; i < n; i++)
		for (j = 0; j < p_dim; j++)
			coord[i][j] = v[j][i];
	double_matrix_free(p_dim, n, v);
	return (coord);
}

/* TREE */

char p_r = 'N';			/* tree branch rotation: 'N' none, 'R' right, 'L' left */

#define MAXPOINTS MAXENTRIES
#define MAXNODES (2 * MAXPOINTS - 1)	/* total number of nodes required for binary tree */

/* BNODE is a node in a binary tree */
typedef struct bnode {
	struct bnode *left;	/* left subtree or NULL */
	struct bnode *right;	/* right subtree or NULL */
	struct bnode *parent;	/* parent node or NULL */
	double left_distance;	/* distance to left node or 0.0 */
	double right_distance;	/* distance to right node or 0.0 */
	double parent_distance;	/* distance to parent node or 0.0 */
	double *pos;		/* position in some embedding or NULL */
	int dim;		/* dimension in some embedding (length of pos) or 0 */
	int index;		/* index >= 0 of leaf node if defined, otherwise -1 */
} BNODE;

BNODE *bnode_alloc(BNODE * parent, int index)
{
	BNODE *B = NULL;
	if ((B = (BNODE *) malloc(sizeof(BNODE))) != NULL) {
		B->left = NULL;
		B->right = NULL;
		B->parent = NULL;
		B->left_distance = 0.0;
		B->right_distance = 0.0;
		B->parent_distance = 0.0;
		B->pos = NULL;
		B->dim = 0;
		B->index = -1;
	}
	if (parent)
		B->parent = parent;
	if (index >= 0)
		B->index = index;
	return (B);
}

/* recursive free BNODE tree */
void bnode_free(BNODE * B)
{
	if (B->left)
		bnode_free(B->left);
	if (B->right)
		bnode_free(B->right);
	if (B->pos)
		free((char *)B->pos);
	free((char *)B);
	B = NULL;
}

void bnode_vec_free(BNODE ** bnode, int n)
{
	int i;
	BNODE *B;
	for (i = 0; i < n; i++) {
		B = bnode[i];
		if (B->pos)
			free((char *)B->pos);
		free((char *)B);
	}
	free((char *)bnode);
	bnode = NULL;
}

BNODE **bnode_vec(int n, int alloc)
/* return vector of 'empty' Bnodes with no parent/left/right/pos/dim/index assignments */
{
	BNODE **vec;
	if ((vec = (BNODE **) malloc(sizeof(BNODE *) * n)) == NULL)
		fprintf(stderr, "Cannot allocate bnode_vec %d\n", n), exit(1);
	int i;
	if (alloc)
		for (i = 0; i < n; i++)
			vec[i] = bnode_alloc(NULL, -1);
	return vec;
}

void bnode_trace(BNODE * A)
{
	BNODE *p;
	int n = 0;
	printf("{trace\n");
	for (p = A; p; p = p->parent, n++) {
		printf("A treebranch %d parent %s %g (%g %g)\n",
		       n,
		       (NULL == p ? "unknowable" :
			(NULL == p->parent ? "null" :
			 (NULL == p->parent->left ? "leftnull" :
			  (p == p->parent->left ? "left" :
			   (NULL == p->parent->right ? "rightnull" :
			    (p == p->parent->right ? "right" :
			     "unknown")))))),
		       p->parent_distance,
		       (p->parent ? p->parent->left_distance : -999.9),
		       (p->parent ? p->parent->right_distance : -999.9));
	}
	printf("}\n\n");
}

/* treebranch dis between two nodes known to be in the same tree = sum of parent->left_distance */
double bnode_dis_treebranch(BNODE * A, BNODE * B)
{
	/* identify the common parent of these two nodes. */
	BNODE *pa, *pb, *C = NULL;
	for (pa = A; pa && !C; pa = pa->parent)
		for (pb = B; pb && !C; pb = pb->parent)
			if (pa == pb) {
				C = pa;
			}
	if (C == NULL)
		fprintf(stderr, "bnode_dis_treebranch, common ancestor not found!\n"), exit(1);

	/* distance from A to common ancestor C */
	int na = 0;
	double disa = 0.0;
	for (pa = A; pa && pa != C; pa = pa->parent, na++) {
		disa += pa->parent_distance;
		if (p_v > 1)
			printf("A %d parent %s %g (%g %g)\n",
			       na,
			       (NULL == pa ? "unknowable" :
				(NULL == pa->parent ? "null" :
				 (NULL == pa->parent->left ? "leftnull" :
				  (pa == pa->parent->left ? "left" :
				   (NULL == pa->parent->right ? "rightnull" :
				    (pa == pa->parent->right ? "right" :
				     "unknown")))))),
			       pa->parent_distance,
			       (pa->parent ? pa->parent->left_distance : -999.9),
			       (pa->parent ? pa->parent->right_distance : -999.9));
	}

	/* distance from B to common ancestor C */
	int nb = 0;
	double disb = 0.0;
	for (pb = B; pb && pb != C; pb = pb->parent, nb++) {
		disb += pb->parent_distance;
		if (p_v > 1)
			printf("B %d parent %s %g (%g %g)\n",
			       nb,
			       (NULL == pb ? "unknowable" :
				(NULL == pb->parent ? "null" :
				 (NULL == pb->parent->left ? "leftnull" :
				  (pb == pb->parent->left ? "left" :
				   (NULL == pb->parent->right ? "rightnull" :
				    (pb == pb->parent->right ? "right" :
				     "unknown")))))),
			       pb->parent_distance,
			       (pb->parent ? pb->parent->left_distance : -999.9),
			       (pb->parent ? pb->parent->right_distance : -999.9));
	}
	double dis = (disa + disb) / 2.0;
	if (p_v > 1)
		printf("Distance from A to (common ancestor) to B is %g\n", dis);
	return (dis);
}

double bnode_dis(BNODE * A, BNODE * B)
/* return Euclidean distance (norm 2) between two N-dimensional BNODEs */
{
	/* special case of positionless binary tree */
	if (A->pos == NULL && B->pos == NULL)
		return 0.0;
	if (A->pos == NULL)
		fprintf(stderr, "A->pos is NULL\n"), exit(1);
	if (B->pos == NULL)
		fprintf(stderr, "B->pos is NULL\n"), exit(1);
	if (A->dim < 1)
		fprintf(stderr, "A->dim < 1\n"), exit(1);
	if (B->dim < 1)
		fprintf(stderr, "B->dim < 1\n"), exit(1);
	if (A->dim != B->dim)
		fprintf(stderr, "A->dim %d != B->dim %d\n", A->dim, B->dim), exit(1);

	double cnt = 0.0, sum = 0.0, dis;
	int i;
	for (i = 0; i < A->dim; i++) {
		dis = A->pos[i] - B->pos[i];
		sum += dis * dis;
		cnt += 1.0;
	}
	dis = sqrt(sum);
	return dis;
}

int bnode_count(BNODE * B)
/* recursive number of defined leaf nodes in tree */
{
	if (B->index >= 0)
		return 1;
	else if (B->left == NULL || B->right == NULL)
		fprintf(stderr, "B->left is NULL or B->right is NULL\n"), exit(1);
	return bnode_count(B->left) + bnode_count(B->right);
}
int bnode_length(BNODE *B)
{
	return bnode_count(B);
}

int PRINTNL = 1;
int PRINTPOS = 0;
int PRINTDIS = 0;

/* assign B->pos to weighted average of left and right trees */
void bnode_average_pos(BNODE * B)
{

	if (B->index >= 0)
		fprintf(stderr, "bnode_average_pos: B->index %d >= 0\n", B->index), exit(1);
	if (!B->left)
		fprintf(stderr, "bnode_average_pos: B->left is NULL\n"), exit(1);
	if (!B->right)
		fprintf(stderr, "bnode_average_pos: B->right is NULL\n"), exit(1);
	if (B->pos)
		fprintf(stderr, "bnode_average_pos: B->pos is already defined\n"), exit(1);
	int dim = B->left->dim;
	if (dim == 0)
		fprintf(stderr, "bnode_average_pos: B->left->dim is 0\n"), exit(1);
	if (B->right->dim != dim)
		fprintf(stderr, "bnode_average_pos: B->right->dim %d != B->left->dim %d\n", B->right->dim, B->left->dim), exit(1);

	/* weighted average of left and right position */
	double lw = (double)bnode_count(B->left);
	double rw = (double)bnode_count(B->right);
	double *pos;
	pos = double_vector(dim);
	int k;
	for (k = 0; k < dim; k++)
		pos[k] = (lw * B->left->pos[k] + rw * B->right->pos[k]) / (lw + rw);
	B->pos = pos;
	B->dim = dim;
}

void printpos(FILE * fp, BNODE * B)
{
	if (!B->pos || !B->dim)
		return;
	int k;
	fprintf(fp, "[%g", B->pos[0]);
	for (k = 1; k < B->dim; k++)
		fprintf(fp, ",%g", B->pos[k]);
	fprintf(fp, "]");
	if (PRINTNL)
		fprintf(fp, "\n");
}

/* recursive assign ith value of index vector to tree-ordered leaf index */
void bnode_indexi(BNODE * B, int *index, int *i)
{
	if (B->index >= 0)
		index[(*i)++] = B->index;
	else if (B->left && B->right) {
		bnode_indexi(B->left, index, &(*i));
		bnode_indexi(B->right, index, &(*i));
	}
	else
		fprintf(stderr, "bnode_indexi: B->left xor B->right in bnode_indexi\n"), exit(1);
}

void avesd(double **mx, int n, int *index, int m, int *jndex, double *ave, double *sd)
/* return average and population standard deviation of values at (n, index) vs (m, jndex) matrix elements */
{
	*ave = *sd = 0.0;
	double val, sum = 0.0, cnt = 0.0, sqr = 0.0;
	int i, j;
	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++) {
			val = mx[index[i]][jndex[j]];
			sum += val;
			cnt += 1.0;
			sqr += val * val;
		}
	if (cnt > 0.0){
		*ave = sum/cnt;
		*sd = sqrt(fabs(sqr/cnt - (*ave)*(*ave)));
	}
}

void bnode_print_metadata(FILE *fp, BNODE *left, BNODE *right)
/* print metadata comparing left and right branches */
{

	int n = bnode_count(left);	/* number of leaves in left subtree */
	int *index = int_vector(n), i = 0;
	bnode_indexi(left, index, &i);

	int m = bnode_count(right);	/* number of leaves in right subtree */
	int *jndex = int_vector(m), j = 0;
	bnode_indexi(right, jndex, &j);

	double pct, sd_pct;
	avesd(global_imx, n, index, m, jndex, &pct, &sd_pct);

	double dis, sd_dis;
	avesd(global_dmx, n, index, m, jndex, &dis, &sd_dis);

	fprintf(fp, "[&pct=%.1f,dis=%.1f,sd_pct=%.1f,sd_dis=%.1f]", pct, dis, sd_pct, sd_dis);

	free((char *)index);
	free((char *)jndex);
}

/* print binary tree */
void bnode_print_tree(FILE * fp, BNODE * B)
{

	if (p_v > 1)
		printf("{parent %s %g (%g %g)}",
		       (NULL == B ? "undef" :
			(NULL == B->parent ? "null" :
			 (B == B->parent->left ? "left" :
			  (B == B->parent->right ? "right" : "unknown")))),
		       B->parent_distance,
		       (B->parent ? B->parent->left_distance : -999.9),
		       (B->parent ? B->parent->right_distance : -999.9));
	if (B->index >= 0)
		fprintf(fp, "%s", facc[B->index]);
	else if (B->left && B->right) {
		/* left */
		fprintf(fp, "(");
		if (PRINTPOS)
			printpos(fp, B->left);
		bnode_print_tree(fp, B->left);
		fprintf(fp, ":%g", B->left_distance);
		if (PRINTNL)
			fprintf(fp, "\n");
		/* right */
		fprintf(fp, ",");
		if (PRINTPOS)
			printpos(fp, B->right);
		bnode_print_tree(fp, B->right);
		fprintf(fp, ":%g", B->right_distance);
		if (PRINTNL)
			fprintf(fp, "\n");
		fprintf(fp, ")");
		if (p_metadata)
			bnode_print_metadata(fp, B->left, B->right);
	}
}

void bnode_print(FILE * fp, BNODE * B)
{
	bnode_print_tree(fp, B);
	fprintf(fp, ":0;");	/* terminate tree */
	if (PRINTNL)
		fprintf(fp, "\n");
}

void bnode_bnodei(BNODE * B, BNODE ** bnode, int *i)
/* recursive assign ith value of bnode vector to tree-ordered leaf node
// DO NOT DELETE the indirectly-allocated bnode in the usual way as it
// might also delete parts of the tree it points to */
{
	if (B->index >= 0)
		bnode[(*i)++] = B;
	else if (B->left && B->right) {
		bnode_bnodei(B->left, bnode, &(*i));
		bnode_bnodei(B->right, bnode, &(*i));
	}
	else
		fprintf(stderr, "bnode_bnodei: B->left xor B->right in bnode_bnodei\n"), exit(1);
}

#define DMX_ONE		0x00000001  /* == 1 */
#define DMX_TWO		0x00000002  /* == 2 */
#define DMX_THREE	0x00000004  /* == 4 */
#define DMX_FOUR	0x00000008  /* == 8 */
#define DMX_FIVE	0x00000010  /* == 16 */
#define DMX_SIX		0x00000020  /* == 32 */
#define DMX_SEVEN	0x00000040  /* == 64 */


/* provide a binary tree from a DISTANCE matrix using nearest-neighbor joining algorithm and DISTANCE averaging (yuck) */
/* dmx_flag should control distance calculation at nnj/node creation */
BNODE *bnode_tree_dmx(int n, int *index, double **dmx, int dmx_flag)
{

	if (dmx_flag & DMX_ONE) { fprintf(stderr, "DMX ONE is coded\n"); }
	else
	if (dmx_flag & DMX_TWO) { printf("DMX TWO\n"), exit(0); }
	else
	if (dmx_flag & DMX_THREE) { printf("DMX THREE\n"), exit(0); }
	else
	if (dmx_flag & DMX_FOUR) { printf("DMX FOUR\n"), exit(0); }
	else
	if (dmx_flag & DMX_FIVE) { printf("DMX FIVE\n"), exit(0); }
	else
	if (dmx_flag & DMX_SIX) { printf("DMX SIX\n"), exit(0); }
	else
	if (dmx_flag & DMX_SEVEN) { printf("DMX SEVEN\n"), exit(0); }
	else { printf("DMX NONE of the above\n"), exit(1); }

	/* algorithm: // allocate BNODE vector of length 2N-1, with N of them with 'natural' index numbers (leaf nodes)
	   // and N-1 of them with index number -1 (tree nodes) without assigning any of them parentage. // allocate
	   integer 'use' vector of length 2N-1, // alternatively double 'weight' vector of length 2N-1, 1.0 for leafs
	   0.0 for others... // allocate dmatrix of length 2N-1 x 2N-1, with first N positions taken by actual
	   distances // find the minimum non-used element in the matrix, join the nodes. */

	int nodes = 2 * n - 1;
	if (nodes > MAXNODES)
		fprintf(stderr, "bnode_tree: nodes %d > MAXNODES %d\n", nodes, MAXNODES), exit(1);

	BNODE **bvec = bnode_vec(nodes, 1);

	/* first N BNODE elements are leaf nodes */
	int i;
	for (i = 0; i < n; i++) {
		bvec[i]->index = (index ? index[i] : i);
		//bvec[i]->dim = dim;
		//bvec[i]->pos = double_vector_copy(dim, pos[i]);
	}

	/* extended copy of dmx, availability vector */
	double **smx = double_matrix(nodes, nodes);
	int *avail = int_vector(nodes);

	int j, k;
	for (i = 0; i < n; i++) {
		avail[i] = 1;
		for (j = i + 1; j < n; j++) {
			smx[j][i] = smx[i][j] = dmx[i][j];
		}
	}

	if (p_v > 1) {
		fprintf(stderr, "SMX n %d nodes %d\n", n, nodes);
		for (i = 0; i < n; i++) {
			fprintf(stderr, "%d:", i);
			for (j = 0; j < n; j++)
				fprintf(stderr, " %7.4g", smx[i][j]);
			fprintf(stderr, "\n");
		}
	}

	/* parent will be one of the allocated nodes */
	BNODE *P;
	P = NULL;

	/* nearest neighbor-joining algorithm */
	int m = n;		/* next available interior node */
	while (1) {
		double mins = 1e36;
		int mini = MAXNODES;
		int minj = MAXNODES;
		int found = 0;
		/* find available (i,j) pair with lowest distance */
		for (i = 0; i < m; i++) {
			if (avail[i]) {
				for (j = i + 1; j < m; j++) {
					if (avail[j]) {
						if (smx[i][j] < mins) {
							mins = smx[i][j];
							mini = i;
							minj = j;
							found++;
						}
					}
				}
			}
		}
		if (!found)
			break;
		if (p_v > 1)
			fprintf(stderr, "found new pair %d %d\n", mini, minj);

		/* point to next available node */
		P = bvec[m];

		/* left and right subtrees */
		BNODE *A = bvec[mini];
		BNODE *B = bvec[minj];

		/* optionally 'rotate' left and right */
		if (p_r == 'N') {	/* no rotation */
			P->left = A;
			P->right = B;
		} else {
			int ac = bnode_count(A);
			int bc = bnode_count(B);
			if (p_r == 'L') {	/* largest subtree to left */
				P->left = (ac >= bc ? A : B);
				P->right = (ac >= bc ? B : A);
			}
			else if (p_r == 'R') {	/* largest subtree to right */
				P->left = (ac >= bc ? B : A);
				P->right = (ac >= bc ? A : B);
			}
		}

		/* and reconnect to parent */
		P->left->parent = P;
		P->right->parent = P;

		/* assign left, right and subtree->parent distances */
		double L = (float)bnode_count(P->left);
		double R = (float)bnode_count(P->right);
		P->left_distance = P->left->parent_distance = smx[mini][minj] * R / (L + R);
		P->right_distance = P->right->parent_distance = smx[mini][minj] * L / (L + R);

		if (p_v > 1) {
			fprintf(stderr, "intermediate tree\n");
			bnode_print(stderr, P);
		}

		if (p_v > 1)
			bnode_trace(P->left);

		if (p_v > 1)
			bnode_trace(P->right);

		/* update distance matrix */
		avail[mini] = 0;
		avail[minj] = 0;

		if (dmx_flag & DMX_ONE) {
			/* Average distance model */
			double del;
			for (i = 0; i < m; i++) {
				if (avail[i]) {
					del = (smx[i][mini] + smx[i][minj]) / 2.0;
					smx[i][m] = smx[m][i] = del;
				}
			}
		}
		else
		if (dmx_flag & DMX_TWO) {
			/* Average leaf-to-leaf across branches distance model */
			fprintf(stderr, "DMX_TWO not implemented yet\n"); exit(1);
			double del;
			for (i = 0; i < m; i++) {
				if (avail[i]) {
					del = (smx[i][mini] + smx[i][minj]) / 2.0;
					smx[i][m] = smx[m][i] = del;
				}
			}
		}
		else {
			printf("dmx_flag value %d apparently not coded\n", dmx_flag);
			fprintf(stderr, "dmx_flag value %d apparently not coded\n", dmx_flag);
			exit(1);
		}

		/* register availability and increment to next available node */
		avail[m] = 1;
		m++;
	}

	if (m != nodes) {
		for (i = 0; i < nodes; i++)
			if (avail[i])
				fprintf(stderr, "node unused:  i %d, left %s, right %s, parent %s, pos %s, index %d, label %s\n",
					i,
					(bvec[i]->left ? "def" : "undef"),
					(bvec[i]->right ? "def" : "undef"),
					(bvec[i]->parent ? "def" : "undef"),
					(bvec[i]->pos ? "def" : "undef"),
					bvec[i]->index, facc[bvec[i]->index]);
		fprintf(stderr, "Fatal: m != nodes : n %d, m %d, nodes %d\n", n, m, nodes), exit(1);
	}
	if (!P)
		fprintf(stderr, "bnode_tree: tree is NULL\n"), exit(1);

	/* free */
	double_matrix_free(nodes, nodes, smx);
	int_vector_free(nodes, avail);
	/* NO bnode_vec_free(bvec, nodes); */

	return P;
}

/* provide a binary tree from a position matrix with dimensions n x dim, using nearest-neighbor joining algorithm */
BNODE *bnode_tree(double **pos, int *index, int n, int dim, double **dmx)
{

	/* algorithm: // allocate BNODE vector of length 2N-1, with N of them with 'natural' index numbers (leaf nodes)
	   // and N-1 of them with index number -1 (tree nodes) without assigning any of them parentage. // allocate
	   integer 'use' vector of length 2N-1, // alternatively double 'weight' vector of length 2N-1, 1.0 for leafs
	   0.0 for others... // allocate dmatrix of length 2N-1 x 2N-1, with first N positions taken by actual
	   distances // find the minimum non-used element in the matrix, join the nodes. */

	int nodes = 2 * n - 1;
	if (nodes > MAXNODES)
		fprintf(stderr, "bnode_tree: nodes %d > MAXNODES %d\n", nodes, MAXNODES), exit(1);

	BNODE **bvec;
	bvec = bnode_vec(nodes, 1);	/* allocate bnode vector */

	/* first N BNODE elements are leaf nodes */
	int i;
	for (i = 0; i < n; i++) {
		bvec[i]->index = (index ? index[i] : i);
		bvec[i]->dim = dim;
		bvec[i]->pos = double_vector_copy(dim, pos[i]);
	}

	/* distsqr matrix and availability vector initially fractionally occupied */
	double **smx;
	smx = double_matrix(nodes, nodes);
	int *avail;
	avail = int_vector(nodes);

	int j, k;
	for (i = 0; i < n; i++) {
		avail[i] = 1;
		for (j = i + 1; j < n; j++) {
			smx[i][j] = 0.0;
			for (k = 0; k < dim; k++) {
				double del = pos[i][k] - pos[j][k];
				smx[i][j] += del * del;
			}
		}
	}

	if (p_v > 1) {
		fprintf(stderr, "SMX n %d nodes %d dim %d\n", n, nodes, dim);
		for (i = 0; i < n; i++) {
			fprintf(stderr, "%d:", i);
			for (j = 0; j < n; j++)
				fprintf(stderr, " %7.4g", smx[i][j]);
			fprintf(stderr, "\n");
		}
	}

	/* parent will be one of the allocated nodes */
	BNODE *P;
	P = NULL;

	/* nearest neighbor-joining algorithm */
	int m = n;		/* next available interior node */
	while (1) {
		double mins = 1e36;
		int mini = MAXNODES;
		int minj = MAXNODES;
		int found = 0;
		/* find available (i,j) pair with lowest distance */
		for (i = 0; i < m; i++) {
			if (avail[i]) {
				for (j = i + 1; j < m; j++) {
					if (avail[j]) {
						if (smx[i][j] < mins) {
							mins = smx[i][j];
							mini = i;
							minj = j;
							found++;
						}
					}
				}
			}
		}
		if (!found)
			break;
		if (p_v > 1)
			fprintf(stderr, "found new pair %d %d\n", mini, minj);

		/* point to next available node */
		P = bvec[m];

		/* left and right subtrees */
		BNODE *A = bvec[mini];
		BNODE *B = bvec[minj];
		int ac = bnode_count(A);
		int bc = bnode_count(B);

		/* optionally 'rotate' left and right */
		if (p_r == 'L') {	/* largest subtree to left */
			P->left = (ac >= bc ? A : B);
			P->right = (ac >= bc ? B : A);
		}
		else if (p_r == 'R') {	/* largest subtree to right */
			P->left = (ac >= bc ? B : A);
			P->right = (ac >= bc ? A : B);
		}
		else {		/* no rotation */
			P->left = A;
			P->right = B;
		}

		/* and reconnect to parent */
		P->left->parent = P;
		P->right->parent = P;

		/* compute P position from weighted average of left and right subtrees */
		bnode_average_pos(P);

		if (p_v > 1) {
			fprintf(stderr, "average pos of new node: [%g", P->pos[0]);
			for (k = 1; k < dim; k++)
				fprintf(stderr, ",%g", P->pos[k]);
			fprintf(stderr, "]\n");
		}

		/* assign left, right and subtree->parent distances */
		P->left_distance = P->left->parent_distance = bnode_dis(P, P->left);
		P->right_distance = P->right->parent_distance = bnode_dis(P, P->right);

		if (p_v > 1) {
			fprintf(stderr, "intermediate tree\n");
			bnode_print(stderr, P);
		}

		if (p_v > 1)
			bnode_trace(P->left);

		if (p_v > 1)
			bnode_trace(P->right);

		/* update sqr distance matrix with new distance information */
		double del;
		for (i = 0; i < m; i++) {
			if (avail[i]) {
				smx[i][m] = 0.0;
				for (k = 0; k < dim; k++) {
					del = bvec[i]->pos[k] - P->pos[k];
					smx[i][m] += del * del;
				}
			}
		}

		/* register availability and increment to next available node */
		avail[mini] = 0;
		avail[minj] = 0;
		avail[m] = 1;
		m++;
	}

	if (m != nodes) {
		for (i = 0; i < nodes; i++)
			if (avail[i])
				fprintf(stderr, "node unused:  i %d, left %s, right %s, parent %s, pos %s, index %d, label %s\n",
					i,
					(bvec[i]->left ? "def" : "undef"),
					(bvec[i]->right ? "def" : "undef"),
					(bvec[i]->parent ? "def" : "undef"),
					(bvec[i]->pos ? "def" : "undef"),
					bvec[i]->index, facc[bvec[i]->index]);
		fprintf(stderr, "Fatal: m != nodes : n %d, m %d, nodes %d\n", n, m, nodes), exit(1);
	}
	if (!P)
		fprintf(stderr, "bnode_tree: tree is NULL\n"), exit(1);

	/* free */
	double_matrix_free(nodes, nodes, smx);
	int_vector_free(nodes, avail);
	/* NO bnode_vec_free(bvec, nodes); */

	return P;
}

void printdmx(FILE * fp, double **dmx, int *index, char **label, int n)
{
	fprintf(fp, "SUB DMX %d\n", n);
	int i, j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++)
			fprintf(fp, " %5s", "");
		for (j = i; j < n; j++)
			fprintf(fp, " %5.2f", dmx[i][j]);
		fprintf(fp, " index %d label %s\n", index[i], label[i]);
	}
}

void bnode_positional_distance_difference(BNODE * B, double **odmx)
{
	int n = bnode_count(B);
	/* number of leaves in this subtree */

	/* assign indices to local vector */
	int *index = int_vector(n), i = 0, j;
	bnode_indexi(B, index, &i);
	if (i != n)
		fprintf(stderr, "bnode_positional_distance_difference: after bnode_indexi i %d != n %d\n", i, n), exit(1);

	/* assign bnode lookup vector */
	BNODE **bnode;
	bnode = bnode_vec(n, 0);/* allocate bnode pointer vector */
	i = 0;
	bnode_bnodei(B, bnode, &i);
	if (i != n)
		fprintf(stderr, "bnode_positional_distance_difference: after bnode_bnodei i %d != n %d\n", i, n), exit(1);
	double sum = 0.0, cnt = 0.0, sqr = 0.0, ave, rms, ddmx, dpos, ddif;
	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			ddmx = odmx[index[i]][index[j]];
			dpos = bnode_dis(bnode[i], bnode[j]);
			ddif = dpos - ddmx;
			if (p_v > 1)
				printf(" ij %d %d dpos %g ddmx %g ddif %g\n",
				       i, j, dpos, ddmx, ddif);
			sum += ddif;
			sqr += ddif * ddif;
			cnt += 1.0;
		}
	}
	if (cnt > 0.0) {
		ave = sum / cnt;
		sqr = sqr / cnt;
		rms = sqrt(sqr - ave * ave);
	}
	if (p_v > 1)
		fprintf(stderr, "bnode_positional_distance_difference n= %d distance (positional - matrix) ave= %g rms= %g\n",
			n, ave, rms);

	free((char *)index);
	free((char *)bnode);
}

void bnode_treebranch_distance_difference(BNODE * B, double **odmx)
{
	/* number of leaves in this subtree */
	int n = bnode_count(B);

	/* assign indices to local vector */
	int *index = int_vector(n), i = 0, j;
	bnode_indexi(B, index, &i);
	if (i != n)
		fprintf(stderr, "bnode_treebranch_distance_difference: after bnode_indexi i %d != n %d\n", i, n), exit(1);

	/* assign bnode lookup vector */
	BNODE **bnode;
	bnode = bnode_vec(n, 0);/* allocate bnode pointer vector */
	i = 0;
	bnode_bnodei(B, bnode, &i);
	if (i != n)
		fprintf(stderr, "bnode_treebranch_distance_difference: after bnode_bnodei i %d != n %d\n", i, n), exit(1);

	double sum = 0.0, cnt = 0.0, sqr = 0.0, ave, rms, ddmx, dbln, ddif;
	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			ddmx = odmx[index[i]][index[j]];
			dbln = bnode_dis_treebranch(bnode[i], bnode[j]);
			ddif = dbln - ddmx;
			if (p_v > 1)
				printf(" ij %d %d dbln %g ddmx %g ddif %g\n",
				       i, j, dbln, ddmx, ddif);
			sum += ddif;
			sqr += ddif * ddif;
			cnt += 1.0;
		}
	}
	if (cnt > 0.0) {
		ave = sum / cnt;
		sqr = sqr / cnt;
		rms = sqrt(sqr - ave * ave);
	}
	free((char *)index);
	free((char *)bnode);
}

/* Re-embed and return new subtree using original distance matrix */
BNODE *bnode_reembed(BNODE * B, char br, double **odmx, int on, int dim)
{

	if (B->index >= 0)
		return B;

	bnode_positional_distance_difference(B, odmx);
	bnode_treebranch_distance_difference(B, odmx);

	int n = bnode_count(B);	/* number of leaves in this subtree */
	if (p_v > 1)
		fprintf(stderr, "bnode_reembed br=%c n= %d on= %d dim %d\n", br, n, on, dim);

	if (p_v > 1)
		bnode_print(stderr, B);

	/* assign indices to local vector */
	int *index = int_vector(n), i = 0, j;
	bnode_indexi(B, index, &i);
	if (i != n)
		fprintf(stderr, "bnode_reembed, i %d != n %d after bnode_indexi\n", i, n), exit(1);

	/* assign local distance matrix */
	double **dmx = double_matrix(n, n);
	for (i = 0; i < n; i++) {
		dmx[i][i] = 0.0;
		for (j = i + 1; j < n; j++)
			dmx[i][j] = dmx[j][i] = odmx[index[i]][index[j]];
	}

	/* embed the new (sub)matrix */
	double **pos = embed_dmx(n, dmx);
	if (!pos)
		fprintf(stderr, "bnode_reembed, pos is NULL\n"), exit(1);

	if (p_v > 1) {
		fprintf(stderr, "POS\n");
		for (i = 0; i < n; i++) {
			fprintf(stderr, "pos :");
			for (j = 0; j < dim; j++)
				fprintf(stderr, " %g", pos[i][j]);
			fprintf(stderr, "\n");
		}
	}

	/* generate a subtree from the new positions */
	BNODE *P;
	P = bnode_tree(pos, index, n, dim, odmx);

	if (p_v > 1) {
		fprintf(stderr, "after reembedn\n");
		bnode_print(stderr, P);
	}

	/* transfer information */
	P->parent = B->parent;
	P->parent_distance = B->parent_distance;

	/* free old subtree */
	bnode_free(B);

	/* connect */
	if (P->parent) {
		if (br == 'L')
			P->parent->left = P;
		else if (br == 'R')
			P->parent->right = P;
		else
			fprintf(stderr, "bnode_reembed, unknown branch rotation br '%c'\n", br), exit(1);
	}
	/* note - do not change P->parent->{left,right}_distance // as they are determined by a 'higher' level
	   embedding // recursively apply this to left and right subbranches */

	P->left = bnode_reembed(P->left, 'L', odmx, on, dim);
	P->right = bnode_reembed(P->right, 'R', odmx, on, dim);

	double_matrix_free(n, n, dmx);
	double_matrix_free(n, dim, pos);
	int_vector_free(n, index);

	/* return new subtree */
	return P;

}

void write_tree(BNODE * P, char *filename)
{
	FILE *fp;
	if ((fp = fopen(filename, "w")) == NULL)
		fprintf(stderr, "Tree file %s cannot be opened for writing\n", filename), exit(1);
	bnode_print(fp, P);
	int n = bnode_count(P);
	fprintf(stderr, "Tree file %s with %d nodes written\n", filename, n);
	fclose(fp);
}

void write_tree_binary(BNODE * P, char *filename)
{
	fprintf(stderr, "write_tree_binary stub bnode %d, filename %s\n", bnode_length(P), filename);
	return ;

	printf("You are attempting to write_tree_binary to %s, unsuccessfully! treelength %d\n", filename, bnode_count(P));
	fprintf(stderr, "You are attempting to write_tree_binary %s, unsuccessfully! treelength %d\n", filename, bnode_count(P));
	exit(1);
}

void write_dmx(double **dmx, char *oprefix)
/* print distance upper half matrix plus diagonal */
{
	char *filename = char_vector(strlen(oprefix) + strlen(".dmx.txt") + 1);
	sprintf(filename, "%s%s", oprefix, ".dmx.txt");
	FILE *fp;
	int i, j, c;
	if ((fp = fopen(filename, "w")) == NULL)
		fprintf(stderr, "Distance file %s cannot be opened for writing\n", filename), exit(1);
	for (c = 0, i = 0; i < g_nent; i++)
		for (j = i; j < g_nent; j++, c++)
			fprintf(fp, "%s %s %f\n", facc[i], facc[j], dmx[i][j]);
	fclose(fp);
	fprintf(stderr, "Wrote dmx %s %d elements, expect N(N+1)/2 %d for N=%d\n",
		filename, c, g_nent*(g_nent+1)/2, g_nent);
	free(filename);
}

BNODE *bnode_distance_tree(int n, double **dmx)
{
	int *index = int_vector_ramp(n);
	int dmx_flag = DMX_ONE;
	BNODE *P = bnode_tree_dmx(n, index, dmx, dmx_flag);
	return(P);
}
/* After
// 1. alignments have now provided a distance matrix
// The the task of tree building can begin
// 2. embed distance matrix into orthogonal coordinates,
// 3. build binary tree based on nearest neighbor-joining algorithm.
// 4. recursively visit each sub-branch and redo embed+tree (steps 2 and 3).
// Garry Paul Gippert.
*/

BNODE *bnode_embed_tree(int n, double **dmx)
{
	double **pos = embed_dmx(n, dmx);
	int *index = int_vector_ramp(n);
	BNODE *P = bnode_tree(pos, index, n, p_dim, dmx);
	int_vector_free(n, index);
	return P;
}

#define COMMAND_LINE_HELP "\n\n\
ACLUST  Generates pairwise alignments, distance matrix and phylogenetic tree from protein FASTA input files.\n\
Input entries may be pre-aligned, for example Fasta format output produced by MAFFT, or unaligned.\n\
\n\
Required:\n\
	-s <path>		filepath and name of substitution score matrix (e.g., '../dat/BLOSUM62.txt')\n\
Optional:\n\
	-p <string>		prefix for all output files (default=name of first input fasta file)\n\
	-e <char>		(D) distance tree only, (S) distance+single embed trees, (F, default) distance+single+full embed trees\n\
	-d <integer>		embed dimension (default 20)\n\
	-dmxfile <my.dmx>	Skips alignment phase and reads directly distance matrix in labelI labelJ DIJ\n\
	-go <float>		Gap open penalty\n\
	-ge <float>		Gap extend penalty\n\
Switches:\n\
	-maln 			read multiple alignment fasta\n\
	-metadata 		switch ON write node metadata in tree file\n\
	-jaln			switch OFF write alignment as JSON file\n\
	-taln			switch ON write alignment as text file\n\
	-baln			switch ON write alignment as binary file (not currently supported)\n\
	-nonself		deactivates self alignments\n\
	-v			activates more verbose output\n\
\n\
Input files are Fasta with >accession on one line and sequences on following until the next > is reached\n\
Output files share a prefix <p>, which is default name of first fasta input file\n\
	<p>_aln.txt	alignments, text\n\
	<p>_aln.js	alignments individual JSON rows (optional, activate using -j) \n\
	<p>_dmx.txt	distance matrix, text\n\
	<p>_dree.txt	distance matrix tree, newick text\n\
	<p>_tree0.txt	single embed tree, newick text\n\
	<p>_tree.txt	fulfull embed tree, newick text\n\
\n\
DETAILS: Score distance matrix based on pairwise local sequence\n\
alignments (Smith & Waterman) OR multiple alignment given in input\n\
Fasta records. Scoredist values (Sonnhammer & Hollich) are normalized\n\
to the shorter sequence length.  Tree computed from distance matrix\n\
by embedding into orthogonal coordinates (metric matrix distance\n\
geometry) and nearest-neighbor joining. Tree refined by re-embedding\n\
and neighbor-joining points in each sub-branch independently, and\n\
recursively. A pure distance NNJ tree is computed also.\n\
\n\
AUTHOR: Garry Paul Gippert, GarryG@dtu.dk, DTU Bioengineering\n\
"

void command_line_help(int c, int argc, char *argv[])
{
	fprintf(stderr, "%s [command_line_parameters_and_flags] my.fasta [another.fasta ...] %s",
		argv[0], COMMAND_LINE_HELP), exit(0);
}
void parameter_value_missing(int cstart, int argc, char *argv[])
{
	int c = cstart - 1;
	fprintf(stderr, "value needed after parameter [%d] '%s', try '%s -h' for HELP\n",
		c, argv[c], argv[0]), exit(1);
}
int pparse(int argc, char *argv[])
/* parse command line and set some input and output filename defaults */
{
	int c = 1;
	while (c < argc) {
		if (strncmp(argv[c], "-dmxfile", 8) == 0) {
			if (++c == argc)
				parameter_value_missing(c, argc, argv);
			f_dmxfilename = char_string(argv[c++]);
			fprintf(stderr, "f_dmxfilename set to '%s'\n", f_dmxfilename);
		}
		else if (strncmp(argv[c], "-h", 2) == 0) {
			++c;
			command_line_help(c, argc, argv);
		}
		else if (strncmp(argv[c], "-s", 2) == 0) {
			if (++c == argc)
				parameter_value_missing(c, argc, argv);
			scorematrixfile = char_string(argv[c++]);
			fprintf(stderr, "scorematrix file set to '%s'\n", scorematrixfile);
		}
		else if (strncmp(argv[c], "-p", 2) == 0) {
			if (++c == argc)
				parameter_value_missing(c, argc, argv);
			oprefix = char_string(argv[c++]);
			fprintf(stderr, "output prefix set to '%s'\n", oprefix);
		}
		else if (strncmp(argv[c], "-d", 2) == 0) {
			if (++c == argc)
				parameter_value_missing(c, argc, argv);
			if (sscanf(argv[c], "%d", &p_dim) == 1)
				fprintf(stderr, "dimension set to %d\n", p_dim);
			else
				fprintf(stderr, "Could not parse argv[%d] '%s'\n", c, argv[c]), exit(1);
			c++;
		}
		else if (strncmp(argv[c], "-e", 2) == 0) {
			if (++c == argc)
				parameter_value_missing(c, argc, argv);
			if (sscanf(argv[c], "%c", &p_e) == 1)
				fprintf(stderr, "Embed set to %c\n", p_e);
			else
				fprintf(stderr, "Could not parse argv[%d] '%s'\n", c, argv[c]), exit(1);
			c++;
		}
		else if (strncmp(argv[c], "-go", 3) == 0) {
			if (++c == argc)
				parameter_value_missing(c, argc, argv);
			if (sscanf(argv[c], "%lf", &p_go) == 1)
				fprintf(stderr, "Gap open set to %g\n", p_go);
			else
				fprintf(stderr, "Could not parse argv[%d] '%s'\n", c, argv[c]), exit(1);
			c++;
		}
		else if (strncmp(argv[c], "-ge", 3) == 0) {
			if (++c == argc)
				parameter_value_missing(c, argc, argv);
			if (sscanf(argv[c], "%lf", &p_ge) == 1)
				fprintf(stderr, "Gap extension set to %g\n", p_ge);
			else
				fprintf(stderr, "Could not parse argv[%d] '%s'\n", c, argv[c]), exit(1);
			c++;
		}
		else if (strncmp(argv[c], "-gx", 3) == 0) {
			if (++c == argc)
				parameter_value_missing(c, argc, argv);
			if (sscanf(argv[c], "%d", &p_gx) == 1)
				fprintf(stderr, "Gap max crossover length set to %d\n", p_gx);
			else
				fprintf(stderr, "Could not parse argv[%d] '%s'\n", c, argv[c]), exit(1);
			c++;
		}
		/* switch ON <-> OFF */
		else if (strncmp(argv[c], "-nonself", 8) == 0) {
			++c;
			p_nonself = (p_nonself + 1) % 2;
			fprintf(stderr, "nonself flag set to %d\n", p_nonself);
		}
		else if (strncmp(argv[c], "-metadata", 9) == 0) {
			++c;
			p_metadata = (p_metadata + 1) % 2;
			fprintf(stderr, "write node metadata %d\n", p_metadata);
		}
		else if (strncmp(argv[c], "-maln", 5) == 0) {
			++c;
			p_maln = (p_maln + 1)%2;
			fprintf(stderr, "Read multiple alignment %d\n", p_maln);
		}
		else if (strncmp(argv[c], "-baln", 5) == 0) {
			++c;
			p_baln = (p_baln + 1)%2;
			fprintf(stderr, "Write align binary %d\n", p_baln);
		}
		else if (strncmp(argv[c], "-jaln", 5) == 0) {
			++c;
			p_jaln = (p_jaln + 1)%2;
			fprintf(stderr, "Write align json %d\n", p_jaln);
		}
		else if (strncmp(argv[c], "-taln", 5) == 0) {
			++c;
			p_taln = (p_taln + 1)%2;
			fprintf(stderr, "Write align text %d\n", p_taln);
		}
		else if (strncmp(argv[c], "-v", 2) == 0) {
			++c;
			p_v++;
			fprintf(stderr, "verbose flag set to %d\n", p_v);
		}

		/* remaining '-' cases HERE */
		else if (strncmp(argv[c], "-", 1) == 0) {
			fprintf(stderr, "assumed parameter [%d] '%s' not recognized\n", c, argv[c]), exit(1);
		}
		/* terminate parameter parsing */
		else {
			fprintf(stderr, "assume termination of parameter stream\n");
			break;
		}

	}
	/* some validation */
	if (strchr("DSF", p_e) == NULL)
		fprintf(stderr, "Parameter -e %c invalid, must be D, S or F\n", p_e), exit(1);
	fprintf(stderr, "pparse returns %d\n", c);
	return (c);
}

void read_fasta_files(int argc, char *argv[], int cstart)
{
	fprintf(stderr, "%s Trying to read multiple fasta files, c start %d of %d\n", argv[0], cstart, argc);
	int c;
	for (c = cstart; c < argc; c++) {
		fprintf(stderr, " argv[%d], %s\n", c, argv[c]);
		read_fasta(argv[c]);
		fprintf(stderr, "read group %d %s, ns now %d\n", g_nsrc, argv[c], g_nent);
		g_nsrc++;
	}
}

void explain_input_better(int argc, char *argv[], int cstart)
{
	fprintf(stderr, "argc %d cstart %d\n", argc, cstart);
	if (cstart == argc)
		fprintf(stderr, "Input file(s) needed, expecting Fasta filenames\n"), exit(1);
}

/* usage
	double **dmx = read_dmx(f_dmxfilename);
*/

int facc_index(char *word)
{
	int i; for (i = 0; i < g_nent; i++) if (strcmp(facc[i], word)==0) return i;
	return -1;
}

double **read_dmx(char *filename)
{
/* remember MAXENTRIES 10000
	int g_nent = 0;
	char *facc[MAXENTRIES];
*/
	g_nent = 0;
	char line[MAXLINELEN], word1[MAXWORDLEN], word2[MAXWORDLEN];
	float value;
	int index1, index2;
	FILE *fp = fopen(filename, "r");
	if (!fp) fprintf(stderr, "Could not open filename %s r mode\n", filename), exit(1);
	while (fgets(line, sizeof(line), fp) != NULL)
	{
		if (strlen(line) == 0 || line[0] == '#')
			continue;
		if (sscanf(line, "%s %s %f", word1, word2, &value) != 3)
			fprintf(stderr, "Could not sscanf word1, word2, value line >>%s<<\n", line), exit(1);
		while((index1 = facc_index(word1)) < 0)
			facc[g_nent++] = char_string(word1);
		while((index2 = facc_index(word2)) < 0)
			facc[g_nent++] = char_string(word2);
		/* fprintf(stderr, "g_nent %d word1 %s index %d word2 %s index %d value %g\n", g_nent, word1, facc_index(word1), word2, facc_index(word2), value); */
	}
	fclose(fp);

	/* that was a very fast way to allocate the labels and dimension the space */
	global_dmx = double_matrix(g_nent, g_nent);
	int i, j;
	for (i = 0; i < g_nent; i++)
		for (j = i;  j < g_nent; j++)
			global_dmx[i][j] = global_dmx[j][i] = -99.0;

	/* Reopen the file and rescan the data. Simply faster than recording it the first time */
	fp = fopen(filename, "r");
	if (!fp) fprintf(stderr, "Could not open filename %s r mode the second time !!\n", filename), exit(1);
	while (fgets(line, sizeof(line), fp) != NULL)
	{
		if (strlen(line) == 0 || line[0] == '#')
			continue;
		if (sscanf(line, "%s %s %f", word1, word2, &value) != 3)
			fprintf(stderr, "Could not sscanf word1, word2, value line >>%s<<\n", line), exit(1);
		if((index1 = facc_index(word1)) < 0)
			fprintf(stderr, "Unexpected that word >>%s<< is not on facc list g_nent %d\n", word1, g_nent), exit(1);
		if((index2 = facc_index(word2)) < 0)
			fprintf(stderr, "Unexpected that word >>%s<< is not on facc list g_nent %d\n", word2, g_nent), exit(1);
		if (global_dmx[index1][index2] > -90.0)
			fprintf(stderr, "Unexpected duplicate word1 %s index1 %d word2 %s index2 %d g_nent %d\n",
				word1, index1, word2, index2, g_nent), exit(1);
		global_dmx[index1][index2] = global_dmx[index2][index1] = value;
	}
	fclose(fp);
	return (global_dmx);
}

int main(int argc, char *argv[])
{
	char *treefile = NULL;
	int c = pparse(argc, argv);
	if (!oprefix) {
		oprefix = char_string("this");
		fprintf(stderr, "oprefix initialized to %s\n", oprefix);
	}
	if (c == 1)
		explain_input_better(argc, argv, c);
	if (!oprefix)
		oprefix = char_string(argv[c]);
	fprintf(stderr, "oprefix set to %s\n", oprefix);

	double **dmx = NULL;

	if (!f_dmxfilename) {
		/* usual business of generating dmx */
		/* Read scorematrix before Fasta  - why ? maybe it was to initialize index lookup ? */
		if (!scorematrixfile)
			scorematrixfile = char_string(DEFAULT_SCORE_MATRIX);
		read_scorematrix(scorematrixfile);
		read_fasta_files(argc, argv, c);
		/* perform heavy lifting all pairs alignment, returning double **dmx */
		dmx = align_fasta();
		write_dmx(dmx, oprefix);
		fprintf(stderr, "Proceeding with distance matrix dimension %d created by alignment\n", g_nent);
	} else {
		/* read_dmx returns double ** but also sets g_nent and point labels facc */
		dmx = read_dmx(f_dmxfilename);
		fprintf(stderr, "Got DMX with %d entries, scanning for possible holes...\n", g_nent);
		int i, j, neg = 0;
		for (i = 0; i < g_nent; i++)
			for (j = 0; j < g_nent; j++)
				if ( dmx[i][j] < -0.0 ) {
					fprintf(stderr, "Negative i, j %d %d v %g\n", i, j, dmx[i][j]);
					neg++;
				}
		if (neg > 0)
			fprintf(stderr, "Cannot proceed without a complete distance matrix\n"), exit(1);
		fprintf(stderr, "Proceeding with distance matrix dimension %d read from file >>%s<<\n", g_nent, f_dmxfilename);
	}

	/* Construct tree based on distance matrix alone.  */
	BNODE *dree = bnode_distance_tree(g_nent, dmx);
	treefile = char_vector(strlen(oprefix) + strlen(".dree.txt") + 1);
	sprintf(treefile, "%s%s", oprefix, ".dree.txt");
	fprintf(stderr, "Distance tree %s\n", treefile);
	write_tree(dree, treefile);
	free(treefile);
	if (p_e == 'D')
		fprintf(stderr, "Halt after distance tree\n"), exit(0);

	/* Construct tree based on single pass embedding.
	*/
	BNODE *tree = bnode_embed_tree(g_nent, dmx);
	treefile = char_vector(strlen(oprefix) + strlen(".tree0.txt") + 1);
	sprintf(treefile, "%s%s", oprefix, ".tree0.txt");
	fprintf(stderr, "Single embed tree %s\n", treefile);
	write_tree(tree, treefile);
	free(treefile);
	if (p_e == 'S')
		fprintf(stderr, "Halt after single embed tree\n"), exit(0);

	/* Continue with full recursive embedding.
	*/
	tree->left = bnode_reembed(tree->left, 'L', dmx, g_nent, p_dim);
	tree->right = bnode_reembed(tree->right, 'R', dmx, g_nent, p_dim);
	treefile = char_vector(strlen(oprefix) + strlen(".tree.txt") + 1);
	sprintf(treefile, "%s%s", oprefix, ".tree.txt");
	fprintf(stderr, "Full recursive embed tree %s\n", treefile);
	write_tree(tree, treefile);
	free(treefile);

	char *binary_treefile = char_string("binary_treefile.btf");
	write_tree_binary(tree, binary_treefile);

	exit(0);
}
