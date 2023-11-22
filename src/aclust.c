/* ACLUST

Aclust generates a phylogenetic tree from FASTA input.

Please contact the author Garry Paul Gippert, GarryG@dtu.dk, DTU Bioengineering,
Danish Technical University, with questions or for more information.

1. For each pair of input sequences, a distance is computed using a modified version of ScoreDist
(Sohnhammer & Hollich, 2005). Input sequences may represent a multiple alignment, OR, if the input
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
	prefix_aln.js	JSON-parsable details of each pairwise alignment
	prefix_dmx.txt	Distance matrix in plain text format <labelI> <labelJ> <distanceIJ>
	prefix_tree.txt	Newick-format tree

ACLUST was developed and written by Garry Paul Gippert, and packaged as a single C source file in 2023.

ACLUST is made available with a GNU General Public License v3.0, and may be used, modified and distributed
freely as long as license and copyleft notices are preserved.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXSEQUENCELEN 10000
#define MAXFILENAME 1000
#define MAXLINELEN 1000
#define DEFAULT_SCORE_MATRIX "dat/BLOSUM62.dat"

#define ALIGN_GAP_CHAR '-'
#define ALIGN_GAP_VAL -99.9
#define ALIGN_GAP_NDX -1
#define ALIGN_PAD_CHAR '+'
#define ALIGN_PAD_VAL -99.9
#define ALIGN_PAD_NDX -2
#define ALIGN_GAP	0x0001	/* enable GAP */
#define ALIGN_PAD	0x0002	/* enable PAD */
#define ALIGN_CROSS	0x0004	/* enable gap cross-over */
#define ALIGN_BODY	0x0008	/* enable gap body-gaps (deprecated) */

#define	EPSILON		1.0e-8
#define NEARZERO(a)	(fabs(a) < 10.0*EPSILON ? 0.0 : (a))
#define SIGN(a)		( (a) < 0.0 ?   (-1.0) :    (1.0) )

int p_v = 0;			/* verbose flag, set to 1 for additional diagnostic output */
double p_fg = 12.0;		/* first gap penalty */
double p_ng = 1.0;		/* next gap penalty */
int p_json = 1;			/* output alignments in pseudo-json format */
int p_nonself = 0;		/* do not align with self (show only off-diagonal elements */


void j_osb(FILE * fp)
/* output jsonified open square bracket */
{
	fprintf(fp, "[");
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
	fprintf(fp, "{\"key\": \"%s\", \"value\": %d", key, value);
	j_cmt(fp, comment);
	j_url(fp, url);
	fprintf(fp, "}");
	if (comma == YES)
		fprintf(fp, ",");
}

void j_str(FILE * fp, int comma, char *key, char *value, char *comment, char *url)
/* output jsonified string */
{
	fprintf(fp, "{\"key\": \"%s\", \"value\": \"%s\"", key, value);
	j_cmt(fp, comment);
	j_url(fp, url);
	fprintf(fp, "}");
	if (comma == YES)
		fprintf(fp, ",");
}

void j_dbl(FILE * fp, int comma, char *key, double value, char *comment, char *url)
/* output jsonified string */
{
	fprintf(fp, "{\"key\": \"%s\", \"value\": %g", key, value);
	j_cmt(fp, comment);
	j_url(fp, url);
	fprintf(fp, "}");
	if (comma == YES)
		fprintf(fp, ",");
}

void j_csb(FILE * fp)
/* output jsonified close square bracket */
{
	fprintf(fp, "]\n");
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
#define MAXENTRIES 3000
char *facc[MAXENTRIES];
char *fseq[MAXENTRIES];
int fsrc[MAXENTRIES];		/* source fasta file names, worst case one sequence per file */
/* int  *frti[MAXENTRIES]; residue index deactivated */

/* Global parameters start with p_ */
int p_m = 0;			/* multiple alignment input */

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
		if (sscanf(line, "%[^\n]", text) != 1)
			fprintf(stderr, "Cannot sscanf %%[^\n]:%s\n", line), exit(1);
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

double **pair_score_matrix(int fi, int fj)
/* provide Blosum62 substitution score matrix for pair of fasta elements fi, fj */
{
	char *si = fseq[fi], *sj = fseq[fj];
	int i, j, ni = strlen(si), nj = strlen(sj);
	double **m = double_matrix(ni, nj);
	for (i = 0; i < ni; i++)
		for (j = 0; j < nj; j++)
			m[i][j] = scorematrix_element(si[i], sj[j]);
	return (m);
}

double align_score(int n1, int n2, double **S, double **M, int *o1, int *o2, int flag)
/* Generate optimal local alignment path using affine gap penalties
// return alignment score, indirectly return sequence offsets, and filled-in match matrix
// attributed to Smith & Waterman, 1981
// Note: Pointers provide an optimization/obfuscation trade-off. */
{
	int i, j;
	double **T, **U, **V, *Ti, *Tp, *Si, *Vi, *Vp, *Ui, *Mi, *Up, t1, t2, t3, Tij, Uij, Vij;
	double ascore = 0.0;
	*o1 = -1;
	*o2 = -1;

	/* Initialize transition matrices */
	T = double_matrix(n1 + 1, n2 + 1);
	U = double_matrix(n1 + 1, n2 + 1);
	V = double_matrix(n1 + 1, n2 + 1);
	for (j = n2 - 1; j >= 0; j--) {
		Tij = T[n1][j + 1];
		T[n1][j] = Tij;
		U[n1][j] = Tij;
		V[n1][j] = Tij - p_fg;
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
		Ui[n2] = Ti[n2] - p_fg;
		Vi[n2] = Ti[n2];

		/* backwards from last column */
		for (j = n2 - 1; j >= 0; j--) {
			Tij = Tp[j + 1] + Si[j];
			Mi[j] = Tij;
			/* insertion in sequence 2 */
			t1 = Vp[j] - p_ng;
			t2 = Tp[j] - p_fg;
			Vij = (t1 > t2 ? t1 : t2);
			if (flag & ALIGN_CROSS) {
				t3 = Up[j] - p_ng;
				Vij = (t3 > Vij ? t3 : Vij);
			}
			if (Vij > Tij)
				Tij = Vij;
			/* insertion in sequence 1 */
			t1 = Ui[j + 1] - p_ng;
			t2 = Ti[j + 1] - p_fg;
			Uij = (t1 > t2 ? t1 : t2);
			if (flag & ALIGN_CROSS) {
				t3 = Vi[j + 1] - p_ng;
				Uij = (t3 > Uij ? t3 : Uij);
			}
			if (Uij > Tij)
				Tij = Uij;
			if (Tij > ascore) {
				if (p_v)
					printf("backprop I %d J %d Ascore now %g, delta %g\n", i, j, Tij, Tij - ascore);
				ascore = Tij;
				(*o1) = i;
				(*o2) = j;
			}
			if (Tij < 0.0)
				Tij = 0.0;
			Ti[j] = Tij;
			Vi[j] = Vij;
			Ui[j] = Uij;
		}
	}
	double_matrix_free(n1 + 1, n2 + 1, T);
	double_matrix_free(n1 + 1, n2 + 1, U);
	double_matrix_free(n1 + 1, n2 + 1, V);
	/* it is possible this ascore is larger that computed from the generated alignment strings. */
	return (ascore);
}

int align_index(int n1, int n2, double **S, double **M, int o1, int o2, int **d1, int **d2, int align_flag)
/* revisit optimal path through alignment match matrix, and deduce sequence indices
// return total length of (optionally padded) alignment, return indirectly sequence indices */
{
	int k, l, nk, nl, bk, bl, tg, *t1, *t2;
	double max, tmp;
	int plen = 0;
	*d1 = NULL;
	*d2 = NULL;
	if (o1 == -1 || o2 == -1)
		return (plen);
	t1 = int_vector(n1 + n2);
	t2 = int_vector(n1 + n2);

	/* pad at alignment start */
	if (align_flag & ALIGN_PAD) {
		for (k = 0; k < o1; k++) {
			t1[plen] = k;
			t2[plen] = ALIGN_PAD_NDX;
			plen++;
		}
		for (l = 0; l < o2; l++) {
			t1[plen] = ALIGN_PAD_NDX;
			t2[plen] = l;
			plen++;
		}
	}

	/* first aligned position */
	k = o1;
	l = o2;

	while (1) {
		/* previous match state */
		t1[plen] = k;
		t2[plen] = l;
		plen++;

		/* next k,l */
		nk = k + 1;
		nl = l + 1;
		if (nk >= n1 || nl >= n2)
			break;

		/* match state */
		bk = nk;
		bl = nl;
		max = M[nk][nl];
		if (max <= 0.0)
			break;

		/* gap in column ? */
		for (k = nk + 1; k < n1; k++) {
			tmp = M[k][nl] - p_fg - p_ng * (k - nk - 1);
			if (tmp > max) {
				max = tmp;
				bk = k;
			}
		}

		/* gap in row ? */
		for (l = nl + 1; l < n2; l++) {
			tmp = M[nk][l] - p_fg - p_ng * (l - nl - 1);
			if (tmp > max) {
				max = tmp;
				bl = l;
			}
		}
		/* gap in both row and column ? */
		if (align_flag & ALIGN_BODY)
			for (k = nk; k < n1; k++) {
				for (l = nl; l < n2; l++) {
					tg = k - nk + l - nl;
					tmp = M[k][l] - (tg ? p_fg - p_ng * (tg - 1) : 0.0);
					if (tmp > max) {
						max = tmp;
						bk = k;
						bl = l;
					}
				}
			}

		if (align_flag & ALIGN_GAP) {
			/* gap in column */
			for (k = nk; k < bk; k++) {
				t1[plen] = k;
				t2[plen] = ALIGN_GAP_NDX;
				plen++;
			}
			/* gap in row */
			for (l = nl; l < bl; l++) {
				t1[plen] = ALIGN_GAP_NDX;
				t2[plen] = l;
				plen++;
			}
		}
		/* last position */
		k = bk;
		l = bl;
	}

	/* optional pad C-terminal with non-aligned regions */
	if (align_flag & ALIGN_PAD) {
		for (k = nk; k < n1; k++) {
			t1[plen] = k;
			t2[plen] = ALIGN_PAD_NDX;
			plen++;
		}
		for (l = nl; l < n2; l++) {
			t1[plen] = ALIGN_PAD_NDX;
			t2[plen] = l;
			plen++;
		}
	}

	/* indirect return condensed list of alignment indices */
	*d1 = int_vector(plen);
	*d2 = int_vector(plen);
	memcpy((*d1), t1, (plen) * sizeof(int));
	memcpy((*d2), t2, (plen) * sizeof(int));
	int_vector_free(n1 + n2, t1);
	int_vector_free(n1 + n2, t2);
	return (plen);
}

void align_strings(char *name1, int n1, char *name2, int n2, double **S, double **M, int o1, int o2, char *s1, char *s2, char **a1, char **a2, int align_flag)
/* produce padded alignment strings *a1 and *a2 */
{
	int k, l, nk, nl, bk, bl, i, plen = 0;
	double max, tmp;
	char *t1, *t2;
	if (p_v)
		printf("Alignment strings between %s (%d) and %s (%d)\n", name1, n1, name2, n2);
	(*a1) = NULL;
	(*a2) = NULL;
	if (o1 == -1 || o2 == -1)
		return;

	/* tmp space for accumulating alignment strings */
	t1 = char_vector(n1 + n2 + 2);
	t2 = char_vector(n1 + n2 + 2);

	/* optional pad N-terminal with non-aligned regions */
	if (align_flag & ALIGN_PAD) {
		for (k = 0; k < o1; k++) {
			t1[plen] = s1[k];
			t2[plen] = ALIGN_PAD_CHAR;
			plen++;
		}
		for (l = 0; l < o2; l++) {
			t1[plen] = ALIGN_PAD_CHAR;
			t2[plen] = s2[l];
			plen++;
		}
	}

	/* first match */
	k = o1;
	l = o2;

	while (1) {
		/* align state, plen is 'padded' aka full length alignment strings */
		t1[plen] = s1[k];
		t2[plen] = s2[l];
		double scoreelement = scorematrix_element(t1[plen], t2[plen]);
		if (p_v)
			printf("aln plen %d K %3d%c  vs L  %3d%c,  score %g\n",
			       plen, k, s1[k], l, s2[l], scoreelement);
		plen++;

		/* next k,l */
		nk = k + 1;
		nl = l + 1;
		if (nk >= n1 || nl >= n2)
			break;

		/* match state */
		bk = nk;
		bl = nl;
		max = M[nk][nl];
		if (max <= 0.0)
			break;

		/* gap in column ? */
		for (i = nk + 1; i < n1; i++) {
			tmp = M[i][nl] - p_fg - p_ng * (i - nk - 1);
			if (tmp > max) {
				max = tmp;
				bk = i;
			}
		}

		/* gap in row ? */
		for (i = nl + 1; i < n2; i++) {
			tmp = M[nk][i] - p_fg - p_ng * (i - nl - 1);
			if (tmp > max) {
				max = tmp;
				bl = i;
			}
		}

		if (align_flag & ALIGN_GAP) {
			/* insert column gaps in strings */
			for (i = nk; i < bk; i++) {
				t1[plen] = s1[i];
				t2[plen] = ALIGN_GAP_CHAR;
				plen++;
			}

			/* insert row gaps in strings */
			for (i = nl; i < bl; i++) {
				t1[plen] = ALIGN_GAP_CHAR;
				t2[plen] = s2[i];
				plen++;
			}
		}

		/* last best position */
		k = bk;
		l = bl;
	}

	/* optional pad C-terminal with non-aligned regions */
	if (align_flag & ALIGN_PAD) {
		for (i = nk; i < n1; i++) {
			t1[plen] = s1[i];
			t2[plen] = ALIGN_PAD_CHAR;
			plen++;
		}
		for (i = nl; i < n2; i++) {
			t1[plen] = ALIGN_PAD_CHAR;
			t2[plen] = s2[i];
			plen++;
		}
	}

	/* politely terminate strings */
	t1[plen] = 0;
	t2[plen] = 0;

	/* return first alignment in string of length plen */
	(*a1) = char_string(t1);
	free(t1);

	/* return second alignment in string of length plen */
	(*a2) = char_string(t2);
	free(t2);
}

#define MAXSCOREDIST 9999.9
double compute_scoredistance(double ma, double mr, double m1, double m2, double scale)
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

FILE *alnfp = NULL;		/* file pointer for writing alignment free text */
FILE *jsnfp = NULL;		/* file pointer for writing alignment JSON line-by-line */

int p_strict = 0;		/* set to 1, and no deviation is allowed in the recomputation of alignment score */

double align_stats(char *name1, char *a1, char *name2, char *a2, int expected_plen, double expected_ascore)
/* produce alignment statistics and return score distance */
{
	//plen = count total padded length of alignment
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
	int plen, alen, mlen, ilen, glen, olen, clen, nlen, o1, o2, minlen;
	double ascore, mscore, mscore1, mscore2, mscorer, scale, sd0, sd1, sd2, sd;

	if (a1 == NULL || a2 == NULL)
		fprintf(stderr, "alignment string a1 is NULL or a2 is NULL\n"), exit(1);
	if ((plen = strlen(a1)) != strlen(a2))
		fprintf(stderr, "strlen(a1) %ld != strlen(a2) %ld\n%s\n%s\n", strlen(a1), strlen(a2), a1, a2), exit(1);
	if (expected_plen > 0 && plen != expected_plen)
		fprintf(stderr, "plen %d != expected_plen %d\n", plen, expected_plen), exit(1);

	/* initialize counting and score varibles */
	alen = mlen = ilen = glen = olen = clen = nlen = 0;
	mscore = mscore1 = mscore2 = mscorer = 0.0;
	o1 = o2 = -1;

	/* compute gaps and match scores from the alignment strings */
	int i, isg = 0, n1 = 0, n2 = 0;	/* sequence positions n1 and n2 */
	double s;
	char a, b;

	double sum_score = 0.0, gap_score;
	for (i = 0; i < plen; i++) {
		a = a1[i];
		b = a2[i];
		if (a != ALIGN_PAD_CHAR && a != ALIGN_GAP_CHAR)
			n1++;
		if (b != ALIGN_PAD_CHAR && b != ALIGN_GAP_CHAR)
			n2++;
		s = scorematrix_element(a, b);
		if (fabs(s) < 100)
			sum_score += s;
		if (p_v)
			printf("p_fg %g x olen %d + p_ng %g x (glen %d - olen %d)\n", p_fg, olen, p_ng, glen, olen);
		gap_score = p_fg * (float)(olen) + p_ng * (float)(glen - olen);
		if (p_v)
			printf("stats plen %d : '%d%c' vs '%d%c' :  score %g, sumscore %g, gapscore %g, total_score %g\n",
			       i, n1, a, n2, b, s, sum_score, gap_score, (sum_score - gap_score));
		/* ignore pad positions (outside of the local alignment) */
		if (a == ALIGN_PAD_CHAR || b == ALIGN_PAD_CHAR)
			continue;
		/* aligned positions */
		/* set offsets */
		if (o1 < 0)
			o1 = n1, o2 = n2;
		alen++;
		if (a == ALIGN_GAP_CHAR || b == ALIGN_GAP_CHAR) {
			glen++;
			if (!isg)
				olen++;
			isg = 1;
			continue;
		}
		/* matched positions */
		mlen++;
		if (a == b)
			ilen++;
		if (s > 0)
			clen++;
		if (s >= 0)
			nlen++;
		mscore += s;
		mscore1 += scorematrix_element(a, a);
		mscore2 += scorematrix_element(b, b);
		mscorer -= 1.0;	/* random match expectation */
		isg = 0;
	}
	int error = 0;
	double gapcost = olen * p_fg + (glen - olen) * p_ng;
	ascore = mscore - gapcost;
	if (expected_ascore > 0.0 && ascore != expected_ascore) {
		if (p_v)
			fprintf(stderr, "Pair %s %s computed ascore %g != expected %g\n", name1, name2, ascore, expected_ascore);
		if (p_strict)
			fprintf(stderr, "Exiting in p_strict mode\n"), exit(1);
		double g = sqrt((ascore - expected_ascore) * (ascore - expected_ascore)) / (0.5 * (ascore + expected_ascore));
		if (p_v)
			fprintf(stderr, "Pair %s %s, g = |D|/A %g, ascore %g, expected %g, D %g, A %g\n",
				name1, name2, g, ascore, expected_ascore, (ascore - expected_ascore), 0.5 * (ascore + expected_ascore));
#define MAXG 0.01
		if (g > MAXG)
			if (p_v)
				fprintf(stderr, "Error g %g > MAXG %g\n", g, MAXG), error++;
	}

	/* compute score distances */
	/* original score distance of Sonhammer & Hollich normalized to alignment length */
	scale = (double)alen / (double)mlen;
	sd0 = compute_scoredistance(mscore, mscorer, mscore1, mscore2, scale);

	/* normalize to first sequence length */
	scale = (double)n1 / (double)mlen;
	sd1 = compute_scoredistance(mscore, mscorer, mscore1, mscore2, scale);

	/* normalize to second sequence length */
	scale = (double)n2 / (double)mlen;
	sd2 = compute_scoredistance(mscore, mscorer, mscore1, mscore2, scale);

	/* preferred score distance to minimum length */
	minlen = (n1 < n2 ? n1 : n2);
	scale = (double)minlen / (double)mlen;
	sd = compute_scoredistance(mscore, mscorer, mscore1, mscore2, scale);

	/* write alignments */
	if (p_json) {
		/* JSON format parsable line-by-line */
		j_osb(jsnfp);
		/* input */
		j_str(jsnfp, YES, "name1", name1, NULL, NULL);
		j_int(jsnfp, YES, "len1", n1, NULL, NULL);
		j_str(jsnfp, YES, "name2", name2, NULL, NULL);
		j_int(jsnfp, YES, "len2", n2, NULL, NULL);
		/* alignment strings */
		j_str(jsnfp, YES, "aln1", a1, NULL, NULL);
		j_int(jsnfp, YES, "o1", o1, NULL, NULL);
		j_str(jsnfp, YES, "aln2", a2, NULL, NULL);
		j_int(jsnfp, YES, "o2", o2, NULL, NULL);
		/* alignment counts and scores */
		j_int(jsnfp, YES, "plen", plen, NULL, NULL);
		j_int(jsnfp, YES, "alen", alen, NULL, NULL);
		j_int(jsnfp, YES, "mlen", mlen, NULL, NULL);
		j_int(jsnfp, YES, "ilen", ilen, NULL, NULL);
		j_int(jsnfp, YES, "glen", glen, NULL, NULL);
		j_int(jsnfp, YES, "olen", olen, NULL, NULL);
		j_int(jsnfp, YES, "clen", clen, NULL, NULL);
		j_int(jsnfp, YES, "nlen", nlen, NULL, NULL);
		j_dbl(jsnfp, YES, "ascore", ascore, NULL, NULL);
		j_dbl(jsnfp, YES, "mscore", mscore, NULL, NULL);
		j_dbl(jsnfp, YES, "mscore1", mscore1, NULL, NULL);
		j_dbl(jsnfp, YES, "mscore2", mscore2, NULL, NULL);
		j_dbl(jsnfp, YES, "mscorer", mscorer, NULL, NULL);
		/* score distances */
		j_dbl(jsnfp, YES, "sd0", sd0, NULL, NULL);
		j_dbl(jsnfp, YES, "sd1", sd1, NULL, NULL);
		j_dbl(jsnfp, YES, "sd2", sd2, NULL, NULL);
		j_dbl(jsnfp, NO, "sd", sd, NULL, NULL);	/* last data element receives a NO to solve a json-related
							   issue */
		j_csb(jsnfp);
		fflush(jsnfp);
	}
	if (alnfp) {
		/* free-form text */
		fprintf(alnfp, "Align %s %d x %s %d", name1, n1, name2, n2);
		fprintf(alnfp, " O1 %d O2 %d", o1, o2);
		fprintf(alnfp, " Plen %d Alen %d Mlen %d Ilen %d Glen %d Olen %d Clen %d Nlen %d",
			plen, alen, mlen, ilen, glen, olen, clen, nlen);
		fprintf(alnfp, " Ascore %f Mscore %f", ascore, mscore);
		fprintf(alnfp, " M1 %g M2 %g MR %g", mscore1, mscore2, mscorer);
		fprintf(alnfp, " SD0 %g SD1 %g SD2 %g SD %g\n", sd0, sd1, sd2, sd);
		fprintf(alnfp, "%s\n%s\n", a1, a2);
		fprintf(alnfp, "\n");	/* extra newline for human readability */
		fflush(alnfp);
	}
	if (error) {
		fprintf(stderr, "Error detected in align_stats - check most recent alignments\n");
		if (p_strict)
			exit(1);
	}
	return (sd);
}

double pair_malign(int fi, int fj)
/* from pre-aligned pair of fasta elements fi, fj
// return ScoreDist normalized to shortest of the two sequences
// write alignment info to stdout in pseudo-JSON */
{
	/* pre-aligned fasta sequences */
	char *a1 = fseq[fi], *a2 = fseq[fj];
	double d = align_stats(facc[fi], a1, facc[fj], a2, -1, -1.0);
	return (d);
}

double pair_align(int fi, int fj)
/* align pair of fasta elements fi, fj
// return ScoreDist normalized to shortest of the two sequences
// write alignment info to stdout in pseudo-JSON */
{
	/* going into the alignment: */
	int align_flag = ALIGN_GAP | ALIGN_PAD | ALIGN_CROSS;
	char *s1 = fseq[fi], *s2 = fseq[fj], *a1 = NULL, *a2 = NULL;
	int o1, o2, n1 = strlen(s1), n2 = strlen(s2), *d1 = NULL, *d2 = NULL;

	double **sx = pair_score_matrix(fi, fj);	/* score matrix */
	double **mx = double_matrix(n1, n2);	/* match matrix */

	/* compute optimal alignment and alignment score */
	double ascore = align_score(n1, n2, sx, mx, &o1, &o2, align_flag);

	/* revisit alignment and compute total padded alignment length */
	int plen = align_index(n1, n2, sx, mx, o1, o2, &d1, &d2, align_flag);

	/* generate alignment strings */
	align_strings(facc[fi], n1, facc[fj], n2, sx, mx, o1, o2, s1, s2, &a1, &a2, align_flag);

	/* compute alignment statistics */
	double d = align_stats(facc[fi], a1, facc[fj], a2, plen, ascore);

	/* memory free */
	free(a1);
	free(a2);
	int_vector_free(plen, d1);
	int_vector_free(plen, d2);
	double_matrix_free(n1, n2, sx);
	double_matrix_free(n1, n2, mx);

	return (d);
}


double **align_fasta()
/* return all-vs-all pairwise distance matrix from multiple alignment or computed pairwise alignments
// return a symmetric matrix containing align scoredistance normalized
// to the shortest sequence length of the pair */
{
	int i, inext, j;
	double **dmx = double_matrix(g_nent, g_nent);
	for (i = 0; i < g_nent; i++)
		for (j = (p_nonself ? i + 1 : i); j < g_nent; j++)
			dmx[i][j] = dmx[j][i] = (p_m ? pair_malign(i, j) : pair_align(i, j));
	return (dmx);
}

/* EMBED */

int p_dim = 20;			/* embed dimension: default 20 */
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
char p_e = 'F';			/* recursive reembed: 'F' full, 'S' single */

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
	BNODE *B;
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
/* formerly known as bnode_length */
{
	if (B->index >= 0)
		return 1;
	else if (B->left == NULL || B->right == NULL)
		fprintf(stderr, "B->left is NULL or B->right is NULL\n"), exit(1);
	return bnode_count(B->left) + bnode_count(B->right);
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
	}
}

void bnode_print(FILE * fp, BNODE * B)
{
	bnode_print_tree(fp, B);
	fprintf(fp, ":0;");	/* terminate tree */
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

/* provide a binary tree from a position matrix with dimensions n x dim, using nearest-neighbor joining algorithm */
BNODE *bnode_tree(double **pos, int *index, int n, int dim)
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
	P = bnode_tree(pos, index, n, dim);

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

BNODE *bnode_recursive_embed(int n, double **dmx)
/* Main Aclust method starts here:
// 1. embed distance matrix into orthogonal coordinates,
// 2. build binary tree based on nearest neighbor-joining algorithm.
// 3. recursively visit each sub-branch and redo embed+tree (steps 1 and 2).
// Garry Paul Gippert.
*/
{
	int dim = p_dim;
	double **pos = embed_dmx(n, dmx);
	int *index = int_vector_ramp(n);
	BNODE *P = bnode_tree(pos, index, n, dim);
	if (p_e == 'F') {
		fprintf(stderr, "Full recursive embed/cluster\n");
		P->left = bnode_reembed(P->left, 'L', dmx, n, dim);
		P->right = bnode_reembed(P->right, 'R', dmx, n, dim);
	}
	else {
		fprintf(stderr, "Single pass embed/cluster\n");
	}
	int_vector_free(n, index);
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
}

void write_tree_binary(BNODE * P, char *filename)
{
	printf("You are attempting to write_tree_binary to %s, unsuccessfully! treelength %d\n", filename, bnode_count(P));
	fprintf(stderr, "You are attempting to write_tree_binary %s, unsuccessfully! treelength %d\n", filename, bnode_count(P));
	exit(1);
}

void write_dmx(double **dmx, char *filename)
/* print distance upper half matrix plus diagonal */
{
	FILE *fp;
	int i, j;
	if ((fp = fopen(filename, "w")) == NULL)
		fprintf(stderr, "Distance file %s cannot be opened for writing\n", filename), exit(1);
	for (i = 0; i < g_nent; i++)
		for (j = i; j < g_nent; j++)
			fprintf(fp, "%s %s %f\n", facc[i], facc[j], dmx[i][j]);
	fclose(fp);
}

char *scorematrixfile = NULL;
char *oprefix = NULL;
char *alnfile = NULL;
char *jsnfile = NULL;
char *dmxfile = NULL;
char *treefile = NULL;


#define COMMAND_LINE_HELP "\n\n\
ACLUST  Computes pairwise alignments, a distance matrix, and a phylogenetic tree from\n\
protein FASTA input files. Input sequences may be pre-aligned, representing a multiple alignment.\n\
Otherwise sequence alignments are computed on the fly.\n\
\n\
Optional parameters:\n\
	-s <input_scorematrix_file>	possibly '$BIA/dat/BLOSUM62.txt'\n\
	-p <output_prefix>		<first_input_filename_including_dotfa>\n\
\n\
Optional flags:\n\
	-m 	activates	M.ultiple alignment mode\n\
	-j  (de)activates	J.SON output file\n\
	-v	activates	V.erbose program output\n\
	-nonself	activates	Non-self, therefore only show off-diagonal alignments\n\
\n\
Output files share a prefix <p>, which is default name of first fasta input file\n\
	<p>.aln.txt	alignments, text\n\
	<p>.aln.js	alignments individual JSON rows (optional, activate using -j) \n\
	<p>.dmx.txt	distance matrix, text\n\
	<p>.tree.txt	Newick-format phylogenetic tree, text\n\
\n\
DETAILS: Score distance matrix based on pairwise local sequence\n\
alignments (Smith & Waterman) OR multiple alignment given in input\n\
Fasta records. Scoredist values (SohnHammer & Hollich) are normalized\n\
to the shorter sequence length.  Tree computed from distance matrix\n\
by embedding into orthogonal coordinates (metric matrix distance\n\
geometry) and nearest-neighbor joining. Tree refined by re-embedding\n\
and neighbor-joining points in each sub-branch independently, and\n\
recursively.\n\
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
		if (strncmp(argv[c], "-h", 2) == 0) {
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
		else if (strncmp(argv[c], "-nonself", 8) == 0) {
			++c;
			p_nonself = (p_nonself + 1) % 2;
			fprintf(stderr, "nonself flag set to %d\n", p_nonself);
		}
		else if (strncmp(argv[c], "-j", 2) == 0) {
			++c;
			p_json = (p_json + 1) % 2;
			fprintf(stderr, "json flag set to %d\n", p_json);
		}
		else if (strncmp(argv[c], "-m", 2) == 0) {
			++c;
			p_m = 1;
			fprintf(stderr, "multiple align set to %d\n", p_m);
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
	fprintf(stderr, "pparse returns %d\n", c);
	return (c);
}

void output_initialization()
/*  default file openings after parsing, but before doing anything else */
{
	if (!oprefix)
		oprefix = char_string("this");
	fprintf(stderr, "oprefix initialized to %s\n", oprefix);

	if (!alnfile) {
		alnfile = char_vector(strlen(oprefix) + strlen(".aln.txt") + 1);
		sprintf(alnfile, "%s%s", oprefix, ".aln.txt");
	}
	fprintf(stderr, "alnfile %s\n", alnfile);
	if ((alnfp = fopen(alnfile, "w")) == NULL)
		fprintf(stderr, "Unable to open aln file %s for writing\n", alnfile), exit(1);

	/* optional write JSON text */
	if (p_json) {
		if (!jsnfile) {
			jsnfile = char_vector(strlen(oprefix) + strlen(".aln.js") + 1);
			sprintf(jsnfile, "%s%s", oprefix, ".aln.js");
		}
		fprintf(stderr, "jsnfile %s\n", jsnfile);
		if ((jsnfp = fopen(jsnfile, "w")) == NULL)
			fprintf(stderr, "Unable to open JSON file %s for writing\n", jsnfile), exit(1);
	}

	if (!dmxfile) {
		dmxfile = char_vector(strlen(oprefix) + strlen(".dmx.txt") + 1);
		sprintf(dmxfile, "%s%s", oprefix, ".dmx.txt");
	}
	fprintf(stderr, "dmxfile %s\n", dmxfile);

	if (!treefile) {
		treefile = char_vector(strlen(oprefix) + strlen(".tree.txt") + 1);
		sprintf(treefile, "%s%s", oprefix, ".tree.txt");
	}
	fprintf(stderr, "treefile %s\n", treefile);

	/* Here is where we initialize additional output files with a standard prefix */
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

int main(int argc, char *argv[])
{
	int c = pparse(argc, argv);
	if (c == 1)
		explain_input_better(argc, argv, c);

	/* Read scorematrix before Fasta */

	if (!scorematrixfile)
		scorematrixfile = char_string(DEFAULT_SCORE_MATRIX);
	read_scorematrix(scorematrixfile);

	read_fasta_files(argc, argv, c);

	if (!oprefix)
		oprefix = char_string(argv[c]);
	fprintf(stderr, "oprefix set to %s\n", oprefix);

	output_initialization();

	double **dmx = align_fasta();

	write_dmx(dmx, dmxfile);

	BNODE *tree = bnode_recursive_embed(g_nent, dmx);

	write_tree(tree, treefile);

	char *binary_treefile = char_string("binary_treefile.btf");

	write_tree_binary(tree, binary_treefile);

	exit(0);
}
