/* ACLUST

COLLAPSE CLUSTERS

Input file is a distance matrix text file with three columns representing N points and N*(N+1)/2 distances.

	labelI labelI 0.0
	labelI labelJ distanceIJ
	labelJ labelJ 0.0
	etc.

Output files produced: (all with shared prefix)
	prefix_dree.txt	Distance matrix NNJ tree, Newick format.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXSEQUENCELEN 10000
#define MAXFILENAME 1000
#define MAXLINELEN 1000
#define MAXWORDLEN 100
#define DEFAULT_SCORE_MATRIX "dat/BLOSUM62.dat"

#define	EPSILON		1.0e-8
#define NEARZERO(a)	(fabs(a) < 10.0*EPSILON ? 0.0 : (a))
#define SIGN(a)		( (a) < 0.0 ?   (-1.0) :    (1.0) )

int p_v = 0;			/* verbose flag, set to 1 for additional diagnostic output */
double p_go = 12.0;		/* gap open = first gap penalty */
double p_ge = 1.0;		/* gap extend = next gap penalty */
int p_json = 1;			/* output alignments in pseudo-json format */
int p_nonself = 0;		/* do not align with self (show only off-diagonal elements */

char *f_dmx = NULL;
char *oprefix = NULL;

FILE *jsnfp = NULL;		/* file pointer for writing alignment JSON line-by-line */
FILE *alnfp = NULL;		/* file pointer for writing alignment free text */

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

/* SECTION FASTA */

int g_nent = 0, g_nsrc = 0;
#define MAXENTRIES 10000
char *facc[MAXENTRIES];
/* int  *frti[MAXENTRIES]; residue index deactivated */

/* Global parameters start with p_ */
int p_m = 0;			/* multiple alignment input */

/* used in parsing text */
char line[MAXLINELEN], text[MAXLINELEN], acc[MAXLINELEN], seq[MAXSEQUENCELEN];

int p_strict = 0;		/* set to 1, and no deviation is allowed in the recomputation of alignment score */

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

/* traverse binary tree
*/

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
	fprintf(stderr, "bnode_distance_tree dmx_flag=%d\n", dmx_flag);
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
	BNODE *P = bnode_tree(pos, index, n, p_dim);
	int_vector_free(n, index);
	return P;
}

#define COMMAND_LINE_HELP "\n\n\
CCLUST  Collapse clusters from distance matrix and NNJ tree.\n\
Required parameters:\n\
	-d <my.dmx>		Filename of pairwise distances dmx.txt 'labelI labelJ DIJ'\n\
Optional parameters:\n\
	-p <string>		prefix for all output files (default=name of first input fasta file)\n\
	-e <char>		(D) distance tree only, (S) distance+single embed trees, (F, default) distance+single+full embed trees\n\
	-go <float>		Gap open penalty\n\
	-ge <float>		Gap extend penalty\n\
Optional flags:\n\
	-m 			activates to interpret input Fasta as MSA\n\
Less important flags:\n\
	-j			deactivates writing of JSON alignment file\n\
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
Fasta records. Scoredist values (SohnHammer & Hollich) are normalized\n\
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
		if (strncmp(argv[c], "-d", 2) == 0) {
			if (++c == argc)
				parameter_value_missing(c, argc, argv);
			f_dmx = char_string(argv[c++]);
			fprintf(stderr, "f_dmx set to '%s'\n", f_dmx);
		}
		else if (strncmp(argv[c], "-h", 2) == 0) {
			++c;
			command_line_help(c, argc, argv);
		}
		else if (strncmp(argv[c], "-p", 2) == 0) {
			if (++c == argc)
				parameter_value_missing(c, argc, argv);
			oprefix = char_string(argv[c++]);
			fprintf(stderr, "output prefix set to '%s'\n", oprefix);
		}
		else if (strncmp(argv[c], "-e", 2) == 0) {
			if (++c == argc)
				parameter_value_missing(c, argc, argv);
			if (sscanf(argv[c], "%c", &p_e) == 1)
				fprintf(stderr, "Embed set to %d\n", p_e);
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
		/* remaining '-' cases */
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

	/* set some variables */
	if (!oprefix) {
		oprefix = char_string("this");
		fprintf(stderr, "oprefix initialized to %s\n", oprefix);
	}
	if (!oprefix)
		oprefix = char_string(argv[c]);
	fprintf(stderr, "oprefix set to %s\n", oprefix);
	return (c);
}

int facc_index(char *word)
{
	int i; for (i = 0; i < g_nent; i++) if (strcmp(facc[i], word)==0) return i;
	return -1;
}

double **read_dmx(char *filename)
/* usage double **dmx = read_dmx(f_dmx); */
/* remember MAXENTRIES 10000 int g_nent = 0; char *facc[MAXENTRIES]; */
{
	double **dmx = NULL;
	int i, j;
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
	fprintf(stderr, "Parsed %d labels from DMX file %s\n", g_nent, filename);

	/* that was a very fast way to allocate the labels and dimension the space */
	dmx = double_matrix(g_nent, g_nent);
	for (i = 0; i < g_nent; i++)
		for (j = i;  j < g_nent; j++)
			dmx[i][j] = dmx[j][i] = -99.0;

	/* Reopen the file and rescan the data. Simply faster than recording it the first time */
	fp = fopen(filename, "r");
	if (!fp) fprintf(stderr, "Could not open filename %s r mode the second time !!\n", filename), exit(1);
	int cnt = 0;
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
		if (dmx[index1][index2] > -90.0)
			fprintf(stderr, "Unexpected duplicate word1 %s index1 %d word2 %s index2 %d g_nent %d\n",
				word1, index1, word2, index2, g_nent), exit(1);
		dmx[index1][index2] = dmx[index2][index1] = value;
		cnt += 1;
	}
	fclose(fp);
	fprintf(stderr, "Parsed %d matrix values from DMX file %s, expected N*(N+1)/2 %d\n", cnt, filename, g_nent*(g_nent+1)/2);

	/* examine negative values - we do not proceed with an incomplete DMX */
	int neg = 0;
	for (i = 0; i < g_nent; i++)
		for (j = 0; j < g_nent; j++)
			if ( dmx[i][j] < -0.0 ) {
				fprintf(stderr, "Negative i, j %d %d v %g\n", i, j, dmx[i][j]);
				neg++;
			}
	if (neg > 0)
		fprintf(stderr, "Cannot proceed without a complete distance matrix\n"), exit(1);
	fprintf(stderr, "DMX ready\n");
	return (dmx);
}

int *bnode_index_vector(BNODE *B, int *n)
/* return vector of leaf node indices and indirect vector length */
{
	*n = bnode_count(B); /* number of leaves in this subtree */
	int *index = int_vector(*n), i = 0;
	bnode_indexi(B, index, &i);
	if (*n != i)
		fprintf(stderr, "*n %d != i %d\n", *n, i), exit(1);
	return(index);
}

#define UNDEF_DIS -99.99
double branch_distance_within(BNODE *A, int n, double **dmx)
/* compute average dmx distance among leaf nodes in branche A */
{
	int a, *indexa, i, j;
	indexa = bnode_index_vector(A, &a);
	double dis, sum = 0.0, sqr = 0.0, cnt = 0.0, ave=UNDEF_DIS, rms=UNDEF_DIS;
	for (i = 0; i < a; i++)
		for (j = i+1; j < a; j++) {
			dis = dmx[indexa[i]][indexa[j]]; 
			sum += dis;
			sqr += dis*dis;
			cnt += 1.0;
		}
	if (cnt > 0.0) {
		ave = sum/cnt;
		rms = sqrt(sqr/cnt - ave*ave);
	}
	fprintf(stderr, "Na %d Cnt %g Ave %g Rms %g\n", a, cnt, ave, rms);
	int_vector_free(a, indexa);
	return ave;
}

double branch_distance_between(BNODE *A, BNODE *B, int n, double **dmx)
/* compute average dmx distance between leaf nodes in branches A vs B */
{
	int a, b, *indexa, *indexb, i, j;
	indexa = bnode_index_vector(A, &a);
	indexb = bnode_index_vector(B, &b);
	double dis, sum = 0.0, sqr = 0.0, cnt = 0.0, ave=UNDEF_DIS, rms=UNDEF_DIS;
	for (i = 0; i < a; i++)
		for (j = 0; j < b; j++) {
			dis = dmx[indexa[i]][indexb[j]]; 
			sum += dis;
			sqr += dis*dis;
			cnt += 1.0;
		}
	if (cnt > 0.0) {
		ave = sum/cnt;
		rms = sqrt(sqr/cnt - ave*ave);
	}
	fprintf(stderr, "Na %d Nb %d Cnt %g Ave %g Rms %g\n", a, b, cnt, ave, rms);
	int_vector_free(a, indexa);
	int_vector_free(b, indexb);
	return ave;
}

char *leftmost(BNODE *B)
{
	if (B->left)
		return(leftmost(B->left));
	return facc[B->index];
}

char *rightmost(BNODE *B)
{
	if (B->right)
		return(rightmost(B->right));
	return facc[B->index];
}

void bnode_split(BNODE *B, int n, double **dmx)
{
	printf("bnode index %d %g : left I %d, D %g, N %d, AD %g ; right I %d, D %g, N %d, AD %g ; distance %g\n",
		B->index,
		B->parent_distance,
		B->left->index,
		B->left_distance,
		bnode_count(B->left),
		branch_distance_within(B->left, n, dmx),
		B->right->index,
		B->right_distance,
		bnode_count(B->right),
		branch_distance_within(B->right, n, dmx),
		branch_distance_between(B->left, B->right, n, dmx));
	printf("extreme %s|%s\n", leftmost(B->left), rightmost(B->left));
	printf("extreme %s|%s\n", leftmost(B->right), rightmost(B->right));
	exit(0);
}

int main(int argc, char *argv[])
{
	int c = pparse(argc, argv);

	/* read_dmx returns double ** but also sets g_nent and point labels facc */
	double **dmx = read_dmx(f_dmx);

	/* binary NNJ tree */
	BNODE *dree = bnode_distance_tree(g_nent, dmx);

	/* write binary NNJ tree */
	char *treefile = NULL;
	treefile = char_vector(strlen(oprefix) + strlen(".dree.txt") + 1);
	sprintf(treefile, "%s%s", oprefix, ".dree.txt");
	fprintf(stderr, "Distance tree %s\n", treefile);
	write_tree(dree, treefile);
	free(treefile);

	bnode_split(dree, g_nent, dmx);

	exit(0);

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
