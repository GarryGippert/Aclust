/* CCLUST - COLLAPSE CLUSTERS

Compute NNJ binary tree from distance matrix and generate iTol datasets that
collapse clusters.

Input file is a distance matrix represented by a text file with three columns:
	labelI labelI 0.0
	labelI labelJ distanceIJ
	labelJ labelJ 0.0
	...
Input should contain N different labels and all N*(N+1)/2 pairwise distances
(including self).

Command line parameters include:
	-p 'prefix'			= prefix used for output files

Output includes
	<prefix>.dree.txt file		= NNJ tree, Newick format.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXSEQUENCELEN 10000
#define MAXFILENAME 1000
#define MAXLINELEN 1000
#define MAXWORDLEN 100
#define	EPSILON		1.0e-8
#define NEARZERO(a)	(fabs(a) < 10.0*EPSILON ? 0.0 : (a))
#define SIGN(a)		( (a) < 0.0 ?   (-1.0) :    (1.0) )

int p_v = 0;			/* verbose flag, set to 1 for additional diagnostic output */
char *oprefix = NULL;

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

/* used in parsing text */
char line[MAXLINELEN], text[MAXLINELEN], acc[MAXLINELEN];

int p_strict = 0;		/* set to 1, and no deviation is allowed in the recomputation of alignment score */

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
	if (! B)
		return 0;
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
		bnode_print_tree(fp, B->left);
		fprintf(fp, ":%g", B->left_distance);
		if (PRINTNL)
			fprintf(fp, "\n");
		/* right */
		fprintf(fp, ",");
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

#define DMX_ONE		0x00000001  /* == 1 */
#define DMX_TWO		0x00000002  /* == 2 */

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

#define COMMAND_LINE_HELP "\n\n\
CCLUST  Collapse clusters from distance matrix and NNJ tree.\n\
Required parameters:\n\
	-d <my.dmx>		Filename of pairwise distances dmx.txt 'labelI labelJ DIJ'\n\
Optional parameters:\n\
	-p <string>		prefix for all output files (default=name of first input fasta file)\n\
	-e <char>		(D) distance tree only, (S) distance+single embed trees, (F, default) distance+single+full embed trees\n\
	-go <float>		Gap open penalty\n\
	-ge <float>		Gap extend penalty\n\
	-v			activates more verbose output\n\
\n\
Output files share a common prefix <p>, which is default name of first fasta input file\n\
	<p>_dree.txt	DMX NNJ binary tree, Newick format\n\
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
		if (strncmp(argv[c], "-h", 2) == 0) {
			++c;
			command_line_help(c, argc, argv);
		}
		else if (strncmp(argv[c], "-p", 2) == 0) {
			if (++c == argc)
				parameter_value_missing(c, argc, argv);
			oprefix = char_string(argv[c++]);
			fprintf(stderr, "output prefix set to '%s'\n", oprefix);
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
	fprintf(stderr, "Attempt read DMX from '%s'\n", filename);
	double **dmx = NULL;
	g_nent = 0;
	int i, j;
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
	if (! A)
		return UNDEF_DIS;
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
	if (p_v)
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
	if (p_v)
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

void bnode_printone(BNODE *B, int n, double **dmx)
{
	fprintf(stderr, "parent   %4d %4d %s %g\n",
		B->index, bnode_count(B), (B->index >= 0 ? facc[B->index] : ""), B->parent_distance);
	fprintf(stderr, "  left   %4d %4d %g <dL> %g\n",
		(B->left ? B->left->index : -1), bnode_count(B->left), B->left_distance, branch_distance_within(B->left, n, dmx));
	fprintf(stderr, " right   %4d %4d %g <dR> %g\n",
		(B->right ? B->right->index : -1), bnode_count(B->right), B->right_distance, branch_distance_within(B->right, n, dmx));
	fprintf(stderr, " between %4d %4d <dLR> %g\n", bnode_count(B->left), bnode_count(B->right), branch_distance_between(B->left, B->right, n, dmx));
}

void bnode_collapse(BNODE *B, int n, double **dmx, FILE *fp)
{
	int minc = 2;
	double mind = 50.0;
	if ( bnode_count(B) <= minc )
		return;
	if ( branch_distance_between(B->left, B->right, n, dmx) <= mind ) {
		fprintf(fp, "%s|%s\n", leftmost(B), rightmost(B));
		return;
	}
	bnode_collapse(B->left, n, dmx, fp);
	bnode_collapse(B->right, n, dmx, fp);
}

void bnode_abrl(BNODE *B, double branch, double *sum, double *cnt)
/* recursively compute sum of total branchlengths to leaf nodes, and count of leaf nodes */
{
	/* termination - not sure when we would ever need this... */
	if (! B)
		return;
	/* count leaf nodes and add parent distance to sum */
	if (B->index > -1) {
		if (p_v)
			fprintf(stderr, "leaf %4d %-20s branch %7.3f\n", B->index, facc[B->index], branch);
		*sum += branch;
		*cnt += 1.0;
		return;
	}
	bnode_abrl(B->left, branch + B->left_distance, &(*sum), &(*cnt));
	bnode_abrl(B->right, branch + B->right_distance, &(*sum), &(*cnt));
}

void bnode_collapse_abrl(BNODE *B, int n, double **dmx, FILE *fp)
{
	double sum = 0.0, cnt = 0.0, ave;
	bnode_abrl(B, 0.0, &sum, &cnt);
	if (cnt <= 1.0)
		return;
	if ((ave = sum/cnt) <= 50.0) {
		fprintf(fp, "# ABRL sum %g cnt %g ave %g\n", sum, cnt, ave);
		fprintf(fp, "%s|%s\n", leftmost(B), rightmost(B));
		return;
	}
	bnode_collapse_abrl(B->left, n, dmx, fp);
	bnode_collapse_abrl(B->right, n, dmx, fp);
}

int main(int argc, char *argv[])
{
	int c = pparse(argc, argv);
	int dmx_flag = DMX_ONE;
	char *f_dmx = char_string(argv[c]);

	/* read double **dmx, but also set global variables containing
	   point labels and count, char **facc and g_nent, respectively */
	double **dmx = read_dmx(f_dmx);

	/* construct DMX NNJ binary tree */
	int *index = int_vector_ramp(g_nent);
	fprintf(stderr, "bnode_tree_dmx dmx_flag=%d\n", dmx_flag);
	BNODE *dree = bnode_tree_dmx(g_nent, index, dmx, dmx_flag);

	/* write DMX NNJ binary tree */
	char *treefile = NULL;
	treefile = char_vector(strlen(oprefix) + strlen(".dree.txt") + 1);
	sprintf(treefile, "%s%s", oprefix, ".dree.txt");
	fprintf(stderr, "Distance tree %s\n", treefile);
	write_tree(dree, treefile);
	free(treefile);

	FILE *fp = stdout;
	fprintf(fp, "COLLAPSE\nDATA\n");
#ifdef METHOD1
	bnode_collapse(dree, g_nent, dmx, fp);
#else
	bnode_collapse_abrl(dree, g_nent, dmx, fp);
#endif
	int_vector_free(g_nent, index);
	exit(0);
}
