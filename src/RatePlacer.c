//  Created by Rasmus Nielsen on 11/30/21.
//  RatePlacer.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include "tools.h"
#include "opt.h"
#include <errno.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <fcntl.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#pragma GCC diagnostic ignored "-Wstringop-overflow"
#pragma GCC diagnostic ignored "-Wunused-result"

#define MINBL 0.0000001
#define MAXBL 5.0
#define NUMCAT 4
#define INF DBL_MAX
#define VERBOSE 0
#define testmax(a, b) \
	({ __typeof__ (a) _a = (a); \
	 __typeof__ (b) _b = (b); \
	 _a > _b ? _a : _b; })

#define INITIAL_BUFFER_SIZE 1000 // Initial buffer size, adjust as needed

// assignAges holds the max possible age for the assignment node of a read! So nodeage + bls of assignment node for that read
// double LRVEC[4][4], RRVEC[4][4], RRVAL[4], PMAT[3][NUMCAT][4][4];
double **LRVEC, **RRVEC, **RRVAL; //, PMAT[3][NUMCAT][4][4];;
_Thread_local double PMAT[3 * NUMCAT * 4 * 4];  // I think
double **statevector, **FRACLIKE, **nodeages, **bls, ***readlike, testAge, **pi, **par, *maxAges, totMaxAge, *assignAges, errorTest, rooted, curAgeBound;
unsigned long int numbase, numquery, queryagesknown;
int ***DATA, **QUERYDATA, *assignments, **nodeOrder, *usedReads, *treeAssign, *usedTrees, toMerge;
unsigned long int *readlength, *startpos;
_Thread_local int onDindic = 0; // hack to avoid passing this indicator around
unsigned long int *numseq, *treeRoots;
long numTrees;

// TO DO: Try to make this not global in the future
//int curRead, curTree;

// temps for processing of reads
int ***readsTreeSorted, *readOrder, **tempAssignments, **numReadsPerAssign;
unsigned long int *numbases, **readLengthTemp, **startposTemp, tempnumquery;
double ****readLikeTemp;

_Thread_local int tip, comma = 0; /*globals used to read in the tree. Old code - don't ask.*/

_Thread_local FILE *infile;
FILE *outfile, *readsfile = NULL;

/* ---- gzip-transparent file helpers ---- */
static _Thread_local int   infile_is_pipe = 0;
static _Thread_local char *infile_membuf  = NULL;

/* Open for reading; .gz files decompressed via gzip -dc pipe (not seekable). */
static FILE *xfopen(const char *path)
{
	size_t n = strlen(path);
	if (n > 3 && memcmp(path + n - 3, ".gz", 3) == 0) {
		char cmd[4096];
		snprintf(cmd, sizeof(cmd), "gzip -dc -- \"%s\"", path);
		infile_is_pipe = 1;
		infile_membuf  = NULL;
		return popen(cmd, "r");
	}
	infile_is_pipe = 0;
	infile_membuf  = NULL;
	return fopen(path, "r");
}

static void xfclose(FILE *f)
{
	if (infile_is_pipe) pclose(f); else fclose(f);
	infile_is_pipe = 0;
}

/* Like xfopen but always seekable: .gz files are decompressed into memory
   and returned as an fmemopen stream. Use for files needing fgetpos/fsetpos. */
static FILE *xfopen_seekable(const char *path)
{
	size_t n = strlen(path);
	if (n > 3 && memcmp(path + n - 3, ".gz", 3) == 0) {
		char cmd[4096];
		snprintf(cmd, sizeof(cmd), "gzip -dc -- \"%s\"", path);
		FILE *gz = popen(cmd, "r");
		if (!gz) return NULL;
		size_t cap = 65536, used = 0;
		char *buf = malloc(cap);
		char tmp[8192]; size_t nr;
		while ((nr = fread(tmp, 1, sizeof(tmp), gz)) > 0) {
			if (used + nr > cap) { cap = (cap + nr) * 2; buf = realloc(buf, cap); }
			memcpy(buf + used, tmp, nr);
			used += nr;
		}
		pclose(gz);
		infile_membuf  = buf;
		infile_is_pipe = 0;
		return fmemopen(buf, used, "r");
	}
	infile_is_pipe = 0;
	infile_membuf  = NULL;
	return fopen(path, "r");
}

static void xfclose_seekable(FILE *f)
{
	fclose(f);
	free(infile_membuf);
	infile_membuf = NULL;
}
/* ---- end gzip helpers ---- */

// Adding tree information so that we can do node uncertainity!!
struct node
{
	int up[2]; // for storing children node number, assuming binary. -1 means no children (current node is leaf)
	int down;  // for storing parent node number
	double bl; // Scaled branch length to parent
};

#define MERGE_MODE_FRACTION 0
#define MERGE_MODE_BP 1

double merge_coverage_threshold = 0.05; // Default threshold
int merge_coverage_mode = MERGE_MODE_FRACTION; // Default mode (fraction)

struct node **trees;

/*subfunction needed by 'getclade*/
void linknodes(int i, int j, int nodee, unsigned long int treeNum) /*linking i down to nodee and j down to nodee*/
{
	trees[treeNum][nodee].up[0] = j;
	trees[treeNum][nodee].up[1] = i;
	trees[treeNum][i].down = nodee;
	trees[treeNum][j].down = nodee;
}

/*subfunction needed by 'getclade*/
/** What does this subfunction do exactly? **/
// is tip used to track number of leaves observed while reading in the file?
// 	if so, should there be a check that the number of tips cannot be larger than numleaves?
int specsearch(int numleaves, int treeNum)
{
	char ch;
	int i = 1;
	ch = fgetc(infile);
	if ((ch != ')') && (ch != '(') && (ch != ',') && (ch != ' ') && (ch != '\t') && (ch != '\n') && (ch != EOF))
	{
		ungetc(ch, infile);
		(void)fscanf(infile, "%d", &tip);
		while ((ch = fgetc(infile)) != ':' && (i < 10))
			i++;
		trees[treeNum][tip + numleaves - 2].up[0] = -1;
		trees[treeNum][tip + numleaves - 2].up[1] = -1;
		while ((ch = (fgetc(infile))) == ' ')
			;
		ungetc(ch, infile);
		(void)fscanf(infile, "%lf", &trees[treeNum][tip + numleaves - 2].bl);
		/*printf("\nbranchlength of node %i =%f",tip+numleaves-1,tree[tip+numleaves-2].bl);*/
		return 1;
	}
	else
	{
		ungetc(ch, infile);
		return 0;
	}
}

/*subfunction needed by 'getclade*/
/** Gets node number... but how? **/
int getnodenumb()
{
	char c;
	int i, j = 0;
	fpos_t position;
	i = 0;
	fgetpos(infile, &position);
	do
	{
		c = fgetc(infile);
		if (c == ',')
			i++;
		if (c == '(')
			j = j - 1;
		if (c == ')')
			j++;
	} while ((j < 0) && (c != EOF));
	fsetpos(infile, &position);
	return (i + comma + 1);
}

/*some old code for reading a Newick tree*/
int getclade(int numleaves, int treeNum)
{
	int n1, n2, n3;
	char ch;

	do
	{
		if (specsearch(numleaves, treeNum) == 1)
		{
			/*tip++;*/
			return tip + (numleaves - 1);
		}
		ch = fgetc(infile);
		if (ch == ',')
		{
			comma++;
		}
		if (ch == ')')
		{
			if ((ch = fgetc(infile)) != ':')
			{
				ungetc(ch, infile);
			}
			else
			{
				do
				{
					ch = (fgetc(infile));
				} while ((ch == '\n') || (ch == ' '));
				ungetc(ch, infile);
				(void)fscanf(infile, "%lf", &trees[treeNum][n3 - 1].bl);
			}
			// returns root node of the (sub)tree
			return n3;
		}
		if (ch == '(')
		{
			n3 = getnodenumb();
			// Recursive function down each subtree
			n1 = getclade(numleaves, treeNum);
			n2 = getclade(numleaves, treeNum);
			linknodes(n1 - 1, n2 - 1, n3 - 1, treeNum);
		}
	} while (ch != ';');

	return -1;
}

void allocatetreememmory(int numleaves, unsigned long int treeNum)
{
	int i;

	// Number of nodes based off full binary tree
	trees[treeNum] = malloc((numleaves * 2 - 1) * (sizeof(struct node)));
}

void freetreememmory(void)
{
	for (unsigned long int i = 0; i < numTrees; i++)
	{
		if (usedTrees[i] != 0)
		{
			free(trees[i]);
		}
	}
	free(trees);
	free(usedTrees);
}

/*old code for printing a tree, treeNum = -1 for all trees*/
void printtree(int numleaves, int root, unsigned long int treeNum)

{
	int i;

	if (treeNum == -1)
	{
		for (unsigned long int j = 0; j < numTrees; j++)
		{
			if (usedTrees != 0)
			{
				printf("\nPRINTING TREE %lu\n", treeNum);
				for (i = 0; i < 2 * numleaves - 1; i++)
				{
					if (trees[j][i].up[0] != -1)
						printf("Node %i: up: (%i, %i) down: %i", i, trees[j][i].up[0], trees[j][i].up[1], trees[j][i].down);
					else
						printf(" Node %i (leaf): up: (%i, %i) down: %i", i, trees[j][i].up[0], trees[j][i].up[1], trees[j][i].down);
					if (i != root)
						printf(" (bl: %f)\n", trees[j][i].bl);
					else
						printf(" (root)\n");
				}
			}
		}
	}
	else
	{
		printf("\nPRINTING TREE %lu\n", treeNum);
		for (i = 0; i < 2 * numleaves - 1; i++)
		{
			if (trees[treeNum][i].up[0] != -1)
				printf("Node %i: up: (%i, %i) down: %i", i, trees[treeNum][i].up[0], trees[treeNum][i].up[1], trees[treeNum][i].down);
			else
				printf(" Node %i (leaf): up: (%i, %i) down: %i", i, trees[treeNum][i].up[0], trees[treeNum][i].up[1], trees[treeNum][i].down);
			if (i != root)
				printf(" (bl: %f)\n", trees[treeNum][i].bl);
			else
				printf(" (root)\n");
		}
	}
}

/*Code getting node to leaf length assuming ultrametric (all paths the same length) tree*/
double getMaxAge(int curNode, unsigned long int treeNum)
{
	if (trees[treeNum][curNode].up[0] == -1)
	{
		return (trees[treeNum][curNode].bl);
	}
	else
	{
		return (trees[treeNum][curNode].bl + getMaxAge(trees[treeNum][curNode].up[0], treeNum));
	}
}

// get gfl node number of children when given ratePlacer assignment node encoding
void getGFLChildren(int curNode, int childNodes[3], unsigned long int treeNum)
{
	// leaf, no children of node to check assignment
	if (curNode < numseq[treeNum])
	{
		childNodes[1] = -1;
		childNodes[2] = -1;
	}
	else // non-leaf
	{
		// convert curNode to gfl node
		curNode = curNode - numseq[treeNum];

		if (trees[treeNum][trees[treeNum][curNode].up[0]].up[0] == -1)
		{
			// Children are leaf
			childNodes[1] = trees[treeNum][curNode].up[0] - numseq[treeNum] + 1;
		}
		else
		{
			// children are internal
			childNodes[1] = trees[treeNum][curNode].up[0] + numseq[treeNum];
		}

		if (trees[treeNum][trees[treeNum][curNode].up[1]].up[0] == -1)
		{
			// Children are leaf
			childNodes[2] = trees[treeNum][curNode].up[1] - numseq[treeNum] + 1;
		}
		else
		{
			// children are internal
			childNodes[2] = trees[treeNum][curNode].up[1] + numseq[treeNum];
		}
	}
}

// get gfl nodes of parent and sibling when given ratePlacer assignment node
void getGFLParSib(int curNode, int parSib[3], unsigned long int treeNum)
{
	if (curNode < numseq[treeNum])
	{
		// leaf node
		curNode = curNode + numseq[treeNum] - 1;
	}
	else
	{
		// internal node
		curNode = curNode - numseq[treeNum];

		// if current node is root, no parent or sibling
		if (curNode == treeRoots[treeNum])
		{
			parSib[1] = -1;
			parSib[2] = -1;
		}
	}

	// get parent and convert, parent is always internal
	parSib[1] = trees[treeNum][curNode].down + numseq[treeNum];

	// get sibling
	int par = trees[treeNum][curNode].down;

	if (trees[treeNum][par].up[0] == curNode)
	{
		// parSib[2] = tree[par].up[1] - numseq + 1;
		if (trees[treeNum][trees[treeNum][par].up[1]].up[0] == -1)
		{
			// Children are leaf
			parSib[2] = trees[treeNum][par].up[1] - numseq[treeNum] + 1;
		}
		else
		{
			// children are internal
			parSib[2] = trees[treeNum][par].up[1] + numseq[treeNum];
		}
	}
	else
	{
		// parSib[2] = tree[par].up[0] - numseq + 1;
		if (trees[treeNum][trees[treeNum][par].up[0]].up[0] == -1)
		{
			// Children are leaf
			parSib[2] = trees[treeNum][par].up[0] - numseq[treeNum] + 1;
		}
		else
		{
			// children are internal
			parSib[2] = trees[treeNum][par].up[0] + numseq[treeNum];
		}
	}
}

int getGFLPar(int curNode, int treeNum)
{
	if (curNode < numseq[treeNum])
	{
		// leaf node
		curNode = curNode + numseq[treeNum] - 1;
	}
	else
	{
		// internal node
		curNode = curNode - numseq[treeNum];

		// Node is root
		if (curNode == treeRoots[treeNum])
		{
			return -1;
		}
	}

	// Parent is always internal node
	return (trees[treeNum][curNode].down + numseq[treeNum]);
}

//******************FUNCTIONS TO DEAL WITH THE GAMMA DISTRIBUION********************************************
double LnGamma(double alpha)
{
	/* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.
	   Stirling's formula is used for the central polynomial part of the procedure.
	   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
	   Communications of the Association for Computing Machinery, 9:684
	   */
	double x = alpha, f = 0, z;

	if (x < 7)
	{
		f = 1;
		z = x - 1;
		while (++z < 7)
			f *= z;
		x = z;
		f = -log(f);
	}
	z = 1 / (x * x);
	return f + (x - 0.5) * log(x) - x + .918938533204673 + (((-.000595238095238 * z + .000793650793651) * z - .002777777777778) * z + .083333333333333) / x;
}

double IncompleteGamma(double x, double alpha, double ln_gamma_alpha)
{
	/* returns the incomplete gamma ratio I(x,alpha) where x is the upper
	   limit of the integration and alpha is the shape parameter.
	   returns (-1) if in error
	   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
	   (1) series expansion     if (alpha>x || x<=1)
	   (2) continued fraction   otherwise
	   RATNEST FORTRAN by
	   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
19: 285-287 (AS32)
*/
	int i;
	double p = alpha, g = ln_gamma_alpha;
	double accurate = 1e-20, overflow = 1e30;
	double factor, gin = 0, rn = 0, a = 0, b = 0, an = 0, dif = 0, term = 0, pn[6];

	if (x == 0)
		return (0);
	if (x < 0 || p <= 0)
		return (-1);

	factor = exp(p * log(x) - x - g);
	if (x > 1 && x >= p)
		goto l30;
	/* (1) series expansion */
	gin = 1;
	term = 1;
	rn = p;
l20:
	rn++;
	term *= x / rn;
	gin += term;

	if (term > accurate)
		goto l20;
	gin *= factor / p;
	goto l50;
l30:
	/* (2) continued fraction */
	a = 1 - p;
	b = a + x + 1;
	term = 0;
	pn[0] = 1;
	pn[1] = x;
	pn[2] = x + 1;
	pn[3] = x * b;
	gin = pn[2] / pn[3];
l32:
	a++;
	b += 2;
	term++;
	an = a * term;
	for (i = 0; i < 2; i++)
		pn[i + 4] = b * pn[i + 2] - an * pn[i];
	if (pn[5] == 0)
		goto l35;
	rn = pn[4] / pn[5];
	dif = fabs(gin - rn);
	if (dif > accurate)
		goto l34;
	if (dif <= accurate * rn)
		goto l42;
l34:
	gin = rn;
l35:
	for (i = 0; i < 4; i++)
		pn[i] = pn[i + 2];
	if (fabs(pn[4]) < overflow)
		goto l32;
	for (i = 0; i < 4; i++)
		pn[i] /= overflow;
	goto l32;
l42:
	gin = 1 - factor * gin;

l50:
	/*printf("Incompletegamma got %f %f %f and returned %f\n",x,  alpha,  ln_gamma_alpha,gin);*/
	// printf("");
	return (gin);
}

double PointNormal(double prob)
{
	/* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
	   returns (-9999) if in error
	   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
	   Applied Statistics 22: 96-97 (AS70)

	   Newer methods:
	   Wichura MJ (1988) Algorithm AS 241: the percentage points of the
	   normal distribution.  37: 477-484.
	   Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage
	   points of the normal distribution.  26: 118-121.

*/
	double a0 = -.322232431088, a1 = -1, a2 = -.342242088547, a3 = -.0204231210245;
	double a4 = -.453642210148e-4, b0 = .0993484626060, b1 = .588581570495;
	double b2 = .531103462366, b3 = .103537752850, b4 = .0038560700634;
	double y, z = 0, p = prob, p1;

	p1 = (p < 0.5 ? p : 1 - p);
	if (p1 < 1e-20)
		return (-9999);

	y = sqrt(log(1 / (p1 * p1)));
	z = y + ((((y * a4 + a3) * y + a2) * y + a1) * y + a0) / ((((y * b4 + b3) * y + b2) * y + b1) * y + b0);
	return (p < 0.5 ? -z : z);
}

double PointChi2(double prob, double v)
{
	/* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
	   returns -1 if in error.   0.000002<prob<0.999998
	   RATNEST FORTRAN by
	   Best DJ & Roberts DE (1975) The percentage points of the
	   Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
	   Converted into C by Ziheng Yang, Oct. 1993.
	   */
	double e = .5e-6, aa = .6931471805, p = prob, g;
	double xx, c, ch, a = 0, q = 0, p1 = 0, p2 = 0, t = 0, x = 0, b = 0, s1, s2, s3, s4, s5, s6;

	if (p < .000002 || p > .999998 || v <= 0)
		return (-1);

	g = LnGamma(v / 2);
	xx = v / 2;
	c = xx - 1;
	if (v >= -1.24 * log(p))
		goto l1;

	ch = pow((p * xx * exp(g + xx * aa)), 1 / xx);
	if (ch - e < 0)
		return (ch);
	goto l4;
l1:
	if (v > .32)
		goto l3;
	ch = 0.4;
	a = log(1 - p);
l2:
	q = ch;
	p1 = 1 + ch * (4.67 + ch);
	p2 = ch * (6.73 + ch * (6.66 + ch));
	t = -0.5 + (4.67 + 2 * ch) / p1 - (6.73 + ch * (13.32 + 3 * ch)) / p2;
	ch -= (1 - exp(a + g + .5 * ch + c * aa) * p2 / p1) / t;
	if (fabs(q / ch - 1) - .01 <= 0)
		goto l4;
	else
		goto l2;

l3:
	x = PointNormal(p);
	p1 = 0.222222 / v;
	ch = v * pow((x * sqrt(p1) + 1 - p1), 3.0);
	if (ch > 2.2 * v + 6)
		ch = -2 * (log(1 - p) - c * log(.5 * ch) + g);
l4:
	q = ch;
	p1 = .5 * ch;
	if ((t = IncompleteGamma(p1, xx, g)) < 0)
	{
		printf("\nerr IncompleteGamma");
		return (-1);
	}
	p2 = p - t;
	t = p2 * exp(xx * aa + g + p1 - c * log(ch));
	b = t / ch;
	a = 0.5 * t - b * c;

	s1 = (210 + a * (140 + a * (105 + a * (84 + a * (70 + 60 * a))))) / 420;
	s2 = (420 + a * (735 + a * (966 + a * (1141 + 1278 * a)))) / 2520;
	s3 = (210 + a * (462 + a * (707 + 932 * a))) / 2520;
	s4 = (252 + a * (672 + 1182 * a) + c * (294 + a * (889 + 1740 * a))) / 5040;
	s5 = (84 + 264 * a + c * (175 + 606 * a)) / 2520;
	s6 = (120 + c * (346 + 127 * c)) / 5040;
	ch += t * (1 + 0.5 * t * s1 - b * c * (s1 - b * (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))));
	if (fabs(q / ch - 1) > e)
		goto l4;

	return (ch);
}

// #define PointGamma(prob, alpha, beta) PointChi2(prob, 2.0 * (alpha)) / (2.0 * (beta))

double CDFfunGamma(double x, double par[2])

{
	return IncompleteGamma(par[1] * x, par[0], LnGamma(par[0]));
}

//******************END OF FUNCTIONS TO DEAL WITH THE GAMMA DISTRIBUION********************************************

/***********************************************************
 *  This eigen() works for eigenvalue/vector analysis
 *         for real general square matrix A
 *         A will be destroyed
 *         rr,ri are vectors containing eigenvalues
 *         vr,vi are matrices containing (right) eigenvectors
 *
 *              A*[vr+vi*i] = [vr+vi*i] * diag{rr+ri*i}
 *
 *  Algorithm: Handbook for Automatic Computation, vol 2
 *             by Wilkinson and Reinsch, 1971double times[3], double parameters[7]
 *             most of source codes were taken from a public domain
 *             solftware called MATCALC.
 *  Credits:   to the authors of MATCALC
 *
 *  return     -1 not converged
 *              0 no complex eigenvalues/vectors
 *              1 complex eigenvalues/vectors
 *  Tianlin Wang at University of Illinois
 *  Thu May  6 15:22:31 CDT 1993
 ***************************************************************/

#define BASE 2	   /* base of floating point arithmetic */
#define DIGITS 40  /* no. of digits to the base BASE in the fraction */
#define MAXITER 30 /* max. no. of iterations to converge */

#define pos(i, j, n) ((i) * (n) + (j))

/*A is the matrix, rr = root real(nx1), ri = root imaginary(nx1), vr = real part of eigenvector(nxn), w = working space (size of 2n), set job to 1 (decides if both eigenvectors and eigenvalues should be calculated*/

int eigen(int job, double A[], int n, double rr[], double ri[],
		  double vr[], double vi[], double w[]);
void balance(double mat[], int n, int *low, int *hi, double scale[]);
void unbalance(int n, double vr[], double vi[], int low, int hi,
			   double scale[]);
int realeig(int job, double mat[], int n, int low, int hi, double valr[],
			double vali[], double vr[], double vi[]);
void elemhess(int job, double mat[], int n, int low, int hi,
			  double vr[], double vi[], int work[]);

int eigen(int job, double A[], int n, double rr[], double ri[],
		  double vr[], double vi[], double work[])
{
	/*  double work[n*2]: working space
	 */
	int low, hi, i, j, k, it, istate = 0;
	double tiny = sqrt(pow((double)BASE, (double)(1 - DIGITS))), t;

	balance(A, n, &low, &hi, work);
	elemhess(job, A, n, low, hi, vr, vi, (int *)(work + n));
	if (-1 == realeig(job, A, n, low, hi, rr, ri, vr, vi))
		return (-1);
	if (job)
		unbalance(n, vr, vi, low, hi, work);

	/* sort, added by Z. Yang */
	for (i = 0; i < n; i++)
	{
		for (j = i + 1, it = i, t = rr[i]; j < n; j++)
			if (t < rr[j])
			{
				t = rr[j];
				it = j;
			}
		rr[it] = rr[i];
		rr[i] = t;
		t = ri[it];
		ri[it] = ri[i];
		ri[i] = t;
		for (k = 0; k < n; k++)
		{
			t = vr[k * n + it];
			vr[k * n + it] = vr[k * n + i];
			vr[k * n + i] = t;
			t = vi[k * n + it];
			vi[k * n + it] = vi[k * n + i];
			vi[k * n + i] = t;
		}
		if (fabs(ri[i]) > tiny)
			istate = 1;
	}

	return (istate);
}

/* complex funcctions
 */

complex compl(double re, double im)
{
	complex r;

	r.re = re;
	r.im = im;
	return (r);
}

/*complex conj (complex a)
  {
  a.im = -a.im;
  */

// #define csize(a) (fabs(a.re) + fabs(a.im))

complex cplus(complex a, complex b)
{
	complex c;
	c.re = a.re + b.re;
	c.im = a.im + b.im;
	return (c);
}

complex cminus(complex a, complex b)
{
	complex c;
	c.re = a.re - b.re;
	c.im = a.im - b.im;
	return (c);
}

complex cby(complex a, complex b)
{
	complex c;
	c.re = a.re * b.re - a.im * b.im;
	c.im = a.re * b.im + a.im * b.re;
	return (c);
}

complex cdiv(complex a, complex b)
{
	double ratio, den;
	complex c;

	if (fabs(b.re) <= fabs(b.im))
	{
		ratio = b.re / b.im;
		den = b.im * (1 + ratio * ratio);
		c.re = (a.re * ratio + a.im) / den;
		c.im = (a.im * ratio - a.re) / den;
	}
	else
	{
		ratio = b.im / b.re;
		den = b.re * (1 + ratio * ratio);
		c.re = (a.re + a.im * ratio) / den;
		c.im = (a.im - a.re * ratio) / den;
	}
	return (c);
}

/*complex cexp (complex a)
  {
  complex c;
  c.re = exp(a.re);
  if (fabs(a.im)==0) c.im = 0;
  else  { c.im = c.re*sin(a.im); c.re*=cos(a.im); }
  return (c);
  }*/

complex cfactor(complex x, double a)
{
	complex c;
	c.re = a * x.re;
	c.im = a * x.im;
	return (c);
}

int cxtoy(complex x[], complex y[], int n)
{
	int i;
	FOR(i, n)
	y[i] = x[i];
	return (0);
}

int cmatby(complex a[], complex b[], complex c[], int n, int m, int k)
/* a[n*m], b[m*k], c[n*k]  ......  c = a*b
 */
{
	int i, j, i1;
	complex t;

	FOR(i, n)
	FOR(j, k)
	{
		for (i1 = 0, t = compl(0, 0); i1 < m; i1++)
			t = cplus(t, cby(a[i * m + i1], b[i1 * k + j]));
		c[i * k + j] = t;
	}
	return (0);
}

int cmatout(FILE *fout, complex x[], int n, int m)
{
	int i, j;
	for (i = 0, FPN(fout); i < n; i++, FPN(fout))
		FOR(j, m)
	fprintf(fout, "%7.3f%7.3f  ", x[i * m + j].re, x[i * m + j].im);
	return (0);
}

int cmatinv(complex x[], int n, int m, double space[])
{
	/* x[n*m]  ... m>=n
	 */
	int i, j, k, *irow = (int *)space;
	double xmaxsize, ee = 1e-20;
	complex xmax, t, t1;

	FOR(i, n)
	{
		xmaxsize = 0.;
		for (j = i; j < n; j++)
		{
			if (xmaxsize < csize(x[j * m + i]))
			{
				xmaxsize = csize(x[j * m + i]);
				xmax = x[j * m + i];
				irow[i] = j;
			}
		}
		if (xmaxsize < ee)
		{
			printf("\nDet goes to zero at %8d!\t\n", i + 1);
			return (-1);
		}
		if (irow[i] != i)
		{
			FOR(j, m)
			{
				t = x[i * m + j];
				x[i * m + j] = x[irow[i] * m + j];
				x[irow[i] * m + j] = t;
			}
		}
		t = cdiv(compl(1, 0), x[i * m + i]);
		FOR(j, n)
		{
			if (j == i)
				continue;
			t1 = cby(t, x[j * m + i]);
			FOR(k, m)
			x[j * m + k] = cminus(x[j * m + k], cby(t1, x[i * m + k]));
			x[j * m + i] = cfactor(t1, -1);
		}
		FOR(j, m)
		x[i * m + j] = cby(x[i * m + j], t);
		x[i * m + i] = t;
	}
	for (i = n - 1; i >= 0; i--)
	{
		if (irow[i] == i)
			continue;
		FOR(j, n)
		{
			t = x[j * m + i];
			x[j * m + i] = x[j * m + irow[i]];
			x[j * m + irow[i]] = t;
		}
	}
	return (0);
}

void balance(double mat[], int n, int *low, int *hi, double scale[])
{
	/* Balance a matrix for calculation of eigenvalues and eigenvectors
	 */
	double c, f, g, r, s;
	int i, j, k, l, done;
	/* search for rows isolating an eigenvalue and push them down */
	for (k = n - 1; k >= 0; k--)
	{
		for (j = k; j >= 0; j--)
		{
			for (i = 0; i <= k; i++)
			{
				if (i != j && fabs(mat[pos(j, i, n)]) != 0)
					break;
			}

			if (i > k)
			{
				scale[k] = j;

				if (j != k)
				{
					for (i = 0; i <= k; i++)
					{
						c = mat[pos(i, j, n)];
						mat[pos(i, j, n)] = mat[pos(i, k, n)];
						mat[pos(i, k, n)] = c;
					}

					for (i = 0; i < n; i++)
					{
						c = mat[pos(j, i, n)];
						mat[pos(j, i, n)] = mat[pos(k, i, n)];
						mat[pos(k, i, n)] = c;
					}
				}
				break;
			}
		}
		if (j < 0)
			break;
	}

	/* search for columns isolating an eigenvalue and push them left */

	for (l = 0; l <= k; l++)
	{
		for (j = l; j <= k; j++)
		{
			for (i = l; i <= k; i++)
			{
				if (i != j && fabs(mat[pos(i, j, n)]) != 0)
					break;
			}
			if (i > k)
			{
				scale[l] = j;
				if (j != l)
				{
					for (i = 0; i <= k; i++)
					{
						c = mat[pos(i, j, n)];
						mat[pos(i, j, n)] = mat[pos(i, l, n)];
						mat[pos(i, l, n)] = c;
					}

					for (i = l; i < n; i++)
					{
						c = mat[pos(j, i, n)];
						mat[pos(j, i, n)] = mat[pos(l, i, n)];
						mat[pos(l, i, n)] = c;
					}
				}

				break;
			}
		}

		if (j > k)
			break;
	}

	*hi = k;
	*low = l;

	/* balance the submatrix in rows l through k */

	for (i = l; i <= k; i++)
	{
		scale[i] = 1;
	}

	do
	{
		for (done = 1, i = l; i <= k; i++)
		{
			for (c = 0, r = 0, j = l; j <= k; j++)
			{
				if (j != i)
				{
					c += fabs(mat[pos(j, i, n)]);
					r += fabs(mat[pos(i, j, n)]);
				}
			}

			if (c != 0 && r != 0)
			{
				g = r / BASE;
				f = 1;
				s = c + r;

				while (c < g)
				{
					f *= BASE;
					c *= BASE * BASE;
				}

				g = r * BASE;

				while (c >= g)
				{
					f /= BASE;
					c /= BASE * BASE;
				}

				if ((c + r) / f < 0.95 * s)
				{
					done = 0;
					g = 1 / f;
					scale[i] *= f;

					for (j = l; j < n; j++)
					{
						mat[pos(i, j, n)] *= g;
					}

					for (j = 0; j <= k; j++)
					{
						mat[pos(j, i, n)] *= f;
					}
				}
			}
		}
	} while (!done);
}

/*
 * Transform back eigenvectors of a balanced matrix
 * into the eigenvectors of the original matrix
 */
void unbalance(int n, double vr[], double vi[], int low, int hi, double scale[])
{
	int i, j, k;
	double tmp;

	for (i = low; i <= hi; i++)
	{
		for (j = 0; j < n; j++)
		{
			vr[pos(i, j, n)] *= scale[i];
			vi[pos(i, j, n)] *= scale[i];
		}
	}

	for (i = low - 1; i >= 0; i--)
	{
		if ((k = (int)scale[i]) != i)
		{
			for (j = 0; j < n; j++)
			{
				tmp = vr[pos(i, j, n)];
				vr[pos(i, j, n)] = vr[pos(k, j, n)];
				vr[pos(k, j, n)] = tmp;

				tmp = vi[pos(i, j, n)];
				vi[pos(i, j, n)] = vi[pos(k, j, n)];
				vi[pos(k, j, n)] = tmp;
			}
		}
	}

	for (i = hi + 1; i < n; i++)
	{
		if ((k = (int)scale[i]) != i)
		{
			for (j = 0; j < n; j++)
			{
				tmp = vr[pos(i, j, n)];
				vr[pos(i, j, n)] = vr[pos(k, j, n)];
				vr[pos(k, j, n)] = tmp;

				tmp = vi[pos(i, j, n)];
				vi[pos(i, j, n)] = vi[pos(k, j, n)];
				vi[pos(k, j, n)] = tmp;
			}
		}
	}
}

/*
 * Reduce the submatrix in rows and columns low through hi of real matrix mat to
 * Hessenberg form by elementary similarity transformations
 */
void elemhess(int job, double mat[], int n, int low, int hi, double vr[],
			  double vi[], int work[])
{
	/* work[n] */
	int i, j, m;
	double x, y;

	for (m = low + 1; m < hi; m++)
	{
		for (x = 0, i = m, j = m; j <= hi; j++)
		{
			if (fabs(mat[pos(j, m - 1, n)]) > fabs(x))
			{
				x = mat[pos(j, m - 1, n)];
				i = j;
			}
		}

		if ((work[m] = i) != m)
		{
			for (j = m - 1; j < n; j++)
			{
				y = mat[pos(i, j, n)];
				mat[pos(i, j, n)] = mat[pos(m, j, n)];
				mat[pos(m, j, n)] = y;
			}

			for (j = 0; j <= hi; j++)
			{
				y = mat[pos(j, i, n)];
				mat[pos(j, i, n)] = mat[pos(j, m, n)];
				mat[pos(j, m, n)] = y;
			}
		}

		if (x != 0)
		{
			for (i = m + 1; i <= hi; i++)
			{
				if ((y = mat[pos(i, m - 1, n)]) != 0)
				{
					y = mat[pos(i, m - 1, n)] = y / x;

					for (j = m; j < n; j++)
					{
						mat[pos(i, j, n)] -= y * mat[pos(m, j, n)];
					}

					for (j = 0; j <= hi; j++)
					{
						mat[pos(j, m, n)] += y * mat[pos(j, i, n)];
					}
				}
			}
		}
	}
	if (job)
	{
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				vr[pos(i, j, n)] = 0.0;
				vi[pos(i, j, n)] = 0.0;
			}
			vr[pos(i, i, n)] = 1.0;
		}

		for (m = hi - 1; m > low; m--)
		{
			for (i = m + 1; i <= hi; i++)
			{
				vr[pos(i, m, n)] = mat[pos(i, m - 1, n)];
			}

			if ((i = work[m]) != m)
			{
				for (j = m; j <= hi; j++)
				{
					vr[pos(m, j, n)] = vr[pos(i, j, n)];
					vr[pos(i, j, n)] = 0.0;
				}
				vr[pos(i, m, n)] = 1.0;
			}
		}
	}
}

/*
 * Calculate eigenvalues and eigenvectors of a real upper Hessenberg matrix
 * Return 1 if converges successfully and 0 otherwise
 */

int realeig(int job, double mat[], int n, int low, int hi, double valr[],
			double vali[], double vr[], double vi[])
{
	complex v;
	double p = 0, q = 0, r = 0, s = 0, t, w, x, y, z = 0, ra, sa, norm, eps;
	int niter, en, i, j, k, l, m;
	double precision = pow((double)BASE, (double)(1 - DIGITS));

	eps = precision;
	for (i = 0; i < n; i++)
	{
		valr[i] = 0.0;
		vali[i] = 0.0;
	}
	/* store isolated roots and calculate norm */
	for (norm = 0, i = 0; i < n; i++)
	{
		for (j = max(0, i - 1); j < n; j++)
		{
			norm += fabs(mat[pos(i, j, n)]);
		}
		if (i < low || i > hi)
			valr[i] = mat[pos(i, i, n)];
	}
	t = 0;
	en = hi;

	while (en >= low)
	{
		niter = 0;
		for (;;)
		{

			/* look for single small subdiagonal element */

			for (l = en; l > low; l--)
			{
				s = fabs(mat[pos(l - 1, l - 1, n)]) + fabs(mat[pos(l, l, n)]);
				if (s == 0)
					s = norm;
				if (fabs(mat[pos(l, l - 1, n)]) <= eps * s)
					break;
			}

			/* form shift */

			x = mat[pos(en, en, n)];

			if (l == en)
			{ /* one root found */
				valr[en] = x + t;
				if (job)
					mat[pos(en, en, n)] = x + t;
				en--;
				break;
			}

			y = mat[pos(en - 1, en - 1, n)];
			w = mat[pos(en, en - 1, n)] * mat[pos(en - 1, en, n)];

			if (l == en - 1)
			{ /* two roots found */
				p = (y - x) / 2;
				q = p * p + w;
				z = sqrt(fabs(q));
				x += t;
				if (job)
				{
					mat[pos(en, en, n)] = x;
					mat[pos(en - 1, en - 1, n)] = y + t;
				}
				if (q < 0)
				{ /* complex pair */
					valr[en - 1] = x + p;
					vali[en - 1] = z;
					valr[en] = x + p;
					vali[en] = -z;
				}
				else
				{ /* real pair */
					z = (p < 0) ? p - z : p + z;
					valr[en - 1] = x + z;
					valr[en] = (z == 0) ? x + z : x - w / z;
					if (job)
					{
						x = mat[pos(en, en - 1, n)];
						s = fabs(x) + fabs(z);
						p = x / s;
						q = z / s;
						r = sqrt(p * p + q * q);
						p /= r;
						q /= r;
						for (j = en - 1; j < n; j++)
						{
							z = mat[pos(en - 1, j, n)];
							mat[pos(en - 1, j, n)] = q * z + p *
																 mat[pos(en, j, n)];
							mat[pos(en, j, n)] = q * mat[pos(en, j, n)] - p * z;
						}
						for (i = 0; i <= en; i++)
						{
							z = mat[pos(i, en - 1, n)];
							mat[pos(i, en - 1, n)] = q * z + p * mat[pos(i, en, n)];
							mat[pos(i, en, n)] = q * mat[pos(i, en, n)] - p * z;
						}
						for (i = low; i <= hi; i++)
						{
							z = vr[pos(i, en - 1, n)];
							vr[pos(i, en - 1, n)] = q * z + p * vr[pos(i, en, n)];
							vr[pos(i, en, n)] = q * vr[pos(i, en, n)] - p * z;
						}
					}
				}
				en -= 2;
				break;
			}
			if (niter == MAXITER)
				return (-1);
			if (niter != 0 && niter % 10 == 0)
			{
				t += x;
				for (i = low; i <= en; i++)
					mat[pos(i, i, n)] -= x;
				s = fabs(mat[pos(en, en - 1, n)]) + fabs(mat[pos(en - 1, en - 2, n)]);
				x = y = 0.75 * s;
				w = -0.4375 * s * s;
			}
			niter++;
			/* look for two consecutive small subdiagonal elements */
			for (m = en - 2; m >= l; m--)
			{
				z = mat[pos(m, m, n)];
				r = x - z;
				s = y - z;
				p = (r * s - w) / mat[pos(m + 1, m, n)] + mat[pos(m, m + 1, n)];
				q = mat[pos(m + 1, m + 1, n)] - z - r - s;
				r = mat[pos(m + 2, m + 1, n)];
				s = fabs(p) + fabs(q) + fabs(r);
				p /= s;
				q /= s;
				r /= s;
				if (m == l || fabs(mat[pos(m, m - 1, n)]) * (fabs(q) + fabs(r)) <=
								  eps * (fabs(mat[pos(m - 1, m - 1, n)]) + fabs(z) + fabs(mat[pos(m + 1, m + 1, n)])) * fabs(p))
					break;
			}
			for (i = m + 2; i <= en; i++)
				mat[pos(i, i - 2, n)] = 0;
			for (i = m + 3; i <= en; i++)
				mat[pos(i, i - 3, n)] = 0;
			/* double QR step involving rows l to en and columns m to en */
			for (k = m; k < en; k++)
			{
				if (k != m)
				{
					p = mat[pos(k, k - 1, n)];
					q = mat[pos(k + 1, k - 1, n)];
					r = (k == en - 1) ? 0 : mat[pos(k + 2, k - 1, n)];
					if ((x = fabs(p) + fabs(q) + fabs(r)) == 0)
						continue;
					p /= x;
					q /= x;
					r /= x;
				}
				s = sqrt(p * p + q * q + r * r);
				if (p < 0)
					s = -s;
				if (k != m)
				{
					mat[pos(k, k - 1, n)] = -s * x;
				}
				else if (l != m)
				{
					mat[pos(k, k - 1, n)] = -mat[pos(k, k - 1, n)];
				}
				p += s;
				x = p / s;
				y = q / s;
				z = r / s;
				q /= p;
				r /= p;
				/* row modification */
				for (j = k; j <= (!job ? en : n - 1); j++)
				{
					p = mat[pos(k, j, n)] + q * mat[pos(k + 1, j, n)];
					if (k != en - 1)
					{
						p += r * mat[pos(k + 2, j, n)];
						mat[pos(k + 2, j, n)] -= p * z;
					}
					mat[pos(k + 1, j, n)] -= p * y;
					mat[pos(k, j, n)] -= p * x;
				}
				j = min(en, k + 3);
				/* column modification */
				for (i = (!job ? l : 0); i <= j; i++)
				{
					p = x * mat[pos(i, k, n)] + y * mat[pos(i, k + 1, n)];
					if (k != en - 1)
					{
						p += z * mat[pos(i, k + 2, n)];
						mat[pos(i, k + 2, n)] -= p * r;
					}
					mat[pos(i, k + 1, n)] -= p * q;
					mat[pos(i, k, n)] -= p;
				}
				if (job)
				{ /* accumulate transformations */
					for (i = low; i <= hi; i++)
					{
						p = x * vr[pos(i, k, n)] + y * vr[pos(i, k + 1, n)];
						if (k != en - 1)
						{
							p += z * vr[pos(i, k + 2, n)];
							vr[pos(i, k + 2, n)] -= p * r;
						}
						vr[pos(i, k + 1, n)] -= p * q;
						vr[pos(i, k, n)] -= p;
					}
				}
			}
		}
	}

	if (!job)
		return (0);
	if (norm != 0)
	{
		/* back substitute to find vectors of upper triangular form */
		for (en = n - 1; en >= 0; en--)
		{
			p = valr[en];
			if ((q = vali[en]) < 0)
			{ /* complex vector */
				m = en - 1;
				if (fabs(mat[pos(en, en - 1, n)]) > fabs(mat[pos(en - 1, en, n)]))
				{
					mat[pos(en - 1, en - 1, n)] = q / mat[pos(en, en - 1, n)];
					mat[pos(en - 1, en, n)] = (p - mat[pos(en, en, n)]) /
											  mat[pos(en, en - 1, n)];
				}
				else
				{
					v = cdiv(compl(0.0, -mat[pos(en - 1, en, n)]),
							 compl(mat[pos(en - 1, en - 1, n)] - p, q));
					mat[pos(en - 1, en - 1, n)] = v.re;
					mat[pos(en - 1, en, n)] = v.im;
				}
				mat[pos(en, en - 1, n)] = 0;
				mat[pos(en, en, n)] = 1;
				for (i = en - 2; i >= 0; i--)
				{
					w = mat[pos(i, i, n)] - p;
					ra = 0;
					sa = mat[pos(i, en, n)];
					for (j = m; j < en; j++)
					{
						ra += mat[pos(i, j, n)] * mat[pos(j, en - 1, n)];
						sa += mat[pos(i, j, n)] * mat[pos(j, en, n)];
					}
					if (vali[i] < 0)
					{
						z = w;
						r = ra;
						s = sa;
					}
					else
					{
						m = i;
						if (vali[i] == 0)
						{
							v = cdiv(compl(-ra, -sa), compl(w, q));
							mat[pos(i, en - 1, n)] = v.re;
							mat[pos(i, en, n)] = v.im;
						}
						else
						{ /* solve complex equations */
							x = mat[pos(i, i + 1, n)];
							y = mat[pos(i + 1, i, n)];
							v.re = (valr[i] - p) * (valr[i] - p) + vali[i] * vali[i] - q * q;
							v.im = (valr[i] - p) * 2 * q;
							if ((fabs(v.re) + fabs(v.im)) == 0)
							{
								v.re = eps * norm * (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(z));
							}
							v = cdiv(compl(x * r - z * ra + q * sa, x * s - z * sa - q * ra), v);
							mat[pos(i, en - 1, n)] = v.re;
							mat[pos(i, en, n)] = v.im;
							if (fabs(x) > fabs(z) + fabs(q))
							{
								mat[pos(i + 1, en - 1, n)] =
									(-ra - w * mat[pos(i, en - 1, n)] +
									 q * mat[pos(i, en, n)]) /
									x;
								mat[pos(i + 1, en, n)] = (-sa - w * mat[pos(i, en, n)] -
														  q * mat[pos(i, en - 1, n)]) /
														 x;
							}
							else
							{
								v = cdiv(compl(-r - y * mat[pos(i, en - 1, n)],
											   -s - y * mat[pos(i, en, n)]),
										 compl(z, q));
								mat[pos(i + 1, en - 1, n)] = v.re;
								mat[pos(i + 1, en, n)] = v.im;
							}
						}
					}
				}
			}
			else if (q == 0)
			{ /* real vector */
				m = en;
				mat[pos(en, en, n)] = 1;
				for (i = en - 1; i >= 0; i--)
				{
					w = mat[pos(i, i, n)] - p;
					r = mat[pos(i, en, n)];
					for (j = m; j < en; j++)
					{
						r += mat[pos(i, j, n)] * mat[pos(j, en, n)];
					}
					if (vali[i] < 0)
					{
						z = w;
						s = r;
					}
					else
					{
						m = i;
						if (vali[i] == 0)
						{
							if ((t = w) == 0)
								t = eps * norm;
							mat[pos(i, en, n)] = -r / t;
						}
						else
						{ /* solve real equations */
							x = mat[pos(i, i + 1, n)];
							y = mat[pos(i + 1, i, n)];
							q = (valr[i] - p) * (valr[i] - p) + vali[i] * vali[i];
							t = (x * s - z * r) / q;
							mat[pos(i, en, n)] = t;
							if (fabs(x) <= fabs(z))
							{
								mat[pos(i + 1, en, n)] = (-s - y * t) / z;
							}
							else
							{
								mat[pos(i + 1, en, n)] = (-r - w * t) / x;
							}
						}
					}
				}
			}
		}
		/* vectors of isolated roots */
		for (i = 0; i < n; i++)
		{
			if (i < low || i > hi)
			{
				for (j = i; j < n; j++)
				{
					vr[pos(i, j, n)] = mat[pos(i, j, n)];
				}
			}
		}
		/* multiply by transformation matrix */

		for (j = n - 1; j >= low; j--)
		{
			m = min(j, hi);
			for (i = low; i <= hi; i++)
			{
				for (z = 0, k = low; k <= m; k++)
				{
					z += vr[pos(i, k, n)] * mat[pos(k, j, n)];
				}
				vr[pos(i, j, n)] = z;
			}
		}
	}
	/* rearrange complex eigenvectors */
	for (j = 0; j < n; j++)
	{
		if (vali[j] != 0)
		{
			for (i = 0; i < n; i++)
			{
				vi[pos(i, j, n)] = vr[pos(i, j + 1, n)];
				vr[pos(i, j + 1, n)] = vr[pos(i, j, n)];
				vi[pos(i, j + 1, n)] = -vi[pos(i, j, n)];
			}
			j++;
		}
	}
	return (0);
}

int matinv(double x[], int n, int m, double space[])
{
	/* x[n*m]  ... m>=n
	 */
	register int i, j, k;
	int *irow = (int *)space;
	double ee = 1.0e-20, t, t1, xmax;
	double det = 1.0;

	FOR(i, n)
	{
		xmax = 0.;
		for (j = i; j < n; j++)
		{
			if (xmax < fabs(x[j * m + i]))
			{
				xmax = fabs(x[j * m + i]);
				irow[i] = j;
			}
		}
		det *= xmax;
		if (xmax < ee)
		{
			printf("\nDet becomes zero at %3d!\t\n", i + 1);
			return (-1);
		}
		if (irow[i] != i)
		{
			FOR(j, m)
			{
				t = x[i * m + j];
				x[i * m + j] = x[irow[i] * m + j];
				x[irow[i] * m + j] = t;
			}
		}
		t = 1. / x[i * m + i];
		FOR(j, n)
		{
			if (j == i)
				continue;
			t1 = t * x[j * m + i];
			FOR(k, m)
			x[j * m + k] -= t1 * x[i * m + k];
			x[j * m + i] = -t1;
		}
		FOR(j, m)
		x[i * m + j] *= t;
		x[i * m + i] = t;
	} /* i  */
	for (i = n - 1; i >= 0; i--)
	{
		if (irow[i] == i)
			continue;
		FOR(j, n)
		{
			t = x[j * m + i];
			x[j * m + i] = x[j * m + irow[i]];
			x[j * m + irow[i]] = t;
		}
	}
	return (0);
}

void make_transition_prob_matrices(double t[3], unsigned long int treeNum)
{
	int i, j, k, v, n;
	double EXPOS[4], T;
	// double sum;
	double *PMAT_ptr = PMAT;
	// int pointerReplace = 0;

	for (v = 0; v < 3; v++)
	{
		// printf("branch %i\n",v);
		for (n = 0; n < NUMCAT; n++)
		{
			// printf("Category %i\n",n);
			T = statevector[treeNum][n] * t[v];
			for (k = 0; k < 4; k++)
				EXPOS[k] = exp(T * RRVAL[treeNum][k]);
			for (i = 0; i < 4; i++)
			{
				// sum = 0.0;
				for (j = 0; j < 4; j++)
				{
					// PMAT[pointerReplace] = 0.0;
					*PMAT_ptr = 0.0;
					for (k = 0; k < 4; k++)
						*PMAT_ptr += RRVEC[treeNum][k * 4 + j] * LRVEC[treeNum][i * 4 + k] * EXPOS[k];
					// PMAT[pointerReplace] += RRVEC[treeNum][k * 4 + j]*LRVEC[treeNum][i * 4 + k]*EXPOS[k];
					// printf("%.16f:",PMAT[pointerReplace]);
					// printf("%.16f:", *PMAT_ptr);
					// sum += PMAT[v][n][i][j];
					// if(PMAT[pointerReplace] <= 0.0)
					if (*PMAT_ptr <= 0.0)
					{
						// PMAT[pointerReplace] = 0.00000001;
						*PMAT_ptr = 0.00000001;
						// printf("%.16f:", *PMAT_ptr);
						// printf("%.16f:",PMAT[pointerReplace]);
					}
					// PMAT[pointerReplace] = log(PMAT[pointerReplace]);
					*PMAT_ptr = log(*PMAT_ptr);
					// printf("%.16f\t", *PMAT_ptr);
					// printf("%.16f\t",PMAT[pointerReplace]);
					// pointerReplace++;
					PMAT_ptr++;
				}
				// if (sum > 1.00001 || sum < 0.99999)
				//	printf("Sum: %.16f\n", sum);
			}
		}
	}
}

// double mydistance(double v1[], double v2[], int start1, int start2, int n)
//{
//	int i;
//	double sum=0.0;
//
//	for (i=0; i<n; i++)
//		sum+=(v1[i+start1]-v2[i+start2])*(v1[i+start1]-v2[i+start2]);
//	return sqrt(sum);
// }

// log sum exp trick sum(exp(X)-maxX) + maxX
// Calculates the log sum exp of log like + log prob for the tree calculations
// This was implemented after scaling at base/pos was underflowing
// always used on dim 4 vectors
static inline double logSumExp(double X[4])
{
	int i;
	double maxX = X[0];
	double sumX = 0.0;

	// Find max
	for (i = 1; i < 4; i++)
	{
		if (X[i] > maxX)
			maxX = X[i];
	}

	//printf("\t\t\tMax: %.16f\n", maxX);

	// sum(exp(x) - maxX)
	for (i = 0; i < 4; i++)
	{
		//printf("\t\t\t\t\tX[%d]: %.16f \t %.16f\t %.16f\n", i, X[i], exp(X[i] - maxX), X[i] - maxX);
		sumX += exp(X[i] - maxX);
		//printf("\t\t\tsumX: %.16f\n", sumX);
	}

	//printf("\t\treturn: %.16f\n", log(sumX) + maxX);

	return (log(sumX) + maxX);
}

// transtion probability matrix has already been diagonalized
// this function assumes that the only thing that changes between function calls is the branch lengths or the query sequence and its placement
// node1 and node2 are the nodes around the edge to which the sequence has been assigned
// seq is the identifer of the sequence
// pi are the nucletoide frequencies - this could be made a global to avoid passing them around
// Alignment data should be stored in DATA with -1 indicating missing data
// upper and lower bounds are parameters 7-10 when 0-counting
double getlike_gamma_root_in_trifurcation(double times[3], double parameters[7])
{
	int i, j, k, b, c, v, po, node, seq, treeNum;
	double Like, t[3], A[4], B[4], C[4];
	double *PMAT_ptr_0, *PMAT_ptr_1, *PMAT_ptr_2, *FRACLIKE_ptr;
	// diagonalizaiton has previously been done

	if (onDindic == 1)
	{ // Learn more about what these do, seem to be hitting these values
		if (times[2] < parameters[3] || times[2] > parameters[4])
		{
			return 1000000000.0; // If an actual likelihood is smaler than this we are screwed
		}
	}
	else
	{
		if ((times[2] < parameters[5] || times[2] > parameters[6]) || (times[1] < parameters[3] || times[1] > parameters[4]))
		{
			return 1000000000.0; // If an actual likelihood is smaler than this we are screwed
		}
	}

	node = (int)parameters[2];
	seq = (int)parameters[0];
	treeNum = (int)parameters[1];

	t[1] = times[2];										// length from node to position where query joins
	t[0] = (times[2] + nodeages[treeNum][node]) * times[1]; // length from age of query node to position where query joins
	t[2] = bls[treeNum][node] - times[2];					// length from position where query joins to parent node

	make_transition_prob_matrices(t, treeNum);

	Like = 0.0;
	//FRACLIKE_ptr = &FRACLIKE[treeNum][node * numbases[treeNum] * NUMCAT * 8 + startpos[seq] * NUMCAT * 8];
	// printf("%d\t%d\n", numbases[treeNum], node * numbases[treeNum] * NUMCAT * 8 + startpos[seq] * NUMCAT * 8);
	// FRACLIKE[treeNum][k * numbase * NUMCAT * 8 + j * NUMCAT * 8 + i * 8 + v] = a;
	unsigned long long int size = ((unsigned long long)node) * ((unsigned long long)numbases[treeNum]) * ((unsigned long long)NUMCAT) * 8ULL
                           + ((unsigned long long)startpos[seq]) * ((unsigned long long)NUMCAT) * 8ULL;
	//printf("FRACLIKE index 4: %d vs %llu\n", node * numbases[treeNum] * NUMCAT * 8 + startpos[seq] * NUMCAT * 8, size);
	FRACLIKE_ptr = &FRACLIKE[treeNum][size];

	// printf("Sequence %d\n", seq);

	for (i = startpos[seq]; i < readlength[seq] + startpos[seq]; i++)
	{
		po = i - startpos[seq]; // i keeps track of the positon in the ref sequences while po is the positoon in the read
		b = QUERYDATA[seq][po];

		// printf("\tBase %d\n", b);

		if (b == -1)
		{
			FRACLIKE_ptr += 32; // NUMCAT (4) * 8;
		}
		// if (b!=-1){
		else
		{
			// PMAT_ptr_0 = PMAT; 		//0*NUMCAT*4*4
			PMAT_ptr_1 = &PMAT[64];	 // 1*NUMCAT(4)*4*4
			PMAT_ptr_2 = &PMAT[128]; // 2*NUMCAT(4)*4*4
			for (j = 0; j < NUMCAT; j++)
			{
				// printf("\t\tj: %d\n", j);
				// PMAT_ptr_0 = &PMAT[j * 4 * 4 + b * 4 + b];	// when no error
				PMAT_ptr_0 = &PMAT[j * 4 * 4];
				for (k = 0; k < 4; k++)
				{
					// printf("\t\t\tk: %d\n", k);
					for (v=0; v<4; v++)
					{
						//A[v] = pi[treeNum][v] + PMAT[0][j][v][k] + readlike[seq][po][v];
						//A[v] = pi[treeNum][v] + *PMAT_ptr_0 + readlike[seq][po][v];
						A[v] = pi[treeNum][v] + *(PMAT_ptr_0 + v * 4 + k) + readlike[seq][po][v];
					//	//printf("\t\t\t\tpi[%d] is %.16f anf PMAT[0] is %.16f anf reaflike is %.16f\n", v, pi[treeNum][v], PMAT[0][j][v][k], readlike[seq][po][v]);
					}
					B[k] = logSumExp(A);
					//B[k] = pi[treeNum][b] + *PMAT_ptr_0;
					//// printf("\t\t\t\t\tB:%lf\tpi:%lf\tPMAT[0]:%lf\n", B[k], pi[treeNum][b], PMAT_ptr_0);
					//PMAT_ptr_0++;
					if (node >= numseq[treeNum])
					{ // If not leaf node. Assumes t,c,g,t)leaf nodes are numbered from 0 to numseq-1
						for (v = 0; v < 4; v++)
						{
							// printf("\t\t\t\t\tPMAT[2] is %.16f and FRACLIKE[%d] is %.16f\n", PMAT[2][j][k][v], v, FRACLIKE[treeNum][i][node][j][v]);
							// A[v] = *PMAT_ptr_2 + FRACLIKE[treeNum][i][node][j][v];
							A[v] = *PMAT_ptr_2 + *FRACLIKE_ptr;
							// printf("%d: %.16f\n", seq * numbases[treeNum] * NUMCAT * 8 + i * NUMCAT * 8 + j * 8 + v, *FRACLIKE_ptr);
							PMAT_ptr_2++;
							FRACLIKE_ptr++;
						}
						B[k] += logSumExp(A);
						// printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
						for (v = 0; v < 4; v++)
						{
							// printf("\t\t\t\t\tPMAT[1] is %.16f and FRACLIKE[%d] is %.16f\n", (PMAT_ptr_1, v+4, FRACLIKE[treeNum][i][node][j][v+4]);
							// A[v] = *PMAT_ptr_1 + FRACLIKE[treeNum][i][node][j][v+4];
							A[v] = *PMAT_ptr_1 + *FRACLIKE_ptr;
							// printf("%d: %.16f\n", seq * numbases[treeNum] * NUMCAT * 8 + i * NUMCAT * 8 + j * 8 + v, *FRACLIKE_ptr);
							PMAT_ptr_1++;
							FRACLIKE_ptr++;
						}
						B[k] += logSumExp(A);
						FRACLIKE_ptr -= 8;
						// printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
					}
					else
					{ // If leaf node
						if ((c = DATA[treeNum][node][i]) > -1)
						{
							B[k] += *(PMAT_ptr_1 + c);
							PMAT_ptr_1 += 4;
						}
						for (v = 0; v < 4; v++) // add position here into fraclike
						{
							// printf("\t\t\t\t\tPMAT[2] is %.16f and FRACLIKE[%d] is %.16f\n", *PMAT_ptr_2, v, FRACLIKE[treeNum][i][node][j][v]);
							// A[v] = *PMAT_ptr_2 + FRACLIKE[treeNum][i][node][j][v];
							A[v] = *PMAT_ptr_2 + *FRACLIKE_ptr;
							// printf("%d: %.16f\n", seq * numbases[treeNum] * NUMCAT * 8 + i * NUMCAT * 8 + j * 8 + v, *FRACLIKE_ptr);
							PMAT_ptr_2++;
							FRACLIKE_ptr++;
						}
						B[k] += logSumExp(A);
						FRACLIKE_ptr -= 4;
						// printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
					}
					// FRACLIKE_ptr -= 8;
					//  printf("\t\t\t\t\tB[%d]: %lf\n", k, B[k]);
				}
				FRACLIKE_ptr += 8;
				// for(int l = 0; l < 4; l++)
				// {
				// 	printf("\t\t\t\tB[%d]: %lf\n", l, B[l]);
				// }
				C[j] = logSumExp(B);
				// printf("\t\t\t\tC:%lf\n", logSumExp(B));
			}
			Like += logSumExp(C);
			// printf("\t\t\t\tLike:%lf\n", logSumExp(C));
		}
	}

	// printf("\tLikelihood: %.16f\n", Like);
	// printf("\t\t%d,%.16f,%.16f,%.16f,%.16f,%d,%d,%d\n", seq, (1.0-times[1])*(nodeages[treeNum][seq]+times[2]), times[1], times[2], Like, node * numbases[treeNum] * NUMCAT * 8 + startpos[seq] * NUMCAT * 8, startpos[seq],node * numbases[treeNum] * NUMCAT * 8 + numbases[treeNum] * NUMCAT * 8);

	return -Like; // Notice: a scaling factor of NUMCAT^(number of sites) is missing
}

double getlike_gamma_root_in_trifucation_single_read_brent(double times[], void *extra_data)
{
	int i, j, k, b, c, v, po, node, seq, treeNum;
	double Like, t[3], A[4], B[4], C[4];
	double *PMAT_ptr_0, *PMAT_ptr_1, *PMAT_ptr_2, *FRACLIKE_ptr;
	// diagonalizaiton has previously been done

	seq = *((int *)extra_data);
	node = assignments[seq];
	treeNum = treeAssign[seq];
	// Check bounds and correct if parameter outside of bounds

	//printf("times[0]: %.16f\t times[1]: %.16f\t", times[0], times[1]);

	if (times[0] < 3e-8)
	{
		times[0] = 3e-8;
		//printf("Likelihood: %.16f\n", 1000000000.0 - times[0]);
		//return (1000000000.0 - times[0]);
	}
	else if (times[0] > 1.0 - 3e-8)
	{
		times[0] = 1.0 - 3e-8;
		//printf("Likelihood: %.16f\n", 1000000000.0 + times[0]);
		//return (1000000000.0 + times[0]);
	}

	if (times[1] < 3e-8)
	{
		times[1] = 3e-8;
		//printf("Likelihood: %.16f\n", 1000000000.0 - times[1]);
		//return (1000000000.0 - times[1]);
	}
	else if (times[1] > bls[treeNum][node] - 3e-8)
	{
		times[1] = bls[treeNum][node] - 3e-8;
		//printf("Likelihood: %.16f\n", 1000000000.0 + times[1]);
		//return (1000000000.0 + times[1]);
	}

	//printf("times[0]: %.16f\t times[1]: %.16f\t", times[0], times[1]);

	t[1] = times[1];										// length from node to position where query joins
	t[0] = (times[1] + nodeages[treeNum][node]) * times[0]; // length from age of query node to position where query joins
	t[2] = bls[treeNum][node] - times[1];					// length from position where query joins to parent node

	make_transition_prob_matrices(t, treeNum);

	Like = 0.0;
	//FRACLIKE_ptr = &FRACLIKE[treeNum][node * numbases[treeNum] * NUMCAT * 8 + startpos[seq] * NUMCAT * 8];
	// printf("%d\t%d\n", numbases[treeNum], node * numbases[treeNum] * NUMCAT * 8 + startpos[seq] * NUMCAT * 8);
	unsigned long long int size = ((unsigned long long)node) * ((unsigned long long)numbases[treeNum]) * ((unsigned long long)NUMCAT) * 8ULL
	                            + ((unsigned long long)startpos[seq]) * ((unsigned long long)NUMCAT) * 8ULL;
	// printf("FRACLIKE index 5: %d vs %llu\n", node * numbases[treeNum] * NUMCAT * 8 + startpos[seq] * NUMCAT * 8, size);
	FRACLIKE_ptr = &FRACLIKE[treeNum][size];

	// printf("Sequence %d\n", seq);

	for (i = startpos[seq]; i < readlength[seq] + startpos[seq]; i++)
	{
		po = i - startpos[seq]; // i keeps track of the positon in the ref sequences while po is the positoon in the read
		b = QUERYDATA[seq][po];

		// printf("\tBase %d\n", b);

		if (b == -1)
		{
			FRACLIKE_ptr += 32; // NUMCAT (4) * 8;
		}
		// if (b!=-1){
		else
		{
			// PMAT_ptr_0 = PMAT; 		//0*NUMCAT*4*4
			PMAT_ptr_1 = &PMAT[64];	 // 1*NUMCAT(4)*4*4
			PMAT_ptr_2 = &PMAT[128]; // 2*NUMCAT(4)*4*4
			for (j = 0; j < NUMCAT; j++)
			{
				// printf("\t\tj: %d\n", j);
				// PMAT_ptr_0 = &PMAT[j * 4 * 4 + b * 4];
				PMAT_ptr_0 = &PMAT[j * 4 * 4];
				for (k = 0; k < 4; k++)
				{
					// printf("\t\t\tk: %d\n", k);
					for (v=0; v<4; v++)
					{
						// A[v] = pi[treeNum][v] + PMAT[0][j][v][k] + readlike[seq][po][v];
						A[v] = pi[treeNum][v] + *(PMAT_ptr_0 + v * 4 + k) + readlike[seq][po][v];
					//	//printf("\t\t\t\tpi[%d] is %.16f anf PMAT[0] is %.16f anf reaflike is %.16f\n", v, pi[treeNum][v], PMAT[0][j][v][k], readlike[seq][po][v]);
					}
					B[k] = logSumExp(A);
					// B[k] = pi[treeNum][b] + *PMAT_ptr_0;
					// // printf("\t\t\t\t\tB:%lf\tpi:%lf\tPMAT[0]:%lf\n", B[k], pi[treeNum][b], PMAT_ptr_0);
					// PMAT_ptr_0++;
					if (node >= numseq[treeNum])
					{ // If not leaf node. Assumes t,c,g,t)leaf nodes are numbered from 0 to numseq-1
						for (v = 0; v < 4; v++)
						{
							// printf("\t\t\t\t\tPMAT[2] is %.16f and FRACLIKE[%d] is %.16f\n", PMAT[2][j][k][v], v, FRACLIKE[treeNum][i][node][j][v]);
							// A[v] = *PMAT_ptr_2 + FRACLIKE[treeNum][i][node][j][v];
							A[v] = *PMAT_ptr_2 + *FRACLIKE_ptr;
							// printf("%d: %.16f\n", seq * numbases[treeNum] * NUMCAT * 8 + i * NUMCAT * 8 + j * 8 + v, *FRACLIKE_ptr);
							PMAT_ptr_2++;
							FRACLIKE_ptr++;
						}
						B[k] += logSumExp(A);
						// printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
						for (v = 0; v < 4; v++)
						{
							// printf("\t\t\t\t\tPMAT[1] is %.16f and FRACLIKE[%d] is %.16f\n", (PMAT_ptr_1, v+4, FRACLIKE[treeNum][i][node][j][v+4]);
							// A[v] = *PMAT_ptr_1 + FRACLIKE[treeNum][i][node][j][v+4];
							A[v] = *PMAT_ptr_1 + *FRACLIKE_ptr;
							// printf("%d: %.16f\n", seq * numbases[treeNum] * NUMCAT * 8 + i * NUMCAT * 8 + j * 8 + v, *FRACLIKE_ptr);
							PMAT_ptr_1++;
							FRACLIKE_ptr++;
						}
						B[k] += logSumExp(A);
						FRACLIKE_ptr -= 8;
						// printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
					}
					else
					{ // If leaf node
						if ((c = DATA[treeNum][node][i]) > -1)
						{
							B[k] += *(PMAT_ptr_1 + c);
							PMAT_ptr_1 += 4;
						}
						for (v = 0; v < 4; v++) // add position here into fraclike
						{
							// printf("\t\t\t\t\tPMAT[2] is %.16f and FRACLIKE[%d] is %.16f\n", *PMAT_ptr_2, v, FRACLIKE[treeNum][i][node][j][v]);
							// A[v] = *PMAT_ptr_2 + FRACLIKE[treeNum][i][node][j][v];
							A[v] = *PMAT_ptr_2 + *FRACLIKE_ptr;
							// printf("%d: %.16f\n", seq * numbases[treeNum] * NUMCAT * 8 + i * NUMCAT * 8 + j * 8 + v, *FRACLIKE_ptr);
							PMAT_ptr_2++;
							FRACLIKE_ptr++;
						}
						B[k] += logSumExp(A);
						FRACLIKE_ptr -= 4;
						// printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
					}
					// FRACLIKE_ptr -= 8;
					//  printf("\t\t\t\t\tB[%d]: %lf\n", k, B[k]);
				}
				FRACLIKE_ptr += 8;
				// for(int l = 0; l < 4; l++)
				// {
				// 	printf("\t\t\t\tB[%d]: %lf\n", l, B[l]);
				// }
				C[j] = logSumExp(B);
				// printf("\t\t\t\tC:%lf\n", logSumExp(B));
			}
			Like += logSumExp(C);
			// printf("\t\t\t\tLike:%lf\n", logSumExp(C));
		}
	}

	//printf("Likelihood: %.16f\n", Like);
	// printf("\t\t%d,%.16f,%.16f,%.16f,%.16f,%d,%d,%d\n", seq, (1.0-times[1])*(nodeages[treeNum][seq]+times[2]), times[1], times[2], Like, node * numbases[treeNum] * NUMCAT * 8 + startpos[seq] * NUMCAT * 8, startpos[seq],node * numbases[treeNum] * NUMCAT * 8 + numbases[treeNum] * NUMCAT * 8);

	return -Like; // Notice: a scaling factor of NUMCAT^(number of sites) is missing
}

double getlike_gamma_root_in_trifucation_single_read_brent_reassign(double times[], void *extra_data)
{
	int i, j, k, b, c, v, po, node, seq, treeNum;
	double Like, t[3], A[4], B[4], C[4];
	double *PMAT_ptr_0, *PMAT_ptr_1, *PMAT_ptr_2, *FRACLIKE_ptr;
	// diagonalizaiton has previously been done

	seq = *((int *)extra_data);
	treeNum = *((int *)extra_data + 1);
	node = *((int *)extra_data + 2);

	//printf("\tseq: %d\t treeNum: %d\t node: %d\n", seq, treeNum, node);

	if (times[0] < 3e-8)
	{
		times[0] = 3e-8;
		//printf("Likelihood: %.16f\n", 1000000000.0 - times[0]);
		//return (1000000000.0 - times[0]);
	}
	else if (times[0] > 1.0 - 3e-8)
	{
		times[0] = 1.0 - 3e-8;
		//printf("Likelihood: %.16f\n", 1000000000.0 + times[0]);
		//return (1000000000.0 + times[0]);
	}

	if (times[1] < 3e-8)
	{
		times[1] = 3e-8;
		//printf("Likelihood: %.16f\n", 1000000000.0 - times[1]);
		//return (1000000000.0 - times[1]);
	}
	else if (times[1] > bls[treeNum][node] - 3e-8)
	{
		times[1] = bls[treeNum][node] - 3e-8;
		//printf("Likelihood: %.16f\n", 1000000000.0 + times[1]);
		//return (1000000000.0 + times[1]);
	}

	//printf("\ttimes[0]: %.16f\ttimes[1]: %.16f\n", times[0], times[1]);

	t[1] = times[1];										// length from node to position where query joins
	t[0] = (times[1] + nodeages[treeNum][node]) * times[0]; // length from age of query node to position where query joins
	t[2] = bls[treeNum][node] - times[1];					// length from position where query joins to parent node

	//printf("\tt[0]: %.16f\tt[1]: %.16f\tt[2]: %.16f\n", t[0], t[1], t[2]);

	make_transition_prob_matrices(t, treeNum);

	Like = 0.0;
	// fraclike pointer starts where the read starts in the alignment at the right node
	//FRACLIKE_ptr = &FRACLIKE[treeNum][node * numbases[treeNum] * NUMCAT * 8 + startposTemp[treeNum][seq] * NUMCAT * 8];

	unsigned long long int size = ((unsigned long long)node) * ((unsigned long long)numbases[treeNum]) * ((unsigned long long)NUMCAT) * 8ULL
	                            + ((unsigned long long)startposTemp[treeNum][seq]) * ((unsigned long long)NUMCAT) * 8ULL;

	// printf("FRACLIKE index 6: %d vs %llu\n", node * numbases[treeNum] * NUMCAT * 8 + startposTemp[treeNum][seq] * NUMCAT * 8, size);
	FRACLIKE_ptr = &FRACLIKE[treeNum][size];

	//for(k = 0; k < 4; k++)
	//{
	//	for(v = 0; v < 4; v++)
	//	{
	//		printf("PMAT[0][0][%d][%d]: %.16f\t", k, v, PMAT[k * 4 + v]);
	//	}
	//	printf("\n");
	//}

	//exit(0);

	for (i = startposTemp[treeNum][seq]; i < readLengthTemp[treeNum][seq] + startposTemp[treeNum][seq]; i++)
	{
		po = i - startposTemp[treeNum][seq]; // i keeps track of the positon in the ref sequences while po is the positoon in the read
		b = readsTreeSorted[treeNum][seq][po];

		//printf("\tBase %d\n", b);

		if (b == -1)
		{
			FRACLIKE_ptr += 32; // NUMCAT(4) * 8;
		}
		else
		{
			// if (b!=-1){
			// PMAT_ptr_0 = PMAT; 		//0*NUMCAT*4*4
			PMAT_ptr_1 = &PMAT[64];	 // 1*NUMCAT(4)*4*4
			PMAT_ptr_2 = &PMAT[128]; // 2*NUMCAT(4)*4*4
			for (j = 0; j < NUMCAT; j++)
			{
				// printf("\t\tj: %d\n", j);
				//PMAT_ptr_0 = &PMAT[j * 4 * 4 + b * 4]; // when no error
				PMAT_ptr_0 = &PMAT[j * 4 * 4];
				for (k = 0; k < 4; k++)
				{
					// printf("\t\t\tk: %d\n", k);
					for (v=0; v<4; v++)
					{
					     //A[v] = pi[treeNum][v] + PMAT[0][j][v][k] + readLikeTemp[treeNum][seq][po][v];
					     A[v] = pi[treeNum][v] + *(PMAT_ptr_0 + v * 4 + k) + readLikeTemp[treeNum][seq][po][v];
					     //printf("PMAT[0][%d][%d][%d]: %.16f\n", j, v, k, *(PMAT_ptr_0 + v * 4 + k));
					//	//printf("\t\t\t\tpi[%d] is %.16f anf PMAT[0] is %.16f anf reaflike is %.16f\n", v, pi[treeNum][v], PMAT[0][j][v][k], readLikeTemp[treeNum][seq][po][v]);
					}
					B[k] = logSumExp(A);
					//B[k] = pi[treeNum][b] + *PMAT_ptr_0;
					//printf("\t\t\t\t\tB:%lf\tpi:%lf\tPMAT[0]:%lf\n", B[k], pi[treeNum][b], *PMAT_ptr_0);
					// printf("\t\t%.16f\n", B[k]);
					//printf("PMAT[0][%d][%d][%d]: %.16f\n", j, b, k, *PMAT_ptr_0);
					//PMAT_ptr_0++;
					if (node >= numseq[treeNum])
					{ // If not leaf node. Assumes t,c,g,t)leaf nodes are numbered from 0 to numseq-1
						for (v = 0; v < 4; v++)
						{
							// printf("\t\t\t\t\tPMAT[2] is %.16f and FRACLIKE[%d] is %.16f\n", *PMAT_ptr_2, v, FRACLIKE[treeNum][i][node][j][v]);
							// A[v] = *PMAT_ptr_2 + FRACLIKE[treeNum][i][node][j][v];
							// printf("\t\tFRACLIKE[%d][%d][%d][%d][%d]:%.16f\t%d\n", treeNum, i, node, j, v, *FRACLIKE_ptr, ptr_pos);
							A[v] = *PMAT_ptr_2 + *FRACLIKE_ptr;
							PMAT_ptr_2++;
							FRACLIKE_ptr++;
							// ptr_pos++;
						}
						B[k] += logSumExp(A);
						//printf("\t\t\t\t\tB:%lf %lf\n", logSumExp(A), B[k]);
						for (v = 0; v < 4; v++)
						{
							// printf("\t\t\t\t\tPMAT[1] is %.16f and FRACLIKE[%d] is %.16f\n", *PMAT_ptr_1, v+4, FRACLIKE[treeNum][i][node][j][v+4]);
							// printf("\t\t\t\t\tPMAT[1][%d][%d][%d] is %.16f\n", j, k, v, *PMAT_ptr_1);
							// A[v] = *PMAT_ptr_1 + FRACLIKE[treeNum][i][node][j][v+4];
							// printf("\t\tFRACLIKE[%d][%d][%d][%d][%d]:%.16f\t%d\n", treeNum, i, node, j, v, *FRACLIKE_ptr, ptr_pos);
							A[v] = *PMAT_ptr_1 + *FRACLIKE_ptr;
							PMAT_ptr_1++;
							FRACLIKE_ptr++;
							// ptr_pos++;
						}
						B[k] += logSumExp(A);
						FRACLIKE_ptr -= 8;
						//printf("\t\t\t\t\tB:%lf %lf\n", logSumExp(A), B[k]);
						// printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
						//printf("\t\t%.16f\n", B[k]);
					}
					else
					{ // If leaf node
						if ((c = DATA[treeNum][node][i]) > -1)
						{
							B[k] += *(PMAT_ptr_1 + c);
							// printf("%.16f\t%d %d\n", *(PMAT_ptr_1 + c), c, 5-c);
							PMAT_ptr_1 += 4;
						}
						for (v = 0; v < 4; v++) // add position here into fraclike
						{
							// printf("\t\t\t\t\tPMAT[2] is %.16f and FRACLIKE[%d] is %.16f\n", *PMAT_ptr_2, v, FRACLIKE[treeNum][i][node][j][v]);
							// A[v] = *PMAT_ptr_2 + FRACLIKE[treeNum][i][node][j][v];
							// printf("\t\tFRACLIKE[%d][%d][%d][%d][%d]:%.16f\t%d\n", treeNum, i, node, j, v, *FRACLIKE_ptr, ptr_pos);
							A[v] = *PMAT_ptr_2 + *FRACLIKE_ptr;
							PMAT_ptr_2++;
							FRACLIKE_ptr++;
							// ptr_pos++;
						}
						B[k] += logSumExp(A);
						//printf("\t\t\t\t\tB:%lf %lf\n", logSumExp(A), B[k]);
						FRACLIKE_ptr -= 4; // because no conditional likelihood
						// ptr_pos+=4;
						//  printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
						//printf("\t\t%.16f\n", B[k]);
					}
					// printf("\t\t\t\t\tB[%d]: %lf\n", k, B[k]);
				}
				FRACLIKE_ptr += 8;
				//for(int l = 0; l < 4; l++)
				//{
				//	printf("\t\t\t\tB[%d]: %lf\n", l, B[l]);
				//}
				C[j] = logSumExp(B);
				//printf("\t\t\t\tC:%lf %lf\n", logSumExp(B), C[j]);
				// printf("\t\t\t\tC:%lf\n", logSumExp(B));
			}
			//for(int l = 0; l < 4; l++)
			//{
			//	printf("\t\t\t\tC[%d]: %lf\n", l, C[l]);
			//}
			Like += logSumExp(C);
			//printf("\t\t\t\tLike:%lf %lf\n", logSumExp(C), Like);
		}
	}

	//printf("\tLikelihood: %.16f\n", Like);
	// printf("\t\t%d,%f,%f,%f,%f\n", seq, (1.0-times[1])*(nodeages[treeNum][seq]+times[2]), times[1], times[2], Like);
	// exit(0);

	return -Like; // Notice: a scaling factor of NUMCAT^(number of sites) is missing
}

// This function is used as the function for Brent's method to estimate the sample age
// Compared to previous versions of estimating the sample age as two layered golden section that
// optimizes the placement of each read for each age, this function will jointly optimize
// the placements of each read and the age
// parameters[0-(n-1)] is the placements of each read
// parameters[n] is the age
// Branch lengths of the read must be calculated based on the age and placement of the reads
double getlike_gamma_root_in_trifucation_sample_age_brent(double parameters[], void *extra_data)
{
	double branch_lengths[numquery];

	double age = *((double *)extra_data);

	//printf("Test age %.16f\n", age);

	// Calculate the branch lengths based on the placement and age bounds
	// Update the placement of the read if branch length < 0 (ie placement
	// and age are incompatible)
	int *tree_ptr = &treeAssign[0];
	int *assign_ptr = &assignments[0];
	double *branch_length_ptr = &branch_lengths[0];
	double *placement_ptr = &parameters[0];

	for (unsigned long int i = 0; i < numquery; i++)
	{
		//printf("\t\tread %lu: %.16f (", i, *placement_ptr);

		// first test if the placement is within the bounds
		if (*placement_ptr < 3e-8)
		{
			*placement_ptr = 3e-8;
		}
		else if (*placement_ptr > bls[*tree_ptr][*assign_ptr] - 3e-8)
		{
			*placement_ptr = bls[*tree_ptr][*assign_ptr] - 3e-8;
		}

		*branch_length_ptr = nodeages[*tree_ptr][*assign_ptr] + *placement_ptr - age;

		// only need to check if placement if below the age to fix since we've
		// already checked the maximum placement above
		if (*branch_length_ptr < 3e-8)
		{
			*placement_ptr = age + 3e-8 - nodeages[*tree_ptr][*assign_ptr];
			*branch_length_ptr = 3e-8;
		}

		//printf("%.16f) %.16f %.16f %.16f\n", *placement_ptr, bls[*tree_ptr][*assign_ptr], nodeages[*tree_ptr][*assign_ptr], *branch_length_ptr);

		// iterate to the next read
		tree_ptr++;
		assign_ptr++;
		branch_length_ptr++;
		placement_ptr++;
	}

	// Re-initialize the pointers
	tree_ptr = &treeAssign[0];
	assign_ptr = &assignments[0];
	branch_length_ptr = &branch_lengths[0];
	placement_ptr = &parameters[0];

	// Iterate through the reads and calculate the likelihood
	int i, j, k, b, c, v, po;
	double Like = 0.0, t[3], A[4], B[4], C[4];
	double *PMAT_ptr_0, *PMAT_ptr_1, *PMAT_ptr_2, *FRACLIKE_ptr;
	unsigned long int *startpos_ptr = &startpos[0];
	unsigned long int *readlength_ptr = &readlength[0];
	//int *QUERYDATA_ptr = &QUERYDATA[0][0]; //UPDATE THIS IF QUERYDATA MADE INTO 2D ARRAY

	// i is seq and j is base
	for (unsigned long int seq = 0; seq < numquery; seq++)
	{
		t[1] = *placement_ptr;									// length from node to position where query joins
		t[0] = *branch_length_ptr;								// length from age of query node to position where query joins
		t[2] = bls[*tree_ptr][*assign_ptr] - *placement_ptr;		// length from position where query joins to parent node

		make_transition_prob_matrices(t, *tree_ptr);

		//FRACLIKE_ptr = &FRACLIKE[*tree_ptr][*assign_ptr * numbases[*tree_ptr] * NUMCAT * 8 + *startpos_ptr * NUMCAT * 8];

		unsigned long long int size = ((unsigned long long)*assign_ptr) * ((unsigned long long)numbases[*tree_ptr]) * ((unsigned long long)NUMCAT) * 8ULL
		                           + ((unsigned long long)*startpos_ptr) * ((unsigned long long)NUMCAT) * 8ULL;
		// printf("FRACLIKE index 7: %d vs %llu\n", *assign_ptr * numbases[*tree_ptr] * NUMCAT * 8 + *startpos_ptr * NUMCAT * 8, size);
		FRACLIKE_ptr = &FRACLIKE[*tree_ptr][size];

		for(unsigned long int i = *startpos_ptr; i < *readlength_ptr + *startpos_ptr; i++)
		{
			po = i - *startpos_ptr;
			b = QUERYDATA[seq][po];

			if(b == -1)
			{
				FRACLIKE_ptr += 32; // NUMCAT (4) * 8;
			}
			else
			{
				// PMAT_ptr_0 = PMAT; 		//0*NUMCAT*4*4
				PMAT_ptr_1 = &PMAT[64];	 // 1*NUMCAT(4)*4*4
				PMAT_ptr_2 = &PMAT[128]; // 2*NUMCAT(4)*4*4
				for (j = 0; j < NUMCAT; j++)
				{
					// PMAT_ptr_0 = &PMAT[j * 4 * 4 + b * 4];
					PMAT_ptr_0 = &PMAT[j * 4 * 4];
					for (k = 0; k < 4; k++)
					{
						// Keep this for when we add errors
						// printf("\t\t\tk: %d\n", k);
						for (v=0; v<4; v++)
						{
							// A[v] = pi[treeNum][v] + PMAT[0][j][v][k] + readlike[seq][po][v];
							A[v] = pi[*tree_ptr][b] + *(PMAT_ptr_0 + v * 4 + k) + readlike[seq][po][v];
						//	//printf("\t\t\t\tpi[%d] is %.16f anf PMAT[0] is %.16f anf reaflike is %.16f\n", v, pi[treeNum][v], PMAT[0][j][v][k], readlike[seq][po][v]);
						}
						B[k] = logSumExp(A);
						//B[k] = pi[*tree_ptr][b] + *PMAT_ptr_0;
						//PMAT_ptr_0++;
						if (*assign_ptr >= numseq[*tree_ptr])
						{ // If not leaf node. Assumes t,c,g,t)leaf nodes are numbered from 0 to numseq-1
							for (v = 0; v < 4; v++)
							{
								A[v] = *PMAT_ptr_2 + *FRACLIKE_ptr;
								PMAT_ptr_2++;
								FRACLIKE_ptr++;
							}
							B[k] += logSumExp(A);
							for (v = 0; v < 4; v++)
							{
								A[v] = *PMAT_ptr_1 + *FRACLIKE_ptr;
								PMAT_ptr_1++;
								FRACLIKE_ptr++;
							}
							B[k] += logSumExp(A);
							FRACLIKE_ptr -= 8;
						}
						else
						{ // If leaf node
							if ((c = DATA[*tree_ptr][*assign_ptr][i]) > -1)
							{
								B[k] += *(PMAT_ptr_1 + c);
								PMAT_ptr_1 += 4;
							}
							for (v = 0; v < 4; v++) // add position here into fraclike
							{
								A[v] = *PMAT_ptr_2 + *FRACLIKE_ptr;
								PMAT_ptr_2++;
								FRACLIKE_ptr++;
							}
							B[k] += logSumExp(A);
							FRACLIKE_ptr -= 4;
						}
					}
					FRACLIKE_ptr += 8;
					C[j] = logSumExp(B);
				}
				Like += logSumExp(C);

			}

		}

		// iterate to the next read
		tree_ptr++;		//treeNum
		assign_ptr++;		//node
		branch_length_ptr++;
		placement_ptr++;
		startpos_ptr++;
		readlength_ptr++;
	}

	//printf("\tLike: %.16f\n", -Like);

	return -Like; // Notice: a scaling factor of NUMCAT^(number of sites) is missing
}

double getlike_gamma_root_in_trifurcation_reassign(double times[3], double parameters[7])
{
	int i, j, k, b, c, v, po, node, seq, treeNum;
	double Like, t[3], A[4], B[4], C[4];
	double *PMAT_ptr_0, *PMAT_ptr_1, *PMAT_ptr_2, *FRACLIKE_ptr;
	// diagonalizaiton has previously been done

	if (onDindic == 1)
	{ // Learn more about what these do, seem to be hitting these values
		if (times[2] < parameters[3] || times[2] > parameters[4])
		{
			return 1000000000.0; // If an actual likelihood is smaler than this we are screwed
		}
	}
	else
	{
		if ((times[2] < parameters[5] || times[2] > parameters[6]) || (times[1] < parameters[3] || times[1] > parameters[4]))
		{
			return 1000000000.0; // If an actual likelihood is smaler than this we are screwed
		}
	}

	node = (int)parameters[2];
	seq = (int)parameters[0];
	treeNum = (int)parameters[1];

	t[1] = times[2];										// length from node to position where query joins
	t[0] = (times[2] + nodeages[treeNum][node]) * times[1]; // length from age of query node to position where query joins
	t[2] = bls[treeNum][node] - times[2];					// length from position where query joins to parent node

	make_transition_prob_matrices(t, treeNum);

	Like = 0.0;
	// fraclike pointer starts where the read starts in the alignment at the right node
	//FRACLIKE_ptr = &FRACLIKE[treeNum][node * numbases[treeNum] * NUMCAT * 8 + startposTemp[treeNum][seq] * NUMCAT * 8];

	unsigned long long int size = ((unsigned long long)node) * ((unsigned long long)numbases[treeNum]) * ((unsigned long long)NUMCAT) * 8ULL
						+ ((unsigned long long)startposTemp[treeNum][seq]) * ((unsigned long long)NUMCAT) * 8ULL;

	// printf("FRACLIKE index 8: %d vs %llu\n", node * numbases[treeNum] * NUMCAT * 8 + startposTemp[treeNum][seq] * NUMCAT * 8, size);
	FRACLIKE_ptr = &FRACLIKE[treeNum][size];
	// printf("%d\t%d\n", numbases[treeNum],node * numbases[treeNum] * NUMCAT * 8 + startposTemp[treeNum][seq] * NUMCAT * 8);
	//  printf("Sequence %d\n", seq);

	for (i = startposTemp[treeNum][seq]; i < readLengthTemp[treeNum][seq] + startposTemp[treeNum][seq]; i++)
	{
		po = i - startposTemp[treeNum][seq]; // i keeps track of the positon in the ref sequences while po is the positoon in the read
		// b = QUERYDATA[seq][po];
		b = readsTreeSorted[treeNum][seq][po];

		// printf("\tBase %d\n", b);

		if (b == -1)
		{
			FRACLIKE_ptr += 32; // NUMCAT(4) * 8;
		}
		else
		{
			// if (b!=-1){
			// PMAT_ptr_0 = PMAT; 		//0*NUMCAT*4*4
			PMAT_ptr_1 = &PMAT[64];	 // 1*NUMCAT(4)*4*4
			PMAT_ptr_2 = &PMAT[128]; // 2*NUMCAT(4)*4*4
			for (j = 0; j < NUMCAT; j++)
			{
				// printf("\t\tj: %d\n", j);
				PMAT_ptr_0 = &PMAT[j * 4 * 4 + b * 4];
				for (k = 0; k < 4; k++)
				{
					// printf("\t\t\tk: %d\n", k);
					for (v=0; v<4; v++)
					{
					     // A[v] = pi[treeNum][v] + PMAT[0][j][v][k] + readLikeTemp[treeNum][seq][po][v];
						A[v] = pi[treeNum][v] + *(PMAT_ptr_0 + v * 4 + k) + readLikeTemp[treeNum][seq][po][v];
					//	//printf("\t\t\t\tpi[%d] is %.16f anf PMAT[0] is %.16f anf reaflike is %.16f\n", v, pi[treeNum][v], PMAT[0][j][v][k], readLikeTemp[treeNum][seq][po][v]);
					}
					B[k] = logSumExp(A);
					//B[k] = pi[treeNum][b] + *PMAT_ptr_0;
					//// printf("\t\t\t\t\tB:%lf\tpi:%lf\tPMAT[0]:%lf\n", B[k], pi[treeNum][b], *PMAT_ptr_0);
					//PMAT_ptr_0++;
					if (node >= numseq[treeNum])
					{ // If not leaf node. Assumes t,c,g,t)leaf nodes are numbered from 0 to numseq-1
						for (v = 0; v < 4; v++)
						{
							// printf("\t\t\t\t\tPMAT[2] is %.16f and FRACLIKE[%d] is %.16f\n", *PMAT_ptr_2, v, FRACLIKE[treeNum][i][node][j][v]);
							// A[v] = *PMAT_ptr_2 + FRACLIKE[treeNum][i][node][j][v];
							// printf("\t\tFRACLIKE[%d][%d][%d][%d][%d]:%.16f\t%d\n", treeNum, i, node, j, v, *FRACLIKE_ptr, ptr_pos);
							A[v] = *PMAT_ptr_2 + *FRACLIKE_ptr;
							PMAT_ptr_2++;
							FRACLIKE_ptr++;
							// ptr_pos++;
						}
						B[k] += logSumExp(A);

						for (v = 0; v < 4; v++)
						{
							// printf("\t\t\t\t\tPMAT[1] is %.16f and FRACLIKE[%d] is %.16f\n", *PMAT_ptr_1, v+4, FRACLIKE[treeNum][i][node][j][v+4]);
							// printf("\t\t\t\t\tPMAT[1][%d][%d][%d] is %.16f\n", j, k, v, *PMAT_ptr_1);
							// A[v] = *PMAT_ptr_1 + FRACLIKE[treeNum][i][node][j][v+4];
							// printf("\t\tFRACLIKE[%d][%d][%d][%d][%d]:%.16f\t%d\n", treeNum, i, node, j, v, *FRACLIKE_ptr, ptr_pos);
							A[v] = *PMAT_ptr_1 + *FRACLIKE_ptr;
							PMAT_ptr_1++;
							FRACLIKE_ptr++;
							// ptr_pos++;
						}
						B[k] += logSumExp(A);
						FRACLIKE_ptr -= 8;
						// printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
						// printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
					}
					else
					{ // If leaf node
						if ((c = DATA[treeNum][node][i]) > -1)
						{
							B[k] += *(PMAT_ptr_1 + c);
							// printf("%.16f\t%d %d\n", *(PMAT_ptr_1 + c), c, 5-c);
							PMAT_ptr_1 += 4;
						}
						for (v = 0; v < 4; v++) // add position here into fraclike
						{
							// printf("\t\t\t\t\tPMAT[2] is %.16f and FRACLIKE[%d] is %.16f\n", *PMAT_ptr_2, v, FRACLIKE[treeNum][i][node][j][v]);
							// A[v] = *PMAT_ptr_2 + FRACLIKE[treeNum][i][node][j][v];
							// printf("\t\tFRACLIKE[%d][%d][%d][%d][%d]:%.16f\t%d\n", treeNum, i, node, j, v, *FRACLIKE_ptr, ptr_pos);
							A[v] = *PMAT_ptr_2 + *FRACLIKE_ptr;
							PMAT_ptr_2++;
							FRACLIKE_ptr++;
							// ptr_pos++;
						}
						B[k] += logSumExp(A);
						FRACLIKE_ptr -= 4; // because no conditional likelihood
										   // ptr_pos+=4;
										   //  printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
					}
					// printf("\t\t\t\t\tB[%d]: %lf\n", k, B[k]);
				}
				FRACLIKE_ptr += 8;
				// for(int l = 0; l < 4; l++)
				// {
				// 	printf("\t\t\t\tB[%d]: %lf\n", l, B[l]);
				// }
				C[j] = logSumExp(B);
				// printf("\t\t\t\tC:%lf\n", logSumExp(B));
			}
			Like += logSumExp(C);
			// printf("\t\t\t\tLike:%lf\n", logSumExp(C));
		}
	}

	// printf("\tLikelihood: %.16f\n", Like);
	// printf("\t\t%d,%f,%f,%f,%f\n", seq, (1.0-times[1])*(nodeages[treeNum][seq]+times[2]), times[1], times[2], Like);

	return -Like; // Notice: a scaling factor of NUMCAT^(number of sites) is missing
}

double getlike_gamma_root_in_trifurcation_Print_Lik(double times[3], double parameters[7])
{
	int i, j, k, b, c, v, po, node, seq, treeNum;
	double Like, t[3], A[4], B[4], C[4];
	double *PMAT_ptr_0, *PMAT_ptr_1, *PMAT_ptr_2, *FRACLIKE_ptr;
	// diagonalizaiton has previously been done

	if (onDindic == 1)
	{ // Learn more about what these do, seem to be hitting these values
		if (times[2] < parameters[3] || times[2] > parameters[4])
		{
			// printf("\t\tError 1: %.16f < %.16f or %.16f > %.16f\n", times[2], parameters[3], times[2], parameters[4]);
			return 1000000000.0; // If an actual likelihood is smaler than this we are screwed
		}
	}
	else
	{
		if ((times[2] < parameters[5] || times[2] > parameters[6]) || (times[1] < parameters[3] || times[1] > parameters[4]))
		{
			// printf("\t\tError 2: (%.16f < %.16f or %.16f > %.16f) or (%.16f < %.16f or %.16f > %.16f)\n", times[2], parameters[5], times[2], parameters[6], times[1], parameters[3], times[1], parameters[4]);
			return 1000000000.0; // If an actual likelihood is smaler than this we are screwed
		}
	}

	node = (int)parameters[2];
	seq = (int)parameters[0];
	treeNum = (int)parameters[1];
	//	printf("Analysing read %d\n", seq);
	// Why copy this? TO DO: Remove this copy and change below code to reflect the change
	// for (i=0; i<4; i++)
	//	printf("pi[%d]: %.16lf\t", i, pi[treeNum][i]);
	//	pi[i]=parameters[i+3];

	// printf("\n");

	t[1] = times[2];										// length from node to position where query joins
	t[0] = (times[2] + nodeages[treeNum][node]) * times[1]; // length from age of query node to position where query joins
	t[2] = bls[treeNum][node] - times[2];					// length from position where query joins to parent node

	// printf("node: %d\tseq: %d\ttree: %d\tt0: %.16f\tt1: %.16f\tt2: %.16f\n", node, seq, treeNum, t[0], t[1], t[2]);
	// printf("t0: %.16f\tt1: %.16f\tt2: %.16f\n", t[0], t[1], t[2]);

	make_transition_prob_matrices(t, treeNum);

	// printf("Errors setting of %d, sequence %d, assignment %d\n", errors, seq, node);

	Like = 0.0;

	// fraclike pointer starts where the read starts in the alignment at the right node
	//FRACLIKE_ptr = &FRACLIKE[treeNum][node * numbases[treeNum] * NUMCAT * 8 + startpos[seq] * NUMCAT * 8];
	unsigned long long int size = ((unsigned long long)node) * ((unsigned long long)numbases[treeNum]) * ((unsigned long long)NUMCAT) * 8ULL
    									+ ((unsigned long long)startpos[seq]) * ((unsigned long long)NUMCAT) * 8ULL;
	// printf("FRACLIKE index 9: %d vs %llu\n", node * numbases[treeNum] * NUMCAT * 8 + startpos[seq] * NUMCAT * 8, size);
	FRACLIKE_ptr = &FRACLIKE[treeNum][size];

	// printf("Sequence %d\n", seq);

	for (i = startpos[seq]; i < readlength[seq] + startpos[seq]; i++)
	{
		po = i - startpos[seq]; // i keeps track of the positon in the ref sequences while po is the positoon in the read
		b = QUERYDATA[seq][po];

		// printf("\tBase %d of pos %d with start %d and readlength %d\n", b, po, startpos[seq], readlength[seq]);

		if (b == -1)
		{
			FRACLIKE_ptr += 32; // NUMCAT (4) * 8
		}
		else
		{
			// PMAT_ptr_0 = PMAT; 		//0*NUMCAT*4*4
			PMAT_ptr_1 = &PMAT[64];	 // 1*NUMCAT*4*4 (4 * 4 * 4)
			PMAT_ptr_2 = &PMAT[128]; // 2*NUMCAT*4*4 (2 * 4 * 4 * 4)
			for (j = 0; j < NUMCAT; j++)
			{
				PMAT_ptr_0 = &PMAT[j * 16 + b * 4];
				// printf("\t\tj: %d\n", j);
				for (k = 0; k < 4; k++)
				{
					// printf("\t\t\tk: %d\n", k);
					for (v=0; v<4; v++)
					{
						//A[v] = pi[treeNum][v] + PMAT[0][j][v][k] + readlike[seq][po][v];
						A[v] = pi[treeNum][v] + *(PMAT_ptr_0 + v * 4 + k) + readlike[seq][po][v];
						//printf("\t\t\t\tpi[%d] is %.16f anf PMAT[0] is %.16f anf reaflike is %.16f\n", v, pi[treeNum][v], PMAT[0][j][v][k], readlike[seq][po][v]);
					}
					B[k] = logSumExp(A);
					// B[k] = pi[treeNum][b] + *PMAT_ptr_0;
					// // printf("\t\t\t\t\tB:%lf\tpi:%lf\tPMAT[0]:%lf\n", B[k], pi[treeNum][b], PMAT[0][j][b][k]);
					// PMAT_ptr_0++;
					if (node >= numseq[treeNum])
					{ // If not leaf node. Assumes the leaf nodes are numbered from 0 to numseq-1
						for (v = 0; v < 4; v++)
						{
							// printf("\t\t\t\t\tPMAT[2] is %.16f and FRACLIKE[%d] is %.16f\n", PMAT[2][j][k][v], v, FRACLIKE[treeNum][i][node][j][v]);
							// A[v] = *PMAT_ptr_2 + FRACLIKE[treeNum][i][node][j][v];
							A[v] = *PMAT_ptr_2 + *FRACLIKE_ptr;
							PMAT_ptr_2++;
							FRACLIKE_ptr++;
						}
						B[k] += logSumExp(A);
						for (v = 0; v < 4; v++)
						{
							// printf("\t\t\t\t\tPMAT[1] is %.16f and FRACLIKE[%d] is %.16f\n", PMAT[1][j][k][v], v+4, FRACLIKE[treeNum][i][node][j][v+4]);
							// A[v] = *PMAT_ptr_1 + FRACLIKE[treeNum][i][node][j][v+4];
							A[v] = *PMAT_ptr_1 + *FRACLIKE_ptr;
							PMAT_ptr_1++;
							FRACLIKE_ptr++;
						}
						B[k] += logSumExp(A);
						FRACLIKE_ptr -= 8;
						// printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
						// printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
					}
					else
					{ // If leaf node
						if ((c = DATA[treeNum][node][i]) > -1)
						{
							B[k] += *(PMAT_ptr_1 + c);
							PMAT_ptr_1 += 4;
						}
						for (v = 0; v < 4; v++) // add position here into fraclike
						{
							// printf("\t\t\t\t\tPMAT[2] is %.16f and FRACLIKE[%d] is %.16f\n", PMAT[2][j][k][v], v, FRACLIKE[treeNum][i][node][j][v]);
							// A[v] = *PMAT_ptr_2 + FRACLIKE[treeNum][i][node][j][v];
							A[v] = *PMAT_ptr_2 + *FRACLIKE_ptr;
							PMAT_ptr_2++;
							FRACLIKE_ptr++;
						}
						B[k] += logSumExp(A);
						FRACLIKE_ptr -= 4;
						// printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
					}
					// printf("\t\t\t\t\tB[%d]: %lf\n", k, B[k]);
				}
				FRACLIKE_ptr += 8;
				// for(int l = 0; l < 4; l++)
				// {
				// 	printf("\t\t\t\tB[%d]: %lf\n", l, B[l]);
				// }
				C[j] = logSumExp(B);
				// printf("\t\t\t\tC:%lf\n", logSumExp(B));
			}
			Like += logSumExp(C);
			// printf("\t\t\t\tLike:%lf\n", logSumExp(C));
		}
	}

	// printf("\tLikelihood: %.16f\n", Like);
	printf("\t\t%d,%d,%d,%.16f,%.16f,%.16f,%.16f\n", seq, node, treeNum, (1.0 - times[1]) * (nodeages[treeNum][seq] + times[2]), times[1], times[2], Like);

	return -Like; // Notice: a scaling factor of NUMCAT^(number of sites) is missing
}

// Same as above, but used to print out site scores
double getlike_gamma_root_in_trifurcation_Print(double times[3], double parameters[3])
{
	int i, j, k, b, c, v, po, node, seq, treeNum;
	double Like, t[3], A[4], B[4], C[4];
	double *PMAT_ptr_0, *PMAT_ptr_1, *PMAT_ptr_2, *FRACLIKE_ptr;

	// Don't need to do the error checking like normally

	node = (int)parameters[2];
	seq = (int)parameters[0];
	treeNum = (int)parameters[1];
	// for (i=0; i<4; i++)
	//	pi[i]=parameters[i+3];

	t[1] = times[2];										// length from node to position where query joins
	t[0] = (times[2] + nodeages[treeNum][node]) * times[1]; // length from age of query node to position where query joins
	t[2] = bls[treeNum][node] - times[2];					// length from position where query joins to parent node

	// printf("node: %d\tseq: %d\ttree: %d\tt0: %.16f\tt1: %.16f\tt2: %.16f\n", node, seq, treeNum, t[0], t[1], t[2]);

	make_transition_prob_matrices(t, treeNum);

	Like = 0.0;
	// fraclike pointer starts where the read starts in the alignment at the right node
	//FRACLIKE_ptr = &FRACLIKE[treeNum][node * numbases[treeNum] * NUMCAT * 8 + startpos[seq] * NUMCAT * 8];
	unsigned long long int size = ((unsigned long long)node) * ((unsigned long long)numbases[treeNum]) * ((unsigned long long)NUMCAT) * 8ULL
 	                         + ((unsigned long long)startpos[seq]) * ((unsigned long long)NUMCAT) * 8ULL;
	// printf("FRACLIKE index 10: %d vs %llu\n", node * numbases[treeNum] * NUMCAT * 8 + startpos[seq] * NUMCAT * 8, size);
	FRACLIKE_ptr = &FRACLIKE[treeNum][size];

	for (i = startpos[seq]; i < readlength[seq] + startpos[seq]; i++)
	{
		po = i - startpos[seq];
		b = QUERYDATA[seq][po];

		if (b == -1)
		{
			FRACLIKE_ptr += 32; // NUMCAT(4) * 8;
		}
		else
		{
			// PMAT_ptr_0 = PMAT; 		//0*NUMCAT*4*4
			PMAT_ptr_1 = &PMAT[64];	 // 1*NUMCAT(4)*4*4
			PMAT_ptr_2 = &PMAT[128]; // 2*NUMCAT(4)*4*4

			for (j = 0; j < NUMCAT; j++)
			{
				PMAT_ptr_0 = &PMAT[j * 4 * 4 + b * 4];
				for (k = 0; k < 4; k++)
				{
					for (v=0; v<4; v++)
					{
						//A[v] = pi[treeNum][v] + PMAT[0][j][v][k] + readlike[seq][po][v];
						A[v] = pi[treeNum][v] + *(PMAT_ptr_0 + v * 4 + k) + readlike[seq][po][v];
					}
					 B[k] = logSumExp(A);
					//B[k] = pi[treeNum][b] + *PMAT_ptr_0;
					//PMAT_ptr_0++;
					if (node >= numseq[treeNum])
					{ // If not leaf node. Assumes the leaf nodes are numbered from 0 to numseq-1
						for (v = 0; v < 4; v++)
						{
							A[v] = *PMAT_ptr_2 + *FRACLIKE_ptr;
							PMAT_ptr_2++;
							FRACLIKE_ptr++;
						}
						B[k] += logSumExp(A);
						for (v = 0; v < 4; v++)
						{
							A[v] = *PMAT_ptr_1 + *FRACLIKE_ptr;
							PMAT_ptr_1++;
							FRACLIKE_ptr++;
						}
						B[k] += logSumExp(A);
						FRACLIKE_ptr -= 8;
					}
					else
					{ // If leaf node
						if ((c = DATA[treeNum][node][i]) > -1)
						{
							B[k] += *(PMAT_ptr_1 + c);
							PMAT_ptr_1 += 4;
						}
						for (v = 0; v < 4; v++) // add position here into fraclike
						{
							A[v] = *PMAT_ptr_2 + *FRACLIKE_ptr;
							PMAT_ptr_2++;
						}
						B[k] += logSumExp(A);
						FRACLIKE_ptr -= 4;
					}
				}
				FRACLIKE_ptr += 8;
				C[j] = logSumExp(B);
			}
			printf("%d, %.16f\n", i, logSumExp(C));
			Like += logSumExp(C);
		}
	}

	return -Like; // Notice: a scaling factor of NUMCAT^(number of sites) is missing
}

double getlike_gamma_root_in_trifurcation_L0(double times[3], double parameters[7])
{

	double T[3];

	T[2] = times[1];
	T[1] = 1.0;
	return getlike_gamma_root_in_trifurcation(T, parameters);
}

double getlike_gamma_root_in_trifurcation_testAge_reassign(double rootPlace, double parameters[8])
{
	double T[3], param[7];

	T[2] = rootPlace;
	T[1] = parameters[3];

	// fix so that it can work
	param[0] = parameters[0];
	param[1] = parameters[1];
	param[2] = parameters[2];

	if (T[1] < 0.0)
	{
		// printf("\t\ttestAge of %.16f: T[1] < 0.0 of %.16f!!\n", testAge, T[1]);
		return 1000000000.0;
	}
	else if (T[1] > 1.0) // This must be between 0 and 1
	{
		// printf("\t\ttestAge of %.16f: T[1] > 1.0 of %.16f!!\n", testAge, T[1]);
		return 1000000000.0;
	}

	for (int i = 4; i < 8; i++)
	{
		param[i - 1] = parameters[i];
	}

	// printf("\t\t\tT[1]: %.16f, T[2]: %.16f\n", T[1], T[2]);

	return getlike_gamma_root_in_trifurcation_reassign(T, param);
}

double getlike_gamma_root_in_trifurcation_testAge(double rootPlace, double parameters[7])
{
	double T[3], returnValue;

	T[2] = rootPlace;
	T[1] = 1.0 - testAge / (nodeages[(int)parameters[1]][(int)parameters[2]] + T[2]);

	// Since its been switched to 2D optimization, allow for the bounds to catch <0 and >1
	if (T[1] < 0.0)
	{
		// This means that the testAge > node + T[2], requiring a "negative" branch. So changing to 0 does not make sense and should not be done
		// printf("\t\ttestAge of %.16f: T[1] < 0.0 of %.16f!!\n", testAge, T[1]);

		return 1000000000.0;
	}
	else if (T[1] > 1.0) // This must be between 0 and 1
	{
		// This only would happen if negative T[2] + node age occurs, which shouldn't? Again, makes more sense to penalize heavily with high likelihood instead of changing to 1
		// printf("\t\ttestAge of %.16f: T[1] > 1.0 of %.16f!!\n", testAge, T[1]);
		return 1000000000.0;
	}

	return getlike_gamma_root_in_trifurcation(T, parameters);
}

double getlike_gamma_root_in_trifurcation_testRoot(double branchLength, double parameters[7])
{
	double T[3], returnValue;

	// printf("Branch length %.16f and rooted %.16f\n", branchLength, rooted);

	T[2] = rooted;
	T[1] = branchLength;

	if (T[1] < 0.0)
	{
		return 1000000000.0;
	}
	else if (T[1] > 1.0) // This must be between 0 and 1
	{
		return 1000000000.0;
	}

	return getlike_gamma_root_in_trifurcation(T, parameters);
}

void inittransitionmatrix()

{
	int i, j;
	double sum, RIVAL[4], RIVEC[4][4], A[4][4], workspace[8], norm;

	// TO DO: Can probably make these a single malloc
	RRVAL = (double **)malloc(numTrees * sizeof(double *));
	RRVEC = (double **)malloc(numTrees * sizeof(double *));
	LRVEC = (double **)malloc(numTrees * sizeof(double *));

	for (unsigned long int treeNum = 0; treeNum < numTrees; treeNum++)
	{
		if (usedTrees[treeNum] == 0)
		{
			// no reads were assigned to this tree, can skip
			continue;
		}
		RRVAL[treeNum] = (double *)malloc(4 * sizeof(double));
		RRVEC[treeNum] = (double *)malloc(4 * 4 * sizeof(double));
		//(double **)malloc(4 * sizeof(double *));
		LRVEC[treeNum] = (double *)malloc(4 * 4 * sizeof(double));
		//(double **)malloc(4 * sizeof(double *));
		// for(i = 0; i < 4; i++)
		//{
		//	RRVEC[treeNum][i] = (double *)malloc(4 * sizeof(double));
		//	LRVEC[treeNum][i] = (double *)malloc(4 * sizeof(double));
		// }

		// 2 (a pia pic + b pia pig + c pia pit + d pic pig + e pic pit + f pig pit)
		norm = 2.0 * (par[treeNum][0] * pi[treeNum][0] * pi[treeNum][1] + par[treeNum][1] * pi[treeNum][0] * pi[treeNum][2] + par[treeNum][2] * pi[treeNum][0] * pi[treeNum][3] +
					  par[treeNum][3] * pi[treeNum][1] * pi[treeNum][2] + par[treeNum][4] * pi[treeNum][1] * pi[treeNum][3] +
					  par[treeNum][5] * pi[treeNum][2] * pi[treeNum][3]);

		A[0][1] = pi[treeNum][1] * par[treeNum][0];
		A[0][2] = pi[treeNum][2] * par[treeNum][1];
		A[0][3] = pi[treeNum][3] * par[treeNum][2];
		A[1][0] = pi[treeNum][0] * par[treeNum][0];
		A[1][2] = pi[treeNum][2] * par[treeNum][3];
		A[1][3] = pi[treeNum][3] * par[treeNum][4];
		A[2][0] = pi[treeNum][0] * par[treeNum][1];
		A[2][1] = pi[treeNum][1] * par[treeNum][3];
		A[2][3] = pi[treeNum][3] * par[treeNum][5];
		A[3][0] = pi[treeNum][0] * par[treeNum][2];
		A[3][1] = pi[treeNum][1] * par[treeNum][4];
		A[3][2] = pi[treeNum][2] * par[treeNum][5];

		// printf("GTR Free Parameter: %lf, %lf, %lf, %lf, %lf, %lf\n", par[treeNum][0], par[treeNum][1], par[treeNum][2], par[treeNum][3], par[treeNum][4],par[treeNum][5]);

		// printf("-q_ii\n");
		for (i = 0; i < 4; i++)
		{
			// Can't you just subtract directly from A[i][i]??
			A[i][i] = 0.0;
			sum = 0.0;
			for (j = 0; j < 4; j++)
			{
				sum = sum + A[i][j];
				A[i][j] = A[i][j] / norm;
			}
			A[i][i] = -sum / norm;
			// printf("%d\t%lf\n", i, A[i][i]);
		}

		// printf("Pi - 0: %lf\t1: %lf\t2: %lf\t3: %lf\n", pi[treeNum][0], pi[treeNum][1], pi[treeNum][2], pi[treeNum][3]);

		// Edit here to get Q (stored as A[][]) and pi/frequencies to calculate the estimated number of differences
		if (eigen(1, A[0], 4, RRVAL[treeNum], RIVAL, RRVEC[treeNum], RIVEC[0], workspace) != 0)
		{
			printf("Transitions matrix did not converge or contained non-real values!\n");
			exit(-1);
		}
		for (i = 0; i < 4; i++)
		{
			for (j = 0; j < 4; j++)
			{
				LRVEC[treeNum][i * 4 + j] = RRVEC[treeNum][i * 4 + j];
				// printf("LRVEC[%d][%d][%d]: %lf\n", treeNum, i, j, LRVEC[treeNum][i * 4 + j]);
			}
		}
		if (matinv(RRVEC[treeNum], 4, 4, workspace) != 0) // TO DO: MAKE SURE THE RRVEC INPUT IS CORRECT FOR THIS FUNCTION
			printf("Could not invert matrix!\nResults may not be reliable!\n");
	}
}

/*some old code for reading a fasta file and storing DNA as ints*/
// Probably needs same fix as get_frac_like.c's readseq
int readseq(int *nb, int treeNum)
{
	int i, j, num, k = 0;
	char c;

	(void)fscanf(infile, "%i %i", &num, nb);
	// printf("there are %i species and %i bases in tree %i\n",num,*nb, treeNum);
	if (num < 3)
	{
		printf("This is for more than two seqs only!\n");
		exit(-1);
	}
	DATA[treeNum] = (int **)malloc(num * (sizeof(int *)));
	for (i = 0; i < num; i++)
		DATA[treeNum][i] = (int *)malloc(*nb * (sizeof(int)));
	while ((c = (fgetc(infile))) != '\n')
		;
	do
	{
		for (i = 0; i < num; i++)
		{
			j = k;
			while ((c = tolower((fgetc(infile)))) != '\n')
			{
				if ((c != ' ') && (c != '\t'))
				{
					if (c == 'a')
						DATA[treeNum][i][j] = 0;
					else if (c == 'c')
						DATA[treeNum][i][j] = 1;
					else if (c == 'g')
						DATA[treeNum][i][j] = 2;
					else if (c == 't')
						DATA[treeNum][i][j] = 3;
					else if (c == 'n' || c == '-' || c == '~')
						DATA[treeNum][i][j] = -1;
					else
					{
						printf("\nBAD BASE (%c) in species %i base %i", c, i + 1, j + 1);
						(void)scanf("%i", &i);
						exit(-1);
					}
					j++;
					if (i == (num - 1))
						k++;
				}
			}
		}
	} while (k < *nb);
	return num;
}

int isblankorreturn(char c)
{
	if (c < 33)
		return 1;
	else
		return 0;
}

/* ---- chunk-buffer helpers ---- */
#define CBUF_SIZE 65536
typedef struct { char buf[CBUF_SIZE + 1]; int pos, len; } ChunkBuf;

static void cbuf_refill(ChunkBuf *cb, FILE *f)
{
	int rem = cb->len - cb->pos;
	if (rem > 0) memmove(cb->buf, cb->buf + cb->pos, rem);
	cb->pos = 0;
	cb->len = rem + (int)fread(cb->buf + rem, 1, CBUF_SIZE - rem, f);
	cb->buf[cb->len] = '\0';
}

static char cbuf_next_nonws(ChunkBuf *cb, FILE *f)
{
	for (;;) {
		while (cb->pos < cb->len) {
			char c = cb->buf[cb->pos++];
			if (!isblankorreturn(c)) return c;
		}
		cbuf_refill(cb, f);
		if (cb->len == 0) return '\0';
	}
}

static int cbuf_next_int(ChunkBuf *cb, FILE *f)
{
	while (cb->pos < cb->len && isblankorreturn(cb->buf[cb->pos])) cb->pos++;
	if (cb->len - cb->pos < 32) cbuf_refill(cb, f);
	char *end;
	int val = (int)strtol(cb->buf + cb->pos, &end, 10);
	cb->pos = (int)(end - cb->buf);
	return val;
}

static double cbuf_next_double(ChunkBuf *cb, FILE *f)
{
	while (cb->pos < cb->len && isblankorreturn(cb->buf[cb->pos])) cb->pos++;
	if (cb->len - cb->pos < 32) cbuf_refill(cb, f);
	char *end;
	double val = strtod(cb->buf + cb->pos, &end);
	cb->pos = (int)(end - cb->buf);
	return val;
}

static unsigned long cbuf_next_ulong(ChunkBuf *cb, FILE *f)
{
	while (cb->pos < cb->len && isblankorreturn(cb->buf[cb->pos])) cb->pos++;
	if (cb->len - cb->pos < 32) cbuf_refill(cb, f);
	char *end;
	unsigned long val = strtoul(cb->buf + cb->pos, &end, 10);
	cb->pos = (int)(end - cb->buf);
	return val;
}
/* ---- end chunk-buffer helpers ---- */

// Probably needs same fix as get_frac_like.c
// Fix so that we read reads into tree sorted datastructure
// Also move final reads into QUERYDATA
int read_query_data(int num)
{
	int i, j, k = 0;
	char c;
	int *numTreeAssigned = calloc(numTrees, sizeof(int));

	ChunkBuf cb;
	cb.pos = 0; cb.len = 0;
	cbuf_refill(&cb, infile);

	for (i = 0; i < num; i++)
	{
		int length = cbuf_next_int(&cb, infile);
		startposTemp[treeAssign[i]][numTreeAssigned[treeAssign[i]]] =
			cbuf_next_ulong(&cb, infile);
		readLengthTemp[treeAssign[i]][numTreeAssigned[treeAssign[i]]] = length;

		// Store reads based on tree assignment
		readsTreeSorted[treeAssign[i]][numTreeAssigned[treeAssign[i]]] = malloc(length * (sizeof(int)));
		readOrder[i] = numTreeAssigned[treeAssign[i]];
		tempAssignments[treeAssign[i]][numTreeAssigned[treeAssign[i]]] = assignments[i];
		numReadsPerAssign[treeAssign[i]][assignments[i]]++;
		numTreeAssigned[treeAssign[i]]++;

		for (j = 0; j < length; j++)
		{
			c = tolower(cbuf_next_nonws(&cb, infile));
			if (c == 'a')
				readsTreeSorted[treeAssign[i]][readOrder[i]][j] = 0;
			else if (c == 'c')
				readsTreeSorted[treeAssign[i]][readOrder[i]][j] = 1;
			else if (c == 'g')
				readsTreeSorted[treeAssign[i]][readOrder[i]][j] = 2;
			else if (c == 't')
				readsTreeSorted[treeAssign[i]][readOrder[i]][j] = 3;
			else if (c == 'n' || c == '-' || c == '~')
				readsTreeSorted[treeAssign[i]][readOrder[i]][j] = -1;
			else
			{
				printf("\nBAD BASE (%c) in query %i base %i", c, i + 1, j + 1);
				exit(-1);
			}
		}
	}

	////print out new reads
	// for (i = 0; i < num; i++)
	//{
	//	for(j = 0; j < readLengthTemp[treeAssign[i]][readOrder[i]]; j++)
	//	{
	//		if (readsTreeSorted[treeAssign[i]][readOrder[i]][j] == 0)
	//			printf("A");
	//		else if (readsTreeSorted[treeAssign[i]][readOrder[i]][j] == 1)
	//			printf("C");
	//		else if (readsTreeSorted[treeAssign[i]][readOrder[i]][j] == 2)
	//			printf("G");
	//		else if (readsTreeSorted[treeAssign[i]][readOrder[i]][j] == 3)
	//			printf("T");
	//		else
	//			printf("-");
	//	}
	//	printf("\n");
	// }

	free(numTreeAssigned);

	return num;
	// printf("Errors setting of %d, sequence %d, assignment %d\n", errors, seq, node);
}

void get_fractionalike(unsigned long int treeNum)
{
	int i, j, k, v, inin;
	char c;
	unsigned long long idx;

	unsigned long long int size = ((long long int)(2 * numseq[treeNum] - 1)) * (long long int)numbases[treeNum] * (long long int)NUMCAT * 8LL;

	FRACLIKE[treeNum] = (double *)calloc(size, sizeof(double));

	ChunkBuf cb;
	cb.pos = 0; cb.len = 0;
	cbuf_refill(&cb, infile);

	for (i = 0; i < NUMCAT; i++)
	{
		c = cbuf_next_nonws(&cb, infile);
		if (c != 'C')
		{
			printf("error reading fractional likelihoods (C%i: %c != C)\n", i, c);
			exit(-1);
		}
		inin = cbuf_next_int(&cb, infile);
		if (inin != i + 1)
		{
			printf("error reading fractional likelihoods (c%i: %i != %i)\n", i, inin, i + 1);
			exit(-1);
		}
		for (j = 0; j < numbases[treeNum]; j++)
		{
			c = cbuf_next_nonws(&cb, infile);
			if (c != 'S')
			{
				printf("error reading fractional likelihoods (s%i: '%c' != S)\n", j, c);
				exit(-1);
			}
			inin = cbuf_next_int(&cb, infile);
			if (inin != j + 1)
			{
				printf("error reading fractional likelihoods(s%i: %i != %i)\n", j, inin, j);
				exit(-1);
			}
			c = cbuf_next_nonws(&cb, infile);
			if (c != ':')
			{
				printf("error reading fractional likelihoods (s%i: '%c' != :)\n", j, c);
				exit(-1);
			}

			// leaf nodes first
			for (k = 0; k < numseq[treeNum]; k++)
			{
				for (v = 0; v < 4; v++)
				{
					idx = ((unsigned long long)k) * numbases[treeNum] * NUMCAT * 8
					    + ((unsigned long long)j) * NUMCAT * 8
					    + ((unsigned long long)i) * 8
					    + v;
					FRACLIKE[treeNum][idx] = cbuf_next_double(&cb, infile);
				}
			}

			// non-leaf nodes with 8
			for (k = numseq[treeNum]; k < 2 * numseq[treeNum] - 1; k++)
			{
				// For fractional likelihood
				for (v = 0; v < 4; v++)
				{
					idx = ((unsigned long long)k) * numbases[treeNum] * NUMCAT * 8
					    + ((unsigned long long)j) * NUMCAT * 8
					    + ((unsigned long long)i) * 8
					    + v;
					FRACLIKE[treeNum][idx] = cbuf_next_double(&cb, infile);
				}

				// For conditional likelihood
				for (v = 0; v < 4; v++)
				{
					idx = ((unsigned long long)k) * numbases[treeNum] * NUMCAT * 8
					    + ((unsigned long long)j) * NUMCAT * 8
					    + ((unsigned long long)i) * 8
					    + v + 4;
					FRACLIKE[treeNum][idx] = cbuf_next_double(&cb, infile);
				}
			}
		}
	}


	// double *test_ptr = &FRACLIKE[treeNum][0];
	// for(int i = 0; i < 2 * numseq[treeNum] - 1; i++)
	//{
	//	//nodes
	//	for(int j = 0; j < numbases[treeNum]; j++)
	//	{
	//		//pos
	//		for(int k = 0; k < NUMCAT; k++)
	//		{
	//			//NUMCAT
	//			for(int l = 0; l < 8; l++)
	//			{
	//				printf("%d %d %d %d %.16f vs %d\n", i, j, k, l, *test_ptr, i * numbases[treeNum] * NUMCAT * 8 + j * NUMCAT * 8 + k * 8 + l);
	//				test_ptr++;
	//			}
	//		}
	//	}
	// }

	// int fracCount = 0;
	// for(i = 0; i < numTrees; i++)
	//{
	//	for(j = 0; j < 2 * numseq[i] - 1; j++)
	//	{
	//		for(int l = 0; l < numbases[i]; l++)
	//		{
	//			for(k = 0; k < NUMCAT; k++)
	//			{
	//				for(v = 0; v < 8; v++)
	//				{
	//					printf("%d: %.16f\n", fracCount,
	//								FRACLIKE[i][fracCount]);
	//					fracCount++;
	//				}
	//			}
	//		}
	//	}

	//}

	// exit(0);
}

// How likely read is not error
// make default error 0.01
void make_readfraclike()
{
	int i;
	double e, ec;

	// TO DO: fix this with more reasonable general errors
	// Should probably make this an option for users, so need to update
	e = log(0.01 / 3.0);	 // Error for parts of sequence not covered by errorProfile, 0.01 probability of error and equal across bases?
	ec = log(1.0000 - 0.01); // Correct probability 0.99

	// readlike = malloc(numquery*(sizeof(double**)));

	// All reads get general error (parallelized: each read's slot is independent)
	#pragma omp parallel for schedule(dynamic, 1)
	for (i = 0; i < numquery; i++)
	{
		int treeNum = treeAssign[i];
		int treeOrder = readOrder[i];
		int j, k;

		readLikeTemp[treeNum][treeOrder] = (double **)malloc(readLengthTemp[treeNum][treeOrder] * (sizeof(double *)));
		// readlike[i] = malloc(readlength[i]*(sizeof(double*)));
		for (j = 0; j < readLengthTemp[treeNum][treeOrder]; j++)
		{
			if (readsTreeSorted[treeNum][treeOrder][j] != -1)
			{
				readLikeTemp[treeNum][treeOrder][j] = malloc(4 * (sizeof(double)));
				for (k = 0; k < 4; k++)
				{
					readLikeTemp[treeNum][treeOrder][j][k] =
						(k == readsTreeSorted[treeNum][treeOrder][j]) ? ec : e;
				}
			}
			else
			{
				// readlike[i][j][k] = e + errorTestLog;
				readLikeTemp[treeNum][treeOrder][j] = NULL;
			}
		}
	}

	// Only rough implementation of new read preprocessing, double check when doing errors again
	// fix error based on error profile
	// Should be more flexible to custom error profiles even if not the most efficient
	// Use fgets + strtol/strtod instead of fscanf: avoids per-call format string parsing overhead
	{
		char linebuf[256];
		while (fgets(linebuf, sizeof(linebuf), infile) != NULL)
		{
			char *p = linebuf, *next;
			int read_id = (int)strtol(p, &next, 10);
			if (next == p) continue;  // blank or malformed line
			p = next;
			int pos_id  = (int)strtol(p, &next, 10);  p = next;
			double err0 = strtod(p, &next);            p = next;
			double err1 = strtod(p, &next);            p = next;
			double err2 = strtod(p, &next);            p = next;
			double err3 = strtod(p, &next);            p = next;

			int treeNum   = treeAssign[read_id];
			int treeOrder = readOrder[read_id];
			if (readsTreeSorted[treeNum][treeOrder][pos_id] != -1)
			{
				readLikeTemp[treeNum][treeOrder][pos_id][0] = err0;
				readLikeTemp[treeNum][treeOrder][pos_id][1] = err1;
				readLikeTemp[treeNum][treeOrder][pos_id][2] = err2;
				readLikeTemp[treeNum][treeOrder][pos_id][3] = err3;
			}
			else
			{
				printf("Warning, error profile includes positions not in query alignment. Please review error profile, but ratePlacer is proceeding. For read %d pos %d\n", read_id, pos_id);
			}
		}
	}

	//for (i = 0; i < numquery; i++)
	//{
	//	treeNum = treeAssign[i];
	//	treeOrder = readOrder[i];
	//
	//	for (j = 0; j < readLengthTemp[treeNum][treeOrder]; j++){
	//		for (k=0; k<4; k++)
	//			printf("%lf,", readLikeTemp[treeNum][treeOrder][j][k]);
	//		printf("\n");
	//	}
	//	printf("\n");
	//}
}

double reassign_singleReadAge(double alpha, double parameters[7])
{
	// TO DO TO DO TO DO!!!!
	int nfun;
	double invector[3], lowbound[3], upbound[3], eh0 = 3e-8, age_like = 0.0, ageIncr, L2;
	double p[4];

	onDindic = 1;

	p[0] = parameters[0]; // read number
	p[1] = parameters[1]; // tree number
	p[2] = parameters[2]; // assignment
	p[3] = alpha;

	// printf("\t\tAlpha %.16f\n", p[3]);

	nfun = 0;
	invector[1] = bls[(int)p[1]][(int)p[2]] / 2.0;
	lowbound[1] = eh0;
	upbound[1] = bls[(int)p[1]][(int)p[2]] - eh0;

	L2 = Brent1D(invector, lowbound, upbound, 1, getlike_gamma_root_in_trifurcation_testAge_reassign, p, 4);

	return (L2);
}

// L is the node we are searching from, L_lik is the associated node
void greedyDown(int extra_data[3], int *L, double *L_lik, int root)
{
	// p[0]: readNum
	// p[1]: treeNum
	// p[2]: curNode
	int nfun, testNodes[3];
	double invector[2], testLik, testLik2;

	// L and L_lik are the curNode and what we are comparing too
	// Get children - childNodes[3]
	getGFLChildren(extra_data[2], testNodes, extra_data[1]);
	// Run on children
	// Test two children
	// Child 1
	extra_data[2] = testNodes[1];
	invector[0] = 0.5;
	invector[1] = bls[extra_data[1]][testNodes[1]] / 2.0;
	// testLik = findmax_amoeba_rand(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign, p, 3);
	// testLik = findmax_amoeba_rand_trans(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign_transform, p, 3);
	// testLik = GoldenSection(invector, lowbound, upbound, 1, reassign_singleReadAge, p, 3);
	testLik = minimize_brent(invector, 2, getlike_gamma_root_in_trifucation_single_read_brent_reassign, 100, extra_data);

	// printf("\t\t\tGreedy down L (%d) L_lik %.16f with child 1 (%d) like %.16f age %.16f and ", *L, *L_lik, testNodes[1], testLik, (1.0-invector[1])*(nodeages[(int)p[1]][(int)p[2]]+invector[2]));
	// printf("\t\t\tGreedy down L (%d) L_lik %.16f with child 1 (%d) like %.16f and ", *L, *L_lik, testNodes[1], testLik);

	// Child 2
	extra_data[2] = testNodes[2];
	invector[0] = 0.5;
	invector[1] = bls[extra_data[1]][testNodes[2]] / 2.0;
	// testLik2 = findmax_amoeba_rand(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign, p, 3);
	// testLik2 = findmax_amoeba_rand_trans(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign_transform, p, 3);
	// testLik2 = GoldenSection(invector, lowbound, upbound, 1, reassign_singleReadAge, p, 3);
	testLik2 = minimize_brent(invector, 2, getlike_gamma_root_in_trifucation_single_read_brent_reassign, 100, extra_data);

	// printf("child 2 (%d) like %.16f age %.16f\n", testNodes[2], testLik2, (1.0-invector[1])*(nodeages[(int)p[1]][(int)p[2]]+invector[2]));
	// printf("child 2 (%d) like %.16f\n", testNodes[2], testLik2);

	// if L_lik best, return
	// if child r, greedyDown with new root of child r if not leaf. if leaf, return with child r
	// if child l, greedyDown with new root of child l if not leaf. if leaf, return with child l

	if (testLik < *L_lik && testLik < testLik2)
	// if (*L_lik - testLik > 1.0 && testLik < testLik2)
	{
		// Child 1 best
		if (testNodes[1] < numseq[extra_data[1]])
		{
			// Leaf node, no more search
			*L = testNodes[1];
			*L_lik = testLik;
		}
		else
		{
			*L = testNodes[1];
			*L_lik = testLik;
			extra_data[2] = testNodes[1];
			greedyDown(extra_data, L, L_lik, root);
		}
	}
	else if (testLik2 < *L_lik && testLik2 < testLik)
	// else if (*L_lik - testLik2 > 1.0 && testLik2 < testLik)
	{
		// Child 2 best
		if (testNodes[2] < numseq[extra_data[1]])
		{
			// Leaf node, no more search
			*L = testNodes[2];
			*L_lik = testLik;
		}
		else
		{
			*L = testNodes[2];
			*L_lik = testLik;
			extra_data[2] = testNodes[2];
			greedyDown(extra_data, L, L_lik, root);
		}
	}
	// Implicit L_lik best, returns without exploring more
}

// L is the node we are searching from, L_lik is the associated node
void greedyUp(int extra_data[3], int *L, double *L_lik, int root)
{
	// p[0]: readNum
	// p[1]: treeNum
	// p[2]: curNode
	int nfun, testNodes[3];
	double invector[2], testLik, testLik2;

	// L and L_lik are the curNode and what we are comparing too
	// Get sibling and parent of "root" - parSib[3], 1 is parent and 2 is sibling
	getGFLParSib(extra_data[2], testNodes, (int)extra_data[1]);

	// Run on sibling and parent (if parent not past root, need to test for that)
	//  Parent
	if (testNodes[1] == root)
	{
		// Can't test, use max likelihood value
		testLik2 = 1000000000.0;
	}
	else
	{
		extra_data[2] = testNodes[1];
		invector[0] = 0.5;
		invector[1] = bls[extra_data[1]][testNodes[1]] / 2.0;
		// testLik2 = findmax_amoeba_rand(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign, p, 3);
		// testLik2 = findmax_amoeba_rand_trans(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign_transform, p, 3);
		// testLik2 = GoldenSection(invector, lowbound, upbound, 1, reassign_singleReadAge, p, 3);
		testLik2 = minimize_brent(invector, 2, getlike_gamma_root_in_trifucation_single_read_brent_reassign, 100, extra_data);
	}

	// printf("\t\t\tGreedy up L (%d) L_lik %.16f with parent (%d) like %.16f age %.16f ", *L, *L_lik, testNodes[1], testLik2, (1.0-invector[1])*(nodeages[(int)p[1]][(int)p[2]]+invector[2]));
	// printf("\t\t\tGreedy up L (%d) L_lik %.16f with parent (%d) like %.16f ", *L, *L_lik, testNodes[1], testLik2);

	// Sibling
	extra_data[2] = testNodes[2];
	invector[0] = 0.5;
	invector[1] = bls[extra_data[1]][testNodes[2]] / 2.0;
	// testLik = findmax_amoeba_rand(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign, p, 3);
	// testLik = findmax_amoeba_rand_trans(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign_transform, p, 3);
	// testLik = GoldenSection(invector, lowbound, upbound, 1, reassign_singleReadAge, p, 3);
	testLik = minimize_brent(invector, 2, getlike_gamma_root_in_trifucation_single_read_brent_reassign, 100, extra_data);

	// printf("and sibling (%d) like %.16f age %.16f\n", testNodes[2], testLik, (1.0-invector[1])*(nodeages[(int)p[1]][(int)p[2]]+invector[2]));
	// printf("and sibling (%d) like %.16f\n", testNodes[2], testLik);

	// if L_lik best, return
	// if sibling, greedyDown with new root of sibling if not leaf. if leaf, return with sibling
	// if parent, greedyUp with new root of parent
	if (testLik < *L_lik && testLik < testLik2)
	// if (*L_lik - testLik > 1.0 && testLik < testLik2)
	{
		// Sibling is best
		if (testNodes[2] < numseq[extra_data[1]])
		{
			// Leaf node, no more search
			*L = testNodes[2];
			*L_lik = testLik;
		}
		else
		{
			*L = testNodes[2];
			*L_lik = testLik;
			extra_data[2] = testNodes[2];
			greedyDown(extra_data, L, L_lik, root);
		}
	}
	else if (testLik2 < *L_lik && testLik2 < testLik)
	// else if (*L_lik - testLik2 > 1.0 && testLik2 < testLik)
	{
		// Parent is best, already tested its not root
		*L = testNodes[1];
		*L_lik = testLik2;
		extra_data[2] = testNodes[1];
		greedyUp(extra_data, L, L_lik, root);
	}
	// Implicit L was best, returns without exploring more
}


// Tronko node assignment testing
// Do not need to run this on leaf nodes
// For internal nodes, test all three edges to the node
// For root, test two edges to the root
void tronkoAssignmentTesting(int root, unsigned long int treeNum)
{
	int i, L1, surNodes[3], originalNode, extra_data[3];

	double L1_lik, testLik, testLik2, invector[2];

	printf("Tronko assignment testing\n");

	// loops for only the number of reads assigned to the tree
	// Test out .1 cutoff for reassignment
#pragma omp parallel for schedule(dynamic, 1) \
	private(surNodes, originalNode, extra_data, L1_lik, testLik, testLik2, invector)
	for (i = 0; i < usedTrees[treeNum]; i++)
	{
		onDindic = 1;
		// printf("Sequence %d of tree %d and node %d\n", i, treeNum, tempAssignments[treeNum][i]);
		// GET INITIAL BRANCH AND PLACEMENT ESTIMATES

		// Skip leaf nodes as there is no testing to do
		if (tempAssignments[treeNum][i] < numseq[treeNum])
		{
			//printf("\tLeaf node, skipping\n");
			continue;
		}

		extra_data[0] = i;
		extra_data[1] = treeNum;
		originalNode = tempAssignments[treeNum][i];
		//L1_lik = INFINITY;

		// Tronko assignment is the root
		if (originalNode == root)
		{
			//printf("\tRoot\n");
			// Get two children
			getGFLChildren(originalNode, surNodes, treeNum);

			// Test two children
			// Child 1
			extra_data[2] = surNodes[1];
			invector[0] = 0.5;
			invector[1] = bls[treeNum][surNodes[1]] / 2.0;
			L1_lik = minimize_brent(invector, 2, getlike_gamma_root_in_trifucation_single_read_brent_reassign, 100, extra_data);

			//if (invector[0] < 3e-8)
			//{
			//	invector[0] = 3e-8;
			//}
			//else if (invector[0] > 1.0 - 3e-8)
			//{
			//	invector[0] = 1.0 - 3e-8;
			//}

			//if (invector[1] < 3e-8)
			//{
			//	invector[1] = 3e-8;
			//}
			//else if (invector[1] > bls[treeNum][surNodes[1]] - 3e-8)
			//{
			//	invector[1] = bls[treeNum][surNodes[1]] - 3e-8;
			//}

			//printf("\tChild 1 (%d) like %.16f age %.16f (%.16f %.16f)\n", surNodes[1], L1_lik, (1.0-invector[0])*(nodeages[treeNum][surNodes[1]]+invector[1]), invector[0], invector[1]);
			//printf("\tChild 1 (%d) like %.16f and ", surNodes[1], L1_lik);

			// Child 2
			extra_data[2] = surNodes[2];
			invector[0] = 0.5;
			invector[1] = bls[treeNum][surNodes[2]] / 2.0;
			testLik = minimize_brent(invector, 2, getlike_gamma_root_in_trifucation_single_read_brent_reassign, 100, extra_data);

			//if (invector[0] < 3e-8)
			//{
			//	invector[0] = 3e-8;
			//}
			//else if (invector[0] > 1.0 - 3e-8)
			//{
			//	invector[0] = 1.0 - 3e-8;
			//}

			//if (invector[1] < 3e-8)
			//{
			//	invector[1] = 3e-8;
			//}
			//else if (invector[1] > bls[treeNum][surNodes[2]] - 3e-8)
			//{
			//	invector[1] = bls[treeNum][surNodes[2]] - 3e-8;
			//}

			//printf("\tChild 2 (%d) like %.16f age %.16f (%.16f %.16f)\n", surNodes[2], testLik, (1.0-invector[0])*(nodeages[treeNum][surNodes[2]]+invector[1]), invector[0], invector[1]);
			//printf("\tChild 2 (%d) like %.16f\n", surNodes[2], L1_lik);

			// Pick max, must pick one as ratePlacer will not work from the root node
			// Child 1 better (remember these are -loglik outputs)
			if (L1_lik < testLik)
			{
				tempAssignments[treeNum][i] = surNodes[1];
				#pragma omp atomic
				numReadsPerAssign[treeNum][surNodes[1]]++;
			}
			else
			{
				// Child 2 better
				tempAssignments[treeNum][i] = surNodes[2];
				#pragma omp atomic
				numReadsPerAssign[treeNum][surNodes[2]]++;
			}

			// remove an assignment count and then add to new
			#pragma omp atomic
			numReadsPerAssign[treeNum][originalNode]--;
		}
		else // Internal node
		{
			//printf("\tInternal\n");
			// Get two children as we place on branch above node
			getGFLChildren(originalNode, surNodes, treeNum);
			// Test
			// Test current and two children
			//  Original Assignment
			extra_data[2] = originalNode;
			invector[0] = 0.5;
			invector[1] = bls[treeNum][originalNode] / 2.0;
			L1_lik = minimize_brent(invector, 2, getlike_gamma_root_in_trifucation_single_read_brent_reassign, 100, extra_data);

			//if (invector[0] < 3e-8)
			//{
			//	invector[0] = 3e-8;
			//}
			//else if (invector[0] > 1.0 - 3e-8)
			//{
			//	invector[0] = 1.0 - 3e-8;
			//}

			//if (invector[1] < 3e-8)
			//{
			//	invector[1] = 3e-8;
			//}
			//else if (invector[1] > bls[treeNum][originalNode] - 3e-8)
			//{
			//	invector[1] = bls[treeNum][originalNode] - 3e-8;
			//}

			//printf("\tOriginal (%d) like %.16f age %.16f (%.16f %.16f)\n", originalNode, L1_lik, (1.0-invector[0])*(nodeages[treeNum][originalNode]+invector[1]), invector[0], invector[1]);
			//printf("\tOriginal (%d) like %.16f ", originalNode, L1_lik);

			// Child 1
			extra_data[2] = surNodes[1];
			invector[0] = 0.5;
			invector[1] = bls[treeNum][surNodes[1]] / 2.0;
			testLik = minimize_brent(invector, 2, getlike_gamma_root_in_trifucation_single_read_brent_reassign, 100, extra_data);

			//if (invector[0] < 3e-8)
			//{
			//	invector[0] = 3e-8;
			//}
			//else if (invector[0] > 1.0 - 3e-8)
			//{
			//	invector[0] = 1.0 - 3e-8;
			//}

			//if (invector[1] < 3e-8)
			//{
			//	invector[1] = 3e-8;
			//}
			//else if (invector[1] > bls[treeNum][surNodes[1]] - 3e-8)
			//{
			//	invector[1] = bls[treeNum][surNodes[1]] - 3e-8;
			//}

			//printf("\tChild 1 (%d) like %.16f age %.16f (%.16f %.16f)\n", surNodes[1], testLik, (1.0-invector[0])*(nodeages[treeNum][surNodes[1]]+invector[1]), invector[0], invector[1]);
			//printf("\tchild 1 (%d) like %.16f and ", surNodes[1], testLik);

			// Child 2
			extra_data[2] = surNodes[2];
			invector[0] = 0.5;
			invector[1] = bls[treeNum][surNodes[2]] / 2.0;
			testLik2 = minimize_brent(invector, 2, getlike_gamma_root_in_trifucation_single_read_brent_reassign, 100, extra_data);

			//if (invector[0] < 3e-8)
			//{
			//	invector[0] = 3e-8;
			//}
			//else if (invector[0] > 1.0 - 3e-8)
			//{
			//	invector[0] = 1.0 - 3e-8;
			//}

			//if (invector[1] < 3e-8)
			//{
			//	invector[1] = 3e-8;
			//}
			//else if (invector[1] > bls[treeNum][surNodes[2]] - 3e-8)
			//{
			//	invector[1] = bls[treeNum][surNodes[2]] - 3e-8;
			//}

			//printf("\tChild 2 (%d) like %.16f age %.16f (%.16f %.16f)\n", surNodes[2], testLik2, (1.0-invector[0])*(nodeages[treeNum][surNodes[2]]+invector[1]), invector[0], invector[1]);
			//printf("\tchild 2 (%d) like %.16f\n", surNodes[2], testLik2);

			// If original node is best, do nothing
			// If child 1 is best, assign to child 1
			// If child 2 is best, assign to child 2
			if (testLik < testLik2 && testLik < L1_lik)
			{
				tempAssignments[treeNum][i] = surNodes[1];
				#pragma omp atomic
				numReadsPerAssign[treeNum][surNodes[1]]++;
				#pragma omp atomic
				numReadsPerAssign[treeNum][originalNode]--;
			}
			else if (testLik2 < testLik && testLik2 < L1_lik)
			{
				tempAssignments[treeNum][i] = surNodes[2];
				#pragma omp atomic
				numReadsPerAssign[treeNum][surNodes[2]]++;
				#pragma omp atomic
				numReadsPerAssign[treeNum][originalNode]--;
			}
		}

		//printf("\tFinal assignment %d\n", tempAssignments[treeNum][i]);
	}

}

void bestAssignment(int root, unsigned long int treeNum)
{
	int i, L1, nfun, testNode, surNodes[3], extra_data[3];
	double L1_lik, testLik, testLik2, invector[2];

	printf("Testing assignments\n");

	// loops for only the number of reads assigned to the tree
	// Test out .1 cutoff for reassignment
#pragma omp parallel for schedule(dynamic, 1) \
	private(L1, nfun, testNode, surNodes, extra_data, L1_lik, testLik, testLik2, invector)
	for (i = 0; i < usedTrees[treeNum]; i++)
	{
		onDindic = 1;
		// printf("Sequence %d of tree %d and node %d\n", i, treeNum, tempAssignments[treeNum][i]);
		// GET INITIAL BRANCH AND PLACEMENT ESTIMATES

		extra_data[0] = i;		// read num does not change
		extra_data[1] = treeNum; // tree num does not change
		extra_data[2] = tempAssignments[treeNum][i];
		L1 = tempAssignments[treeNum][i];
		L1_lik = INFINITY;

		// Tronko assignment is the root
		if (extra_data[2] == root)
		{
			// printf("\tRoot\n");
			// Get two children
			getGFLChildren(extra_data[2], surNodes, treeNum);

			// Test two children
			// Child 1
			testNode = surNodes[1];
			extra_data[2] = testNode;
			invector[0] = 0.5;
			invector[1] = bls[treeNum][testNode] / 2.0;
			// L1_lik = findmax_amoeba_rand(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign, p, 3);
			// L1_lik = findmax_amoeba_rand_trans(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign_transform, p, 3);
			//L1_lik = GoldenSection(invector, lowbound, upbound, 1, reassign_singleReadAge, p, 3);
			L1_lik = minimize_brent(invector, 2, getlike_gamma_root_in_trifucation_single_read_brent_reassign, 100, extra_data);

			// printf("\t\tChild 1 (%d) like %.16f age %.16f and ", surNodes[1], L1_lik, (1.0-invector[1])*(nodeages[treeNum][testNode]+invector[2]));
			// printf("\t\tChild 1 (%d) like %.16f and ", surNodes[1], L1_lik);

			// Child 2
			testNode = surNodes[2];
			extra_data[2] = testNode;
			invector[0] = 0.5;
			invector[1] = bls[treeNum][testNode] / 2.0;
			// testLik = findmax_amoeba_rand(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign, p, 3);
			// testLik = findmax_amoeba_rand_trans(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign_transform, p, 3);
			//testLik = GoldenSection(invector, lowbound, upbound, 1, reassign_singleReadAge, p, 3);
			testLik = minimize_brent(invector, 2, getlike_gamma_root_in_trifucation_single_read_brent_reassign, 100, extra_data);

			// Pick max
			// Traverse down based on max
			// Child 1 better (remember these are -loglik outputs)

			// printf("child 2 (%d) like %.16f age %.16f\n",  surNodes[2], testLik, (1.0-invector[1])*(nodeages[treeNum][testNode]+invector[2]));
			// printf("child 2 (%d) like %.16f\n",  surNodes[2], testLik);

			if (L1_lik < testLik)
			{
				// Search down child 1
				if (surNodes[1] < numseq[treeNum])
				{
					// Leaf node, no more search
					L1 = surNodes[1];
				}
				else
				{
					L1 = surNodes[1];
					extra_data[2] = surNodes[1];
					greedyDown(extra_data, &L1, &L1_lik, root);
				}
			}
			else
			{
				// Child 2 better, search down child 2
				if (surNodes[2] < numseq[treeNum])
				{
					// Leaf node, no more search
					L1 = surNodes[2];
				}
				else
				{
					L1 = surNodes[2];
					L1_lik = testLik;
					extra_data[2] = surNodes[2];
					greedyDown(extra_data, &L1, &L1_lik, root);
				}
			}

			// remove an assignment count and then add to new
			#pragma omp atomic
			numReadsPerAssign[treeNum][tempAssignments[treeNum][i]]--;
			#pragma omp atomic
			numReadsPerAssign[treeNum][L1]++;
			tempAssignments[treeNum][i] = L1;
			// assignAges[i] =  nodeages[treeNum][L1] + bls[treeNum][L1]; //done after merging reads now
		}
		// Tronko assignment is a leaf
		else if (extra_data[2] < numseq[treeNum])
		{
			// printf("\tLeaf\n");
			// Get current branch, parent, and sibling
			getGFLParSib(extra_data[2], surNodes, treeNum);
			// Test
			//  Original Assignment
			testNode = extra_data[2];
			invector[0] = 0.5;
			invector[1] = bls[treeNum][testNode] / 2.0;
			// L1_lik = findmax_amoeba_rand(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign, p, 3);
			// L1_lik = findmax_amoeba_rand_trans(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign_transform, p, 3);
			//L1_lik = GoldenSection(invector, lowbound, upbound, 1, reassign_singleReadAge, p, 3);
			L1_lik = minimize_brent(invector, 2, getlike_gamma_root_in_trifucation_single_read_brent_reassign, 100, extra_data);

			// printf("\t\tOriginal (%d) like %.16f age %.16f ", L1, L1_lik, (1.0-invector[1])*(nodeages[treeNum][testNode]+invector[2]));
			// printf("\t\tOriginal (%d) like %.16f ", L1, L1_lik);

			// Sibling
			testNode = surNodes[2];
			extra_data[2] = testNode;
			invector[0] = 0.5;
			invector[1] = bls[treeNum][testNode] / 2.0;
			// testLik = findmax_amoeba_rand(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign, p, 3);
			// testLik = findmax_amoeba_rand_trans(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign_transform, p, 3);
			//testLik = GoldenSection(invector, lowbound, upbound, 1, reassign_singleReadAge, p, 3);
			testLik = minimize_brent(invector, 2, getlike_gamma_root_in_trifucation_single_read_brent_reassign, 100, extra_data);

			// printf("sibling (%d) like %.16f age %.16f and ", surNodes[2], testLik, (1.0-invector[1])*(nodeages[treeNum][testNode]+invector[2]));
			// printf("sibling (%d) like %.16f and ", surNodes[2], testLik);

			// Parent
			if (surNodes[1] == root)
			{
				// Can't test, use max likelihood value
				testLik2 = 1000000000.0;
			}
			else
			{
				testNode = surNodes[1];
				extra_data[2] = testNode;
				invector[0] = 0.5;
				invector[1] = bls[treeNum][testNode] / 2.0;
				// testLik2 = findmax_amoeba_rand(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign, p, 3);
				// testLik2 = findmax_amoeba_rand_trans(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign_transform, p, 3);
				//testLik2 = GoldenSection(invector, lowbound, upbound, 1, reassign_singleReadAge, p, 3);
				testLik2 = minimize_brent(invector, 2, getlike_gamma_root_in_trifucation_single_read_brent_reassign, 100, extra_data);
			}

			// printf("parent (%d) like %.16f age %.16f\n", surNodes[1], testLik2, (1.0-invector[1])*(nodeages[treeNum][testNode]+invector[2]));
			// printf("parent (%d) like %.16f\n", surNodes[1], testLik2);

			// If current, end
			// if sibling > parent, go down
			// if parent > sibling, go up
			// if(testLik2 < L1_lik)
			//{
			//	L1 = surNodes[1];
			//	L1_lik = testLik2;
			//	p[2] = surNodes[1];
			//	greedyUp(p, &L1, &L1_lik, root);
			//
			//	numReadsPerAssign[treeNum][tempAssignments[treeNum][i]]--;
			//	numReadsPerAssign[treeNum][L1]++;
			//	tempAssignments[treeNum][i] = L1;
			// }
			if (testLik < L1_lik || testLik2 < L1_lik)
			// if (L1_lik - testLik > 1.0 || L1_lik - testLik2 > 1.0)
			{
				if (testLik < testLik2)
				{
					// Sibling is best
					if (surNodes[2] < numseq[treeNum])
					{
						// Leaf node, no more search
						L1 = surNodes[2];
					}
					else
					{
						L1 = surNodes[2];
						L1_lik = testLik;
						extra_data[2] = surNodes[2];
						greedyDown(extra_data, &L1, &L1_lik, root);
					}
				}
				else
				{
					// Parent is best, already tested its not root
					L1 = surNodes[1];
					L1_lik = testLik2;
					extra_data[2] = surNodes[1];
					greedyUp(extra_data, &L1, &L1_lik, root);
				}
				// remove an assignment count and then add to new
				#pragma omp atomic
				numReadsPerAssign[treeNum][tempAssignments[treeNum][i]]--;
				#pragma omp atomic
				numReadsPerAssign[treeNum][L1]++;
				tempAssignments[treeNum][i] = L1;
				// assignAges[i] =  nodeages[treeNum][L1] + bls[treeNum][L1];
			}
			// Implicit L1 is the best, we don't need to do anything
		}
		else
		{
			// printf("\tInternal\n");
			// Get two children as we place on branch above node
			getGFLChildren(extra_data[2], surNodes, treeNum);
			// Test
			// Test current and two children
			//  Original Assignment
			testNode = extra_data[2];
			invector[0] = 0.5;
			invector[1] = bls[treeNum][testNode] / 2.0;
			// L1_lik = findmax_amoeba_rand(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign, p, 3);
			// L1_lik = findmax_amoeba_rand_trans(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign_transform, p, 3);
			//L1_lik = GoldenSection(invector, lowbound, upbound, 1, reassign_singleReadAge, p, 3);
			L1_lik = minimize_brent(invector, 2, getlike_gamma_root_in_trifucation_single_read_brent_reassign, 100, extra_data);

			// printf("\t\tOriginal (%d) like %.16f age %.16f ", L1, L1_lik, (1.0-invector[1])*(nodeages[treeNum][testNode]+invector[2]));
			// printf("\t\tOriginal (%d) like %.16f ", L1, L1_lik);

			// Child 1
			testNode = surNodes[1];
			extra_data[2] = testNode;
			invector[0] = 0.5;
			invector[1] = bls[treeNum][testNode] / 2.0;
			// testLik = findmax_amoeba_rand(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign, p, 3);
			// testLik = findmax_amoeba_rand_trans(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign_transform, p, 3);
			// testLik = GoldenSection(invector, lowbound, upbound, 1, reassign_singleReadAge, p, 3);
			testLik = minimize_brent(invector, 2, getlike_gamma_root_in_trifucation_single_read_brent_reassign, 100, extra_data);

			// printf("child 1 (%d) like %.16f age %.16f and ", surNodes[1], testLik, (1.0-invector[1])*(nodeages[treeNum][testNode]+invector[2]));
			// printf("child 1 (%d) like %.16f and ", surNodes[1], testLik);

			// Child 2
			testNode = surNodes[2];
			extra_data[2] = testNode;
			invector[0] = 0.5;
			invector[1] = bls[treeNum][testNode] / 2.0;
			// testLik2 = findmax_amoeba_rand(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign, p, 3);
			// testLik2 = findmax_amoeba_rand_trans(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_reassign_transform, p, 3);
			// testLik2 = GoldenSection(invector, lowbound, upbound, 1, reassign_singleReadAge, p, 3);
			testLik2 = minimize_brent(invector, 2, getlike_gamma_root_in_trifucation_single_read_brent_reassign, 100, extra_data);

			// printf("child 2 (%d) like %.16f age %.16f\n", surNodes[2], testLik2, (1.0-invector[1])*(nodeages[treeNum][testNode]+invector[2]));
			// printf("child 2 (%d) like %.16f\n", surNodes[2], testLik2);

			// If parent > children, go up
			// If children > parent, go down depending on which child was maximum
			if (L1_lik < testLik && L1_lik < testLik2)
			{
				// Original assignment best, go up
				L1 = tempAssignments[treeNum][i];
				extra_data[2] = tempAssignments[treeNum][i];
				greedyUp(extra_data, &L1, &L1_lik, root);
			}
			else if (testLik < testLik2 && testLik < L1_lik)
			// else if (L1_lik - testLik > 1.0 && testLik < testLik2)
			{
				// Child 1 is best
				if (surNodes[1] < numseq[treeNum])
				{
					// Leaf node, no more search
					L1 = surNodes[1];
				}
				else
				{
					L1 = surNodes[1];
					L1_lik = testLik;
					extra_data[2] = surNodes[1];
					greedyDown(extra_data, &L1, &L1_lik, root);
				}
			}
			else
			// else if (L1_lik - testLik2 > 1.0 && testLik2 < testLik)
			{
				// Child 2 is best
				if (surNodes[2] < numseq[treeNum])
				{
					// Leaf node, no more search
					L1 = surNodes[2];
				}
				else
				{
					L1 = surNodes[2];
					L1_lik = testLik2;
					extra_data[2] = surNodes[2];
					greedyDown(extra_data, &L1, &L1_lik, root);
				}
			}

			// remove an assignment count and then add to new
			#pragma omp atomic
			numReadsPerAssign[treeNum][tempAssignments[treeNum][i]]--;
			#pragma omp atomic
			numReadsPerAssign[treeNum][L1]++;
			tempAssignments[treeNum][i] = L1;
			// assignAges[i] =  nodeages[treeNum][L1] + bls[treeNum][L1];
		}

		// GET ESTIMATES ABOUT BRANCH AND PLACEMENT HERE FOR FINAL CHOICE
		// printf("\tNow node %d\n", tempAssignments[treeNum][i]);
		// printf("\tBest assignment of %d\n", tempAssignments[treeNum][i]);
	}
}

// This function merges all the reads that have been assigned to the same node
void mergeReads(unsigned long int treeNum, unsigned long int refBases)
{
	printf("Merging Reads\n");
	unsigned long int numNodes = 2 * numseq[treeNum] - 1;
	int **perEdgeReads = (int **)malloc(numNodes * sizeof(int *));
	int *perEdgeIndex = (int *)calloc(numNodes, sizeof(int));
	unsigned long int numSeqs = 0, index;

	for (unsigned long int i = 0; i < numNodes; i++)
	{
		if (numReadsPerAssign[treeNum][i] != 0)
		{
			//printf("numReadsPerAssign[treeNum][%d]: %d\n", i, numReadsPerAssign[treeNum][i]);
			perEdgeReads[i] = (int *)malloc(numReadsPerAssign[treeNum][i] * sizeof(int));
			numSeqs++;
		}
	}

	// can I change this?
	index = tempnumquery;
	tempnumquery += numSeqs;

	// printf("%d\n", tempnumquery);

	// Not safe, should be testing with temp pointer for NULL to be able to free...
	int** tempptr = (int **)realloc(QUERYDATA, tempnumquery * (sizeof(int *)));
	if (tempptr != NULL)
	{
		QUERYDATA = tempptr;
	} else {
		printf("Failed to reallocate QUERYDATA in mergeReads\n");
		exit(0);
	}

	unsigned long int* tempptr2 = (unsigned long int *)realloc(startpos, tempnumquery * (sizeof(unsigned long int)));
	if (tempptr2 != NULL)
	{
		startpos = tempptr2;
	} else {
		printf("Failed to reallocate startpos in mergeReads\n");
		exit(0);
	}

	tempptr2 = (unsigned long int *)realloc(readlength, tempnumquery * (sizeof(unsigned long int)));
	if (tempptr2 != NULL)
	{
		readlength = tempptr2;
	} else {
		printf("Failed to reallocate readlength in mergeReads\n");
		exit(0);
	}

	int *tempptr_int = (int *)realloc(assignments, tempnumquery * (sizeof(int)));
	if (tempptr_int != NULL)
	{
		assignments = tempptr_int;
	} else {
		printf("Failed to reallocate assignments in mergeReads\n");
		exit(0);
	}

	tempptr_int = (int *)realloc(treeAssign, tempnumquery * (sizeof(int)));
	if (tempptr_int != NULL)
	{
		treeAssign = tempptr_int;
	} else {
		printf("Failed to reallocate treeAssign in mergeReads\n");
		exit(0);
	}

	double* tempptr3 = (double *)realloc(assignAges, tempnumquery * (sizeof(double)));
	if (tempptr3 != NULL)
	{
		assignAges = tempptr3;
	} else {
		printf("Failed to reallocate assignAges in mergeReads\n");
		exit(0);
	}

	double*** tempptr4 = (double ***)realloc(readlike, tempnumquery * sizeof(double **));
	if (tempptr4 != NULL)
	{
		readlike = tempptr4;
	} else {
		printf("Failed to reallocate readlike in mergeReads\n");
		exit(0);
	}
	//QUERYDATA = (int **)realloc(QUERYDATA, tempnumquery * (sizeof(int *)));
	//startpos = realloc(startpos, tempnumquery * (sizeof(int)));
	//readlength = realloc(readlength, tempnumquery * (sizeof(int)));
	//assignAges = (double *)realloc(assignAges, tempnumquery * (sizeof(double)));
	//assignments = realloc(assignments, tempnumquery * (sizeof(int)));
	//treeAssign = realloc(treeAssign, tempnumquery * (sizeof(int)));
	//readlike = (double ***)realloc(readlike, tempnumquery * sizeof(double **));

	for (unsigned long int i = 0; i < usedTrees[treeNum]; i++)
	{
		perEdgeReads[tempAssignments[treeNum][i]][perEdgeIndex[tempAssignments[treeNum][i]]] = i;
		perEdgeIndex[tempAssignments[treeNum][i]]++;
	}

	// Per-node result storage for pass 1.
	// Pointers are NULL for nodes with no reads or that fail the coverage filter.
	int **node_querydata          = (int **)calloc(numNodes, sizeof(int *));
	double ***node_readlike       = (double ***)calloc(numNodes, sizeof(double **));
	unsigned long int *node_readlength   = (unsigned long int *)malloc(numNodes * sizeof(unsigned long int));
	unsigned long int *node_startpos     = (unsigned long int *)malloc(numNodes * sizeof(unsigned long int));
	double *node_assignage        = (double *)malloc(numNodes * sizeof(double));
	unsigned long int *node_refcoverage  = (unsigned long int *)malloc(numNodes * sizeof(unsigned long int));
	int *node_keep                = (int *)calloc(numNodes, sizeof(int));

	// Pass 1 (parallel): each node independently accumulates its reads, applies coverage
	// filter, builds consensus sequence and error likelihoods into private storage.
	// drand48_r() is used instead of drand48() to avoid a data race on the global RNG state.
#pragma omp parallel for schedule(dynamic)
	for (unsigned long int i = 0; i < numNodes; i++)
	{
		if (numReadsPerAssign[treeNum][i] == 0)
			continue;

		// Private baseCounts and errorSums for this node (calloc = already zeroed)
		int *baseCounts = (int *)calloc(refBases * 4, sizeof(int));
		double *errorSums = (double *)calloc(refBases * 4, sizeof(double));
		unsigned long int refCoverage = 0;

		for (unsigned long int j = 0; j < numReadsPerAssign[treeNum][i]; j++)
		{
			// count bases used perEdgeReads as read index in readsTreeSorted
			unsigned long int readIndex = perEdgeReads[i][j];
			unsigned long int refPos = startposTemp[treeNum][readIndex];

			for (unsigned long int pos = 0; pos < readLengthTemp[treeNum][readIndex]; pos++)
			{
				// baseCounts[(refPos + pos) * 4 + base] == baseCounts[refPos + pos][base]
				if (readsTreeSorted[treeNum][readIndex][pos] != -1)
				{
					errorSums[(refPos + pos) * 4]     += readLikeTemp[treeNum][readIndex][pos][0];
					errorSums[(refPos + pos) * 4 + 1] += readLikeTemp[treeNum][readIndex][pos][1];
					errorSums[(refPos + pos) * 4 + 2] += readLikeTemp[treeNum][readIndex][pos][2];
					errorSums[(refPos + pos) * 4 + 3] += readLikeTemp[treeNum][readIndex][pos][3];

					if (baseCounts[(refPos + pos) * 4] == 0 && baseCounts[(refPos + pos) * 4 + 1] == 0 && baseCounts[(refPos + pos) * 4 + 2] == 0 && baseCounts[(refPos + pos) * 4 + 3] == 0)
						refCoverage++;

					int base = readsTreeSorted[treeNum][readIndex][pos];
					baseCounts[(refPos + pos) * 4 + base]++;
				}
			}
		}

		// printf("Reference length of %d and coverage of %d (%lf) with %d reads\n", refBases, refCoverage, (double)refCoverage/(double)refBases, numReadsPerAssign[treeNum][i]);

		// UNCOMMENT THIS IN NORMAL VERSIONS< JUST FOR TESTING SOME CASES
		// TO DO: DISREGARD -/N FOR refBases (mayeb it is, I should check)
		int keepAssignment = 1; // Assume we keep the assignment initially
		//printf("merge_coverage_mode: %d, merge_coverage_threshold: %f\n", merge_coverage_mode, merge_coverage_threshold);
		if (merge_coverage_mode == MERGE_MODE_FRACTION) {
			if ((double)refCoverage / (double)refBases < merge_coverage_threshold) {
				keepAssignment = 0;
			}
		} else { // MERGE_MODE_BP
			if (refCoverage < (int)merge_coverage_threshold) { // Cast threshold to int for BP comparison
				keepAssignment = 0;
			}
		}

		if (!keepAssignment)
		{
			free(baseCounts);
			free(errorSums);
			continue; // node_keep[i] stays 0
		}

		// printf("Assignment %d has %d bases covered\n", i, refCoverage);

		unsigned long int firstPos = 0, lastPos = 0;

		for (unsigned long int pos = 0; pos < refBases; pos++)
		{
			if (baseCounts[pos * 4] > 0 || baseCounts[pos * 4 + 1] > 0 || baseCounts[pos * 4 + 2] > 0 || baseCounts[pos * 4 + 3] > 0)
			{
				firstPos = pos;
				break;
			}
		}

		lastPos = firstPos;

		for (unsigned long int pos = firstPos + 1; pos < refBases; pos++)
		{
			if (baseCounts[pos * 4] > 0 || baseCounts[pos * 4 + 1] > 0 || baseCounts[pos * 4 + 2] > 0 || baseCounts[pos * 4 + 3] > 0)
			{
				lastPos = pos;
			}
		}

		unsigned long int rdlen = lastPos - firstPos + 1;
		int *qdata = (int *)malloc(rdlen * sizeof(int));
		double **rlike = (double **)malloc(rdlen * sizeof(double *));

		// Per-node RNG state seeded from node index for reproducibility across parallelizations
		struct drand48_data rng_state;
		srand48_r((long int)(i + 1), &rng_state);

		for (unsigned long int pos = 0; pos < rdlen; pos++)
		{
			unsigned long int refPos = pos + firstPos;
			// no base info here
			if (baseCounts[refPos * 4] == 0 && baseCounts[refPos * 4 + 1] == 0 && baseCounts[refPos * 4 + 2] == 0 && baseCounts[refPos * 4 + 3] == 0)
			{
				qdata[pos] = -1;
				rlike[pos] = NULL;
			}
			else
			{
				rlike[pos] = (double *)malloc(4 * sizeof(double));
				rlike[pos][0] = errorSums[refPos * 4];
				rlike[pos][1] = errorSums[refPos * 4 + 1];
				rlike[pos][2] = errorSums[refPos * 4 + 2];
				rlike[pos][3] = errorSums[refPos * 4 + 3];

				//// Revisit when using error versions
				// baseMax = baseCounts[pos * 4];
				// baseIndex = 0;
				////is biased for later bases... what is best approach for ties?
				// for (int k = 1; k < 4; k++)
				//{
				//	if (baseCounts[refPos * 4 + k] > baseMax)
				//	{
				//		baseMax = baseCounts[refPos * 4 + k];
				//		baseIndex = k;
				//	}
				// }
				// qdata[pos] = baseIndex;

				double base_select;
				drand48_r(&rng_state, &base_select);
				// explicity caste everything to float
				// Test uniform
				double base_sum = (double)(baseCounts[refPos * 4] + baseCounts[refPos * 4 + 1] + baseCounts[refPos * 4 + 2] + baseCounts[refPos * 4 + 3]);
				double base_a = (double)baseCounts[refPos * 4] / base_sum;
				double base_c = (double)(baseCounts[refPos * 4] + baseCounts[refPos * 4 + 1]) / base_sum;
				double base_g = (double)(baseCounts[refPos * 4] + baseCounts[refPos * 4 + 1] + baseCounts[refPos * 4 + 2]) / base_sum;
				double base_t = (double)(baseCounts[refPos * 4] + baseCounts[refPos * 4 + 1] + baseCounts[refPos * 4 + 2] + baseCounts[refPos * 4 + 3]) / base_sum;

				// printf("Sum: %lf, A: %lf (%d), C: %lf (%d), G: %lf (%d), T: %lf (%d), rand: %lf\n", base_sum, base_a, baseCounts[refPos * 4], base_c, baseCounts[refPos * 4 + 1], base_g, baseCounts[refPos * 4 + 2], base_t, baseCounts[refPos * 4 + 3], base_select);

				if (base_t != 1.0 || base_select < 0.0 || base_select > 1.0)
				{
					printf("Error in selecting base: %lf, %lf, %lf, %lf, %lf\n", base_a, base_c, base_g, base_t, base_select);
					exit(0);
				}

				if (base_a > base_select)       qdata[pos] = 0;
				else if (base_c > base_select)  qdata[pos] = 1;
				else if (base_g > base_select)  qdata[pos] = 2;
				else if (base_t >= base_select) qdata[pos] = 3;
				else { printf("Error in base selection\n"); exit(0); }
			}
		}

		node_querydata[i]   = qdata;
		node_readlike[i]    = rlike;
		node_readlength[i]  = rdlen;
		node_startpos[i]    = firstPos;
		node_assignage[i]   = nodeages[treeNum][i] + bls[treeNum][i];
		node_refcoverage[i] = refCoverage;
		node_keep[i]        = 1;

		free(baseCounts);
		free(errorSums);
	}

	// Pass 2 (serial): assign sequential index values and transfer pointers to global arrays.
	// readsfile output happens here to preserve node ordering.
	static const char base_chars[] = "ACGT";
	for (unsigned long int i = 0; i < numNodes; i++)
	{
		if (numReadsPerAssign[treeNum][i] == 0)
			continue;

		if (!node_keep[i])
		{
			// do I need to do anything else if I abort?
			tempnumquery--;
			continue;
		}

		readlength[index] = node_readlength[i];
		startpos[index]   = node_startpos[i];
		QUERYDATA[index]  = node_querydata[i];
		readlike[index]   = node_readlike[i];
		treeAssign[index] = treeNum;
		assignments[index] = i;
		assignAges[index] = node_assignage[i];

		if (readsfile)
		{
			fprintf(readsfile, "Sequence %lu, assignment %lu of length %lu starting at %lu containing %d reads covering %lu:\n",
			        index, i, readlength[index], node_startpos[i], numReadsPerAssign[treeNum][i], node_refcoverage[i]);
			for (unsigned long int pos = 0; pos < readlength[index]; pos++)
			{
				if (QUERYDATA[index][pos] == -1) fprintf(readsfile, "-");
				else fprintf(readsfile, "%c", base_chars[QUERYDATA[index][pos]]);
			}
			fprintf(readsfile, "\n");
		}

		index++;
	}

	// Free per-node pointer arrays; pointed-to memory is now owned by global arrays
	free(node_querydata);
	free(node_readlike);
	free(node_readlength);
	free(node_startpos);
	free(node_assignage);
	free(node_refcoverage);
	free(node_keep);

	free(perEdgeIndex);
	for (unsigned long int i = 0; i < numNodes; i++)
	{
		if (numReadsPerAssign[treeNum][i] != 0) // Only free if it was allocated
		{
			free(perEdgeReads[i]);
		}
	}
	free(perEdgeReads);

	// TO DO: Free per tree temps
}

// Do not merge reads and do analysis separately
// void keepSeparateReads(int treeNum, int refBases)
//{
//	//numquery = tempnumquery;
//	// Transfer data to used pointers and make temp pointers null
//
//
//	QUERYDATA = readsTreeSorted;
//	readsTreeSorted = NULL;
//
//	readlength = readLengthTemp;
//	readLengthTemp = NULL;
//
//	startpos = startposTemp;
//	startposTemp = NULL;
//
//	assignments = tempAssignments;
//	tempAssignments = NULL;
//
//	readlike = readLikeTemp;
//	readLikeTemp = NULL;
//}

// Keep reads separate instead of merging them based on assignment
void keepSeparateReads(unsigned long int treeNum)
{
	unsigned long int numReadsInTree = usedTrees[treeNum];
	unsigned long int index;

	if (numReadsInTree == 0)
		return;

	// Calculate the new total number of query sequences
	// Note: This assumes keepSeparateReads is called sequentially for each tree
	index = tempnumquery; // Start indexing from the current total
	tempnumquery += numReadsInTree; // Increment total by the number of reads in this tree

	// Reallocate global arrays to accommodate the new reads
	// Not safe, should be testing with temp pointer for NULL to be able to free...
	int **tempptr = (int **)realloc(QUERYDATA, tempnumquery * (sizeof(int *)));
	if (tempptr != NULL)
	{
		QUERYDATA = tempptr;
	} else {
		printf("Failed to reallocate QUERYDATA in keepSeparateReads\n");
		exit(0);
	}

	unsigned long int *tempptr2 = (unsigned long int *)realloc(startpos, tempnumquery * (sizeof(unsigned long int)));
	if (tempptr2 != NULL)
	{
		startpos = tempptr2;
	} else {
		printf("Failed to reallocate startpos in keepSeparateReads\n");
		exit(0);
	}

	tempptr2 = (unsigned long int *)realloc(readlength, tempnumquery * (sizeof(unsigned long int)));
	if (tempptr2 != NULL)
	{
		readlength = tempptr2;
	} else {
		printf("Failed to reallocate readlength in keepSeparateReads\n");
		exit(0);
	}

	int *tempptr_int = (int *)realloc(assignments, tempnumquery * (sizeof(int)));
	if (tempptr_int != NULL)
	{
		assignments = tempptr_int;
	} else {
		printf("Failed to reallocate assignments in keepSeparateReads\n");
		exit(0);
	}

	tempptr_int = (int *)realloc(treeAssign, tempnumquery * (sizeof(int)));
	if (tempptr_int != NULL)
	{
		treeAssign = tempptr_int;
	} else {
		printf("Failed to reallocate treeAssign in keepSeparateReads\n");
		exit(0);
	}

	double *tempptr3 = (double *)realloc(assignAges, tempnumquery * (sizeof(double)));
	if (tempptr3 != NULL)
	{
		assignAges = tempptr3;
	} else {
		printf("Failed to reallocate assignAges in keepSeparateReads\n");
		exit(0);
	}

	double ***tempptr4 = (double ***)realloc(readlike, tempnumquery * sizeof(double **));
	if (tempptr4 != NULL)
	{
		readlike = tempptr4;
	} else {
		printf("Failed to reallocate readlike in keepSeparateReads\n");
		exit(0);
	}

	// Iterate through each read for the current tree and copy its data
	for (unsigned long int i = 0; i < numReadsInTree; i++)
	{
		unsigned long int currentReadLength = readLengthTemp[treeNum][i];
		unsigned long int currentEdgeAssignment = tempAssignments[treeNum][i];

		// Copy basic info
		startpos[index] = startposTemp[treeNum][i];
		readlength[index] = currentReadLength;
		assignments[index] = currentEdgeAssignment;
		treeAssign[index] = treeNum;
		assignAges[index] = nodeages[treeNum][currentEdgeAssignment] + bls[treeNum][currentEdgeAssignment];

		// Allocate and copy sequence data (QUERYDATA)
		QUERYDATA[index] = (int *)malloc(currentReadLength * sizeof(int));
		if (QUERYDATA[index] == NULL)
		{
			printf("Failed to allocate QUERYDATA[%lu] in keepSeparateReads\n", index);
			exit(0);
		}
		memcpy(QUERYDATA[index], readsTreeSorted[treeNum][i], currentReadLength * sizeof(int));

		// Allocate and copy likelihood data (readlike)
		readlike[index] = (double **)malloc(currentReadLength * sizeof(double *));
		if (readlike[index] == NULL)
		{
			printf("Failed to allocate readlike[%lu] in keepSeparateReads\n", index);
			exit(0);
		}
		for (unsigned long int pos = 0; pos < currentReadLength; pos++)
		{
			if (readLikeTemp[treeNum][i][pos] != NULL)
			{
				readlike[index][pos] = (double *)malloc(4 * sizeof(double));
				if (readlike[index][pos] == NULL)
				{
					printf("Failed to allocate readlike[%lu][%lu] in keepSeparateReads\n", index, pos);
					exit(0);
				}
				memcpy(readlike[index][pos], readLikeTemp[treeNum][i][pos], 4 * sizeof(double));
			}
			else
			{
				readlike[index][pos] = NULL; // Keep NULL if original was NULL
			}
		}

		index++; // Move to the next global index
	}

	// Note: Freeing of the temporary per-tree arrays (readsTreeSorted[treeNum], readLengthTemp[treeNum], etc.)
	// should ideally happen after all trees have been processed, likely in freeData or similar.
}


// numseq: number of sequences in reference data set
// numquery: number of sequences in the query data
// numbase: total length of alignment of reference sequences
// NUMCAT: number of categories in the gamma distribution
// TO DO: Consider which memory allocated data structures should be single continuous chunk (if possible) for best speed
void read_data(char *assignfile, char *fraclikefile, char *querydatafile, char *referencedatafile, char *errorfile, int mode, int reassign_mode)
{
	double MINLIKE = -INF;
	int i, fd;
	char tempFileName[500], strTree[500];

	usedTrees = (int *)calloc(numTrees, sizeof(int));
	if (usedTrees == NULL)
	{
		printf("ERROR: Malloc failure for usedTrees\n");
		exit(0);
	}

	// This may be a good one to make into single malloc at some point since its a known size
	pi = (double **)malloc(numTrees * sizeof(double *));
	par = (double **)malloc(numTrees * sizeof(double *));

	DATA = (int ***)malloc(numTrees * (sizeof(int **)));
	bls = (double **)malloc(numTrees * sizeof(double *));
	nodeages = (double **)malloc(numTrees * sizeof(double *));
	FRACLIKE = (double **)malloc(numTrees * (sizeof(double *)));
	maxAges = (double *)malloc(numTrees * sizeof(double));
	statevector = (double **)malloc(numTrees * sizeof(double *));
	nodeOrder = (int **)malloc(numTrees * sizeof(int *));

	totMaxAge = -INF;

	// READING IN ASSIGNMENTS
	if (VERBOSE)
		printf("Reading in species assignments\n");
	if (NULL == (infile = xfopen(assignfile)))
	{
		puts("Cannot open infile with assignments: ");
		exit(-1);
	}
	fd = fileno(infile);
	posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);

	(void)fscanf(infile, "%lu\n", &numquery);

	assignments = (int *)malloc(numquery * (sizeof(int)));
	treeAssign = (int *)malloc(numquery * (sizeof(int)));
	// assignAges = (double *)malloc(numquery * (sizeof(double)));

	// Determine how nodes are stored and determine if script is needed to convert node names etc
	{
		ChunkBuf acb;
		acb.pos = 0; acb.len = 0;
		cbuf_refill(&acb, infile);
		for (i = 0; i < numquery; i++)
		{
			// read in tree assignment first
			int v = cbuf_next_int(&acb, infile);
			treeAssign[i] = v;
			usedTrees[v]++; // Not only to know if tree is used, but also how many reads
			if (treeAssign[i] < 0 || treeAssign[i] >= numTrees)
			{
				printf("Error reading assignments for trees with %d\n", treeAssign[i]);
				exit(-1);
			}
			// node assignment
			v = cbuf_next_int(&acb, infile);
			assignments[i] = v - 1; // notice that we here convert from counting from 1 to counting from 0
		}
	}

	posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);
	xfclose(infile);

	// stores reads grouped by tree
	readsTreeSorted = (int ***)malloc(numTrees * sizeof(int **));

	// stores length of reference in each tree
	numbases = (unsigned long int *)malloc(numTrees * sizeof(unsigned long int));

	// stores the length of each read grouped by tree
	readLengthTemp = (unsigned long int **)malloc(numTrees * sizeof(unsigned long int *));

	// stores index of read in tree order data structure
	readOrder = (int *)malloc(numquery * sizeof(int));

	// stores temp error profile
	readLikeTemp = (double ****)malloc(numTrees * sizeof(double ***));

	// stores temp start pos
	startposTemp = (unsigned long int **)malloc(numTrees * sizeof(unsigned long int *));

	// stores assignments in tree sort
	tempAssignments = (int **)malloc(numTrees * sizeof(int *));

	//// --- Debug Print 1 ---
	//printf("Outer tempAssignments allocated at: %p\n", (void *)tempAssignments);
	//fflush(stdout);
	//// --- End Debug Print 1 ---

	// stores number of assignments per edge in a tree
	numReadsPerAssign = (int **)malloc(numTrees * sizeof(int *));

#pragma omp parallel for schedule(dynamic) reduction(max: totMaxAge)
	for (unsigned long int treeNum = 0; treeNum < numTrees; treeNum++)
	{
		int i, j, k, v, numcat, fd;
		double a, b, checksum;
		char tempFileName[500], strTree[500];
		unsigned long int numbase_local;

		if (usedTrees[treeNum] == 0)
		{
			// no reads were assigned to this tree, don't need to read in or allocate more memory for this tree
			continue;
		}
		// READING IN OUTPUT FROM GET_FRAC_LIKE WITH RACTIONAL LIKELIHOODS, NODE AGES, AND MORE

		readsTreeSorted[treeNum] = (int **)malloc(usedTrees[treeNum] * sizeof(int *));
		readLengthTemp[treeNum] = (unsigned long int *)malloc(usedTrees[treeNum] * sizeof(unsigned long int));
		readLikeTemp[treeNum] = (double ***)malloc(usedTrees[treeNum] * sizeof(double **));
		startposTemp[treeNum] = (unsigned long int *)malloc(usedTrees[treeNum] * sizeof(unsigned long int));
		tempAssignments[treeNum] = (int *)malloc(usedTrees[treeNum] * sizeof(int));

		sprintf(strTree, "%01lu", treeNum);
		tempFileName[0] = '\0';
		strcat(tempFileName, fraclikefile);
		strcat(tempFileName, "/");
		strcat(tempFileName, strTree);
		strcat(tempFileName, "_likelihood.txt");
		if (NULL == (infile = xfopen(tempFileName)))
		{
			printf("Cannot open infile with fractional likelihoods: %s\n", tempFileName);
			exit(-1);
		}
		fd = fileno(infile);
		posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);

		// line 1
		//  THERE IS AN ASSUMPTION THAT NUMBASE AND NUMCAT ALWAYS THE SAME, SHOULD MAKE A CHECK FOR THAT
		(void)fscanf(infile, "%lu %lu %i\n", &numseq[treeNum], &numbase_local, &numcat);
		numbases[treeNum] = numbase_local;
		if (VERBOSE)
			printf("Reading in fractional likelihoods for %lu sequences\n", numseq[treeNum]);
		if (NUMCAT != numcat)
		{
			printf("Wrong number of categories for the discretization of the gamma distribution");
			exit(-1);
		}
		statevector[treeNum] = malloc(NUMCAT * (sizeof(double)));
		// line 2
		for (i = 0; i < NUMCAT; i++)
			(void)fscanf(infile, "%lf ", &statevector[treeNum][i]);

		nodeages[treeNum] = (double *)malloc((2 * numseq[treeNum] - 1) * (sizeof(double))); // nodeages: ages of internal nodes in reference data
		bls[treeNum] = (double *)malloc((2 * numseq[treeNum] - 1) * (sizeof(double)));		// bls: branch lengths associated with each node
		nodeOrder[treeNum] = (int *)malloc((2 * numseq[treeNum] - 1) * (sizeof(int)));		// nodeOrder: order of nodes in likelihood file reflect order based on nodeage + bls

		numReadsPerAssign[treeNum] = (int *)calloc(2 * numseq[treeNum] - 1, sizeof(int));

		// line 3 - numseq + 3, node ages
		for (i = 0; i < 2 * numseq[treeNum] - 1; i++)
		{
			(void)fscanf(infile, "%i", &j);
			if (j < 0 || j > 2 * numseq[treeNum] - 2)
			{
				printf("Error reading node ages");
				exit(-1);
			}
			if (fgetc(infile) != ':')
			{
				printf("Error reading node ages");
				exit(-1);
			}
			else
			{
				(void)fscanf(infile, "%lf", &a);
				(void)fscanf(infile, "%lf", &b);
				nodeages[treeNum][j] = a;
				bls[treeNum][j] = b;
				nodeOrder[treeNum][i] = j;
				// printf("%d\t%.16f\t%.16f\n", j, a, b);
			}
		}

		(void)fscanf(infile, "%lf", &maxAges[treeNum]);
		// printf("%.16f\n", maxAges[treeNum]);
		// if(treeNum > 0 && maxAges[treeNum] > totMaxAge)
		if (maxAges[treeNum] > totMaxAge)
		{
			totMaxAge = maxAges[treeNum];
		}

		get_fractionalike(treeNum); // reads in all the fractional likelihoods
		posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);
		xfclose(infile);

		// READING IN REFRENCE SEQUENCE DATA
		tempFileName[0] = '\0';
		strcat(tempFileName, referencedatafile);
		strcat(tempFileName, "/");
		strcat(tempFileName, strTree);
		strcat(tempFileName, "_reference.txt");
		if (VERBOSE)
			printf("Reading in reference sequences: ");
		if (NULL == (infile = xfopen(tempFileName)))
		{
			printf("Cannot open infile with reference sequence data");
			exit(-1);
		}
		fd = fileno(infile);
		posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);

		// readseq returns num nodes and alters input with number of bases
		i = readseq(&j, treeNum);
		if (i != numseq[treeNum] || j != numbase_local)
		{
			printf("Number of sequences and sequence lengths do no match in input files\n");
			exit(-1);
		}

		posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);
		xfclose(infile);

		pi[treeNum] = (double *)malloc(4 * sizeof(double));
		par[treeNum] = (double *)malloc(6 * sizeof(double));

		// READING IN GTR parameters DATA
		tempFileName[0] = '\0';
		strcat(tempFileName, referencedatafile);
		strcat(tempFileName, "/");
		strcat(tempFileName, strTree);
		strcat(tempFileName, "_parameter.txt");
		if (VERBOSE)
			printf("Reading in GTR+Gamma parameters\n");
		if (NULL == (infile = xfopen(tempFileName)))
		{
			printf("Cannot open infile with GTR+Gamma parameters");
			exit(-1);
		}
		fd = fileno(infile);
		posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);

		checksum = 0.0;
		// line 1
		(void)fscanf(infile, "%lf", &a);

		// line 2: nucleotide frequencies
		for (i = 0; i < 4; i++)
		{
			(void)fscanf(infile, "%lf", &pi[treeNum][i]);
			checksum += pi[treeNum][i];
		}

		if (checksum > 1.0 + 3e-7 || checksum < 1.0 - 3e-7)
		{ // TO DO: bring back to -8 and figure out why tolerances are not being met with COI database
			printf("WARNING: Nucleotide frequencies not properly scaled or not read correctly (checksum: %.15f)\n", checksum);
			exit(-1);
		}

		for (i = 0; i < 6; i++)
		{
			(void)fscanf(infile, "%lf", &par[treeNum][i]);
		}

		posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);
		xfclose(infile);

		// printf("Printing node ages for tree %d\n", treeNum);

		// for(int i = 0; i < numseq[treeNum]; i++)
		// {
		// 	printf("\tNode %d: %.16f\n", i, nodeages[treeNum][i]);
		// }
	}

	doNRinits(2);
	inittransitionmatrix();
	// unrolled
	//  is there a reason I don't add the following to inittransitionmatrix? or into the next forloop?
	for (unsigned long int i = 0; i < numTrees; i++)
	{
		if (usedTrees[i] == 0)
		{
			// no reads were assigned to this tree, skip
			continue;
		}

		pi[i][0] = log(pi[i][0]);
		pi[i][1] = log(pi[i][1]);
		pi[i][2] = log(pi[i][2]);
		pi[i][3] = log(pi[i][3]);
	}

	tempnumquery = 0;
	startpos = NULL;
	readlength = NULL;
	assignAges = NULL;
	readlike = NULL;
	treeRoots = (unsigned long int *)malloc(sizeof(unsigned long int) * numTrees);
	trees = (struct node **)malloc(sizeof(struct node *) * numTrees);
	unsigned long int *refBasesPerTree = (unsigned long int *)calloc(numTrees, sizeof(unsigned long int));

	// READING IN QUERY DATA + ERROR PROFILE in parallel with TREE STRUCTURE READS.
	// infile is _Thread_local so each section has its own file handle.
#pragma omp parallel sections private(fd, i)
	{
#pragma omp section
		{
			// --- Query data ---
			if (VERBOSE)
				printf("Reading in query data: ");
			if (NULL == (infile = xfopen(querydatafile)))
			{
				puts("Cannot open infile with query data\n");
				exit(-1);
			}
			fd = fileno(infile);
			posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);

			// Line 1, number of queries
			(void)fscanf(infile, "%lu", &numquery);
			printf("There are %lu query sequences\n", numquery);

			// This expression also reads in data to QUERYDATA
			if (read_query_data(numquery) != numquery)
			{
				printf(" Different number of query sequences found in query data file and in assignment file\n");
				exit(-1);
			}

			posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);
			xfclose(infile);

			// Testing if read assignments are valid
			// Here because we need to know number of references in the tree
			if (VERBOSE)
				printf("Testing validity of node assignments.\n");
			for (i = 0; i < numquery; i++)
			{
				if (assignments[i] < 0 || assignments[i] > 2 * numseq[treeAssign[i]] - 1)
				{
					printf("Error reading assignments for read %d with assignment %d from tree %d with %lu references\n", i, assignments[i], treeAssign[i], numseq[treeAssign[i]]);
					exit(-1);
				}
			}

			// --- Error profile ---
			if (VERBOSE)
				printf("Reading in error profile\n");
			if (NULL == (infile = xfopen(errorfile)))
			{
				puts("Cannot open infile with error profile: ");
				exit(-1);
			}
			fd = fileno(infile);
			posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);

			make_readfraclike();

			posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);
			xfclose(infile);
		}

#pragma omp section
		{
			// --- Tree structure reads (independent of query/error data) ---
			for (unsigned long int treeNum = 0; treeNum < numTrees; treeNum++)
			{
				if (usedTrees[treeNum] == 0)
					continue;

				char tfn[500], st[500];
				sprintf(st, "%01lu", treeNum);
				tfn[0] = '\0';
				strcat(tfn, referencedatafile);
				strcat(tfn, "/");
				strcat(tfn, st);
				strcat(tfn, "_reference.txt");

				if (NULL == (infile = xfopen_seekable(tfn)))
				{
					printf("Cannot open infile with tree data\n");
					exit(-1);
				}
				int fd_t = fileno(infile);

				unsigned long int numRef_t, refBases_t;
				(void)fscanf(infile, "%lu %lu", &numRef_t, &refBases_t);
				refBasesPerTree[treeNum] = refBases_t;

				fseek(infile, 0, SEEK_END);
				// Move the file pointer to the beginning of the last line
				fseek(infile, 0, SEEK_SET);

				char ch; int line_length = 0, last_line_length = 0;
				while ((ch = fgetc(infile)) != EOF)
				{
					line_length++;
					if (ch == '\n')
					{
						last_line_length = line_length;
						line_length = 0;
					}
				}
				fseek(infile, -last_line_length, SEEK_END);

				// Globals from Rasmus code, may want to fix...
				tip = 0; comma = 0;

				allocatetreememmory(numseq[treeNum], treeNum);
				treeRoots[treeNum] = getclade(numseq[treeNum], treeNum) - 1 + numseq[treeNum]; // converts to ratePlacer node

				posix_fadvise(fd_t, 0, 0, POSIX_FADV_DONTNEED);
				xfclose_seekable(infile);
			}
		}
	} /* end omp parallel sections */

	// Serial: reassign + merge (depends on both query/error data and tree structures)
	for (unsigned long int treeNum = 0; treeNum < numTrees; treeNum++)
	{
		if (usedTrees[treeNum] == 0)
			continue;

		//printtree(numseq[treeNum], treeRoots[treeNum], treeNum);
		if (reassign_mode == 1) {
			bestAssignment(treeRoots[treeNum], treeNum);
		} else if (reassign_mode == 2) {
			tronkoAssignmentTesting(treeRoots[treeNum], treeNum);
		}

		// Depending on mode, we don't need to merge! Should make a function/edit mergeReads that just transfers to the correct datastructures for the rest of ratePlacer
		if (toMerge)
		{
			mergeReads(treeNum, refBasesPerTree[treeNum]);
		}
		else
		{
			keepSeparateReads(treeNum);
		}

		for (unsigned long int i = 0; i < usedTrees[treeNum]; i++)
		{
			free(readsTreeSorted[treeNum][i]);
		}
		free(readsTreeSorted[treeNum]);

		// No need to keep in memory
		// freetreememmory();
	}

	free(refBasesPerTree);

	numquery = tempnumquery;

	printf("After merging, there are now %lu sequences\n", numquery);

	if (numquery < 1)
	{
		for (unsigned long int i = 0; i < numTrees; i++)
		{
			if (usedTrees[i] == 0)
			{
				// no reads were assigned to this tree, don't need to free memory here
				continue;
			}
			free(numReadsPerAssign[i]);
			for (unsigned long int j = 0; j < usedTrees[i]; j++)
			{
				for (unsigned long int k = 0; k < readLengthTemp[i][j]; k++)
				{
					free(readLikeTemp[i][j][k]);
				}
				free(readLikeTemp[i][j]);
			}
			free(readLikeTemp[i]);
			free(readLengthTemp[i]);
			printf("Freeing tempAssignments[%lu] %p\n", i, tempAssignments[i]);
			free(tempAssignments[i]);
			printf("Freed tempAssignments[%lu] %p\n", i, tempAssignments[i]);
		}
		free(readLikeTemp);

		free(readsTreeSorted);
		free(readOrder);
		free(readLengthTemp);
		free(startposTemp);
		free(tempAssignments);
		free(numReadsPerAssign);
		printf("No reads pass required reference coverage of 5%%. Exiting.\n");
		exit(0);
	}

	// deallocate extra memory if needed/tree
	int* tempptr = (int *)realloc(treeAssign, tempnumquery * (sizeof(int)));
	if (tempptr != NULL)
	{
		treeAssign = tempptr;
	} else {
		printf("Failed to reallocate treeAssign in read_data\n");
		exit(0);
	}

	tempptr = (int *)realloc(assignments, tempnumquery * (sizeof(int)));
	if (tempptr != NULL)
	{
		assignments = tempptr;
	} else {
		printf("Failed to reallocate assignments in read_data\n");
		exit(0);
	}

	unsigned long int *tempptr_ul = (unsigned long int *)realloc(startpos, tempnumquery * (sizeof(unsigned long int)));
	if (tempptr_ul != NULL)
	{
		startpos = tempptr_ul;
	} else {
		printf("Failed to reallocate startpos in read_data\n");
		exit(0);
	}

	tempptr_ul = (unsigned long int *)realloc(readlength, tempnumquery * (sizeof(unsigned long int)));
	if (tempptr_ul != NULL)
	{
		readlength = tempptr_ul;
	} else {
		printf("Failed to reallocate readlength in read_data\n");
		exit(0);
	}

	int **tempptr2 = (int **)realloc(QUERYDATA, tempnumquery * (sizeof(int *)));
	if (tempptr2 != NULL)
	{
		QUERYDATA = tempptr2;
	} else {
		printf("Failed to reallocate QUERYDATA in read_data\n");
		exit(0);
	}

	double *tempptr3 = (double *)realloc(assignAges, tempnumquery * (sizeof(double)));
	if (tempptr3 != NULL)
	{
		assignAges = tempptr3;
	} else {
		printf("Failed to reallocate assignAges in read_data\n");
		exit(0);
	}

	double ***tempptr4 = (double ***)realloc(readlike, tempnumquery * sizeof(double **));
	if (tempptr4 != NULL)
	{
		readlike = tempptr4;
	} else {
		printf("Failed to reallocate readlike in read_data\n");
		exit(0);
	}
	//QUERYDATA = (int **)realloc(QUERYDATA, tempnumquery * (sizeof(int *)));
	//startpos = (int *)realloc(startpos, tempnumquery * (sizeof(int)));
	//readlength = (int *)realloc(readlength, tempnumquery * (sizeof(int)));
	//assignAges = (double *)realloc(assignAges, tempnumquery * (sizeof(double)));
	//assignments = (int *)realloc(assignments, tempnumquery * (sizeof(int)));
	//treeAssign = (int *)realloc(treeAssign, tempnumquery * (sizeof(int)));
	//readlike = (double ***)realloc(readlike, tempnumquery * sizeof(double **));

	//for(int read = 0; read < numquery; read++)
	//{
	//	printf("Query %d is length %d starting as pos %d assigned to tree %d and edge %d\n", read, readlength[read], startpos[read], treeAssign[read], assignments[read]);
	//	for(int pos = 0; pos < readlength[read]; pos++)
	//	{
	//		if(QUERYDATA[read][pos] == -1)
	//			printf("-");
	//		else if (QUERYDATA[read][pos] == 0)
	//			printf("A");
	//		else if (QUERYDATA[read][pos] == 1)
	//			printf("C");
	//		else if (QUERYDATA[read][pos] == 2)
	//			printf("G");
	//		else if (QUERYDATA[read][pos] == 3)
	//			printf("T");
	//		else
	//		{
	//			printf("%d is not valid", QUERYDATA[read][pos]);
	//			exit(0);
	//		}

	//	}
	//	printf("\n");
	//	//for(int pos = 0; pos < readlength[read]; pos++)
	//	//{
	//	//	if (QUERYDATA[read][pos] != -1)
	//	//	{
	//	//		printf("pos %d:", pos);
	//	//		for(int k = 0; k < 4; k++)
	//	//		{
	//	//			printf("%lf,", readlike[read][pos][k]);
	//	//		}
	//	//		printf("\n");
	//	//	}
	//	//}
	// }

	//// --- Debug Print 2 ---
	//printf("Outer tempAssignments before cleanup loop: %p\n", (void *)tempAssignments);
	//fflush(stdout);
	//// --- End Debug Print 2 ---

	for (unsigned long int i = 0; i < numTrees; i++)
	{
		if (usedTrees[i] == 0)
		{
			// no reads were assigned to this tree, don't need to free memory here
			continue;
		}
		free(numReadsPerAssign[i]);
		for (unsigned long int j = 0; j < usedTrees[i]; j++)
		{
			for (unsigned long int k = 0; k < readLengthTemp[i][j]; k++)
			{
				free(readLikeTemp[i][j][k]);
			}
			free(readLikeTemp[i][j]);
		}
		free(startposTemp[i]);
		free(readLikeTemp[i]);
		free(readLengthTemp[i]);
		//printf("Freeing tempAssignments[%d] %p\n", i, tempAssignments[i]);
		free(tempAssignments[i]);
		//printf("Freed tempAssignments[%d] %p\n", i, tempAssignments[i]);
	}
	free(readLikeTemp);

	free(readsTreeSorted);
	free(readOrder);
	free(readLengthTemp);
	free(startposTemp);
	free(tempAssignments);
	free(numReadsPerAssign);

	//exit(0);
}

void maximize_like_seperately_for_all2D_Print(double **par)
{
	printf("Starting maximize seperately for all\n");
	int i, k, v, nfun;
	double p[3];
	double L1, invector[3], lowbound[3], upbound[3], eh0 = 3e-8;

	onDindic = 0;

	// assignmentMode of 0 means single assignment given
	for (i = 0; i < numquery; i++)
	{
		p[0] = i;
		p[1] = treeAssign[i];
		p[2] = assignments[i];

		invector[1] = 0.5;
		invector[2] = bls[treeAssign[i]][assignments[i]] / 2.0;
		lowbound[1] = eh0;
		lowbound[2] = eh0;
		upbound[1] = 1.0 - eh0;
		upbound[2] = bls[treeAssign[i]][assignments[i]] - eh0;
		nfun = 0;

		printf("Sequence %d\n", i);
		L1 = findmax_amoeba(invector, lowbound, upbound, 2, getlike_gamma_root_in_trifurcation_Print_Lik, p, 3);
		printf("\tParameter estimates %.16f %.16f: %.16f\n", invector[1], invector[2], L1);
		printf("\tAssignment %d of tree %d age: %.16f\n", assignments[i], treeAssign[i], nodeages[treeAssign[i]][assignments[i]]);
		printf("\tsequence age: %.16f\n", (1.0 - invector[1]) * (nodeages[treeAssign[i]][assignments[i]] + invector[2]));
		printf("Site scores:\n");
		getlike_gamma_root_in_trifurcation_Print(invector, p);
	}
}

void maximize_like_seperately_for_all2D(double **par)
{
	printf("Starting maximize seperately for all\n");
	int i, k, v, nfun;
	double p[3];
	double L1, invector[3], lowbound[3], upbound[3], eh0 = 3e-8;

	onDindic = 0;

	// assignmentMode of 0 means single assignment given
	for (i = 0; i < numquery; i++)
	{
		p[0] = i;
		p[1] = treeAssign[i];
		p[2] = assignments[i];

		printf("Sequence %d\n", i);

		invector[1] = 0.5;
		invector[2] = bls[treeAssign[i]][assignments[i]] / 2.0;
		lowbound[1] = eh0;
		lowbound[2] = eh0;
		upbound[1] = 1.0 - eh0;
		upbound[2] = bls[treeAssign[i]][assignments[i]] - eh0;
		nfun = 0;

		// printf("\tParameter initial %.16f %.16f\n", invector[1],invector[2]);

		L1 = findmax_amoeba(invector, lowbound, upbound, 2, getlike_gamma_root_in_trifurcation, p, 3);

		printf("\tParameter estimates %.16f %.16f: %.16f\n", invector[1], invector[2], L1);
		printf("\tAssignment %d age: %.16f\n", assignments[i], nodeages[treeAssign[i]][assignments[i]]);
		printf("\tsequence age: %.16f\n", (1.0 - invector[1]) * (nodeages[treeAssign[i]][assignments[i]] + invector[2]));

		invector[0] = 0.5;
		invector[1] = bls[treeAssign[i]][assignments[i]] / 2.0;

		L1 = minimize_brent(invector, 2, getlike_gamma_root_in_trifucation_single_read_brent, 100, &i);

		if (invector[0] < 3e-8)
		{
			invector[0] = 3e-8;
		}
		else if (invector[0] > 1.0 - 3e-8)
		{
			invector[0] = 1.0 - 3e-8;
		}

		if (invector[1] < 3e-8)
		{
			invector[1] = 3e-8;
		}
		else if (invector[1] > bls[treeAssign[i]][assignments[i]] - 3e-8)
		{
			invector[1] = bls[treeAssign[i]][assignments[i]] - 3e-8;
		}

		printf("\tParameter estimates brent %.16f %.16f: %.16f\n", invector[0], invector[1], L1);
        printf("\tAssignment %d age: %.16f\n", assignments[i], nodeages[treeAssign[i]][assignments[i]]);
        printf("\tsequence age: %.16f\n",(1.0-invector[0])*(nodeages[treeAssign[i]][assignments[i]]+invector[1]));


	}
}

void likelihoodratiotest_for_all(double **par, double compareAge)
{
	printf("Starting maximize seperately with likelihood ratio test for all\n");
	int i, k, v, nfun;
	double p[3], L1, L2;
	double invector[3], lowbound[3], upbound[3], eh0 = 3e-8;

	// assignmentMode of 0 means single assignment given
	for (i = 0; i < numquery; i++)
	{
		p[0] = i;
		p[1] = treeAssign[i];
		p[2] = assignments[i];

		invector[1] = 0.5;
		invector[2] = bls[treeAssign[i]][assignments[i]] / 2.0;
		lowbound[1] = eh0;
		lowbound[2] = eh0;
		upbound[1] = 1.0 - eh0;
		upbound[2] = bls[treeAssign[i]][assignments[i]] - eh0;
		// printf("Checking bls %.16f and node age %.16f\n", upbound[2], upbound[1]);
		nfun = 0;
		onDindic = 0;
		printf("Sequence %d\n", i);
		L1 = findmax_amoeba(invector, lowbound, upbound, 2, getlike_gamma_root_in_trifurcation, p, 3);

		printf("\tParameter estimates %.16f %.16f (%.16f): %.16f\n", invector[1], invector[2], bls[treeAssign[i]][assignments[i]], L1);
		printf("\tAssignment %d age of tree %d: %.16f\n", assignments[i], treeAssign[i], nodeages[treeAssign[i]][assignments[i]]);
		printf("\tsequence age: %.16f\n", (1.0 - invector[1]) * (nodeages[treeAssign[i]][assignments[i]] + invector[2]));

		invector[1] = bls[treeAssign[i]][assignments[i]] / 2.0;
		lowbound[1] = eh0;
		upbound[1] = bls[treeAssign[i]][assignments[i]] - eh0;
		invector[2] = -1;
		nfun = 0;
		onDindic = 1;

		// NEEDS TO BE SWITCHED TO GOLDEN SECTION
		// L2 = findmax_amoeba(invector,lowbound, upbound, 1, getlike_gamma_root_in_trifurcation_L0, p, 3);
		if (testAge < eh0)
		{
			testAge = eh0;
		}

		testAge = compareAge;
		L2 = Brent1D(invector, lowbound, upbound, 1, getlike_gamma_root_in_trifurcation_testAge, p, 3);

		printf("\tParameter estimates llr %.16f: %.16f\n", invector[1], L2);
		printf("\tLikelihood ratio statistic to %.16f: %.16f\n", testAge, 2.0 * (L1 - L2));
	}
}

void readContour(double **par)
{
	printf("Likelihood contour points being generated\n");

	int i, k, v;
	double L, parameters[7], times[3], eh0 = 3e-8, bl_increment, alpha_increment = 0.001;

	onDindic = 0;

	// assignmentMode of 0 means single assignment given
	// parameters[0] = sequence
	// parameters[1] = treeNumber
	// parameters[2] = assignment
	// parameters[3] = lower bound of alpha
	// parameters[4] = upper bound of alpha
	// parameters[5] = lower bound of bl
	// parameters[6] = upper bound of bl
	// times[0] = ?
	// times[1] = alpha (in_alpha)
	// times[2] = root placement (in_bl)

	parameters[3] = eh0;
	parameters[4] = 1.0 - eh0;
	parameters[5] = eh0;

	for (i = 0; i < numquery; i++)
	{
		parameters[0] = i;
		parameters[1] = treeAssign[i];
		parameters[2] = assignments[i];
		parameters[6] = bls[treeAssign[i]][assignments[i]];

		times[1] = eh0;

		bl_increment = bls[treeAssign[i]][assignments[i]] / 1000;

		for (k = 0; k < 1000; k++)
		{
			times[2] = eh0;

			for (v = 0; v < 1000; v++)
			{
				L = getlike_gamma_root_in_trifurcation(times, parameters);

				printf("%d,%.16f,%.16f,%.16f,%.16f\n", i, times[1], times[2], L, (1.0 - times[1]) * (nodeages[treeAssign[i]][assignments[i]] + times[2]));

				times[2] += bl_increment;
			}
			times[1] += alpha_increment;
		}

		// Max not always calculated otherwise
		times[1] = eh0;
		times[2] = bls[treeAssign[i]][assignments[i]] - eh0;

		for (k = 0; k < 1000; k++)
		{
			L = getlike_gamma_root_in_trifurcation(times, parameters);

			//printf("%d,%.16f,%.16f,%.16f\n", i, times[1], times[2], L);
			printf("%d,%.16f,%.16f,%.16f,%.16f\n", i, times[1], times[2], L, (1.0 - times[1]) * (nodeages[treeAssign[i]][assignments[i]] + times[2]));

			times[1] += alpha_increment;
		}
	}
}

// Will "drop reads" by changing readStart based on timeInc and test age
// New update: Will not be based on the node ages themselves as we're now making
// this for across trees
double dropReads(int *readStart, double timeInc)
{
	int readDropped = 0;
	double test = testAge;

	do
	{
		test += timeInc;
		// TO DO: Try to make assignAges work here since that would be faster!!
		while (*readStart < numquery - 1 && nodeages[treeAssign[usedReads[*readStart]]][assignments[usedReads[*readStart]]] + bls[treeAssign[usedReads[*readStart]]][assignments[usedReads[*readStart]]] <= test) // similar -2 as above
		{
			*readStart = *readStart + 1;
			readDropped = 1;
			// printf("Read dropped!\n");
		}
	} while (*readStart < numquery - 1 && readDropped == 0);

	return nodeages[treeAssign[usedReads[*readStart]]][assignments[usedReads[*readStart]]] + bls[treeAssign[usedReads[*readStart]]][assignments[usedReads[*readStart]]];
}

// merge functionality, using in-place merge sort to reduce memory load
//  Will base merge on nodeAge[i]+nodebl[i] to get age in which reads assigned to i will need to cut in ratePlacer
//  for reference: https://www.geeksforgeeks.org/merge-sort/ and https://www.geeksforgeeks.org/in-place-merge-sort/
void merge(int left, int right, int mid)
{
	// If already merged, are there cases where this may be true but not sorted?
	if (assignAges[usedReads[mid]] <= assignAges[usedReads[mid + 1]])
	{
		return;
	}

	unsigned long int midstart = mid + 1;

	while (left <= mid && midstart <= right)
	{
		// if elements are in right place (ie left < right)
		if (assignAges[usedReads[left]] <= assignAges[usedReads[midstart]])
		{
			left++;
		}
		else
		{
			unsigned long int index = midstart, tempOrder = usedReads[midstart];
			// double tempAge = assignAges[midstart];

			// shift all elements between left and midstart to right by 1
			while (index != left)
			{
				usedReads[index] = usedReads[index - 1];
				// assignAges[index] = assignAges[index - 1];
				index--;
			}
			usedReads[left] = tempOrder;
			// assignAges[left] = tempAge;

			// Update
			left++;
			mid++;
			midstart++; // may be redundant, but removes some extra calc?
		}
	}
}

// merge sort implementation for sort_ages
void mergeSort(int left, int right)
{
	if (left < right)
	{
		// Mid point calculation to hopefully avoid overflow
		int mid = left + (right - left) / 2;

		// recurse down
		mergeSort(left, mid);
		mergeSort(mid + 1, right);

		// merge
		merge(left, right, mid);
	}
}

// Initialization where reads are ordered by their assignment node order
void orderReads()
{
	// usedReads = (int *)malloc(sizeof(int) * numquery);
	//  TO DO: Determine if this is the most efficient way
	//  or if I should have merge sort work on all data associated with
	//  the reads instead
	for (unsigned long int i = 0; i < numquery; i++)
	{
		usedReads[i] = i;
	}

	mergeSort(0, numquery - 1);

	// printf("Used reads: ");
	// for(int i = 0; i < numquery; i++)
	//{
	//	printf("%d(%d) ", usedReads[i], assignments[usedReads[i]]);
	//	//printf("%d is age %.16f and %.16f\n", i, assignAges[usedReads[i]], nodeages[treeAssign[usedReads[i]]][assignments[usedReads[i]]] + bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]]);
	// }
	// printf("\n");
}

// getlike_ages but for a single tree
double getlike_ages_tree(double times, double parameters[7])
{
	// age for optimization should be times, is redudant and should fix?
	testAge = times;

	if (testAge < 0.0)
	{
		return 1000000000.0;
	}

	// The following should be thought about because it would be inconvienent to redo for each age
	// May want to make userReads a global, but think after this is implemented and working
	// TO DO TO DO TO DO!!!!
	int readStart = parameters[0], tree = parameters[1], i, nfun;
	double invector[3], lowbound[3], upbound[3], eh0 = 3e-8, age_like = 0.0, ageIncr, L2;
	double p[3];

	onDindic = 1;

	// Drop any reads that cannot be of test age
	for (i = readStart; i < numquery; i++)
	{
		if (treeAssign[usedReads[i]] != tree)
		{
			continue;
		}

		p[0] = usedReads[i];
		p[1] = treeAssign[usedReads[i]];
		p[2] = assignments[usedReads[i]];

		nfun = 0;
		invector[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]] / 2.0;
		lowbound[1] = eh0;
		upbound[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]] - eh0;

		L2 = Brent1D(invector, lowbound, upbound, 1, getlike_gamma_root_in_trifurcation_testAge, p, 3);

		age_like += L2; // sum of log likelihoods
	}

	return (age_like);
}

double getlike_ages_rough(double times, double parameters[7])
{
	// age for optimization should be times, is redudant and should fix?
	// printf("Testing age %.16f in getlike\n", times);
	testAge = times;

	// printf("get like test age of %.16f\n", testAge);

	if (testAge < 0.0)
	{
		return 1000000000.0;
	}

	// The following should be thought about because it would be inconvienent to redo for each age
	// May want to make userReads a global, but think after this is implemented and working
	// TO DO TO DO TO DO!!!!
	int readStart = parameters[0], i, nfun;
	double invector[3], lowbound[3], upbound[3], eh0 = 3e-8, age_like = 0.0, ageIncr, L2;
	double p[3];

	onDindic = 1;

	// Drop any reads that cannot be of test age
	for (i = readStart; i < numquery; i++)
	{
		// if(usedReads[i] != 84)
		//	continue;
		p[0] = usedReads[i];
		p[1] = treeAssign[usedReads[i]];
		p[2] = assignments[usedReads[i]];

		nfun = 0;
		invector[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]] / 2.0;
		lowbound[1] = eh0;
		upbound[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]] - eh0;

		L2 = GoldenSection_rough(invector, lowbound, upbound, 1, getlike_gamma_root_in_trifurcation_testAge, p, 3);

		// printf("Read %d\tTime: %.16f\tLikelihood: %.16f\n", usedReads[i], testAge, L2);

		// printf("\tRead %d assigned to tree %d and node %d with likelihood %.16f with root placement of %.16f\n", usedReads[i], treeAssign[usedReads[i]], assignments[usedReads[i]], L2, nodeages[treeAssign[usedReads[i]]][assignments[usedReads[i]]]+invector[1]);

		// printf("\tRead %d assigned to tree %d and node %d with likelihood %.16f with root placement of %.16f\n", usedReads[i], treeAssign[usedReads[i]], assignments[usedReads[i]], L2, nodeages[treeAssign[usedReads[i]]][assignments[usedReads[i]]]+invector[1]);
		age_like += L2; // sum of log likelihoods
						// printf("\t\tassignment %d like contribution %.16f\n", assignments[usedReads[i]], L2);
	}
	// printf("Likelihood of age %.16f: %.16f\n", testAge, age_like);

	// printf("\tlikelihood sum of %.16f\n", age_like);

	return (age_like);
}

double getlike_ages(double times, double parameters[7])
{
	// age for optimization should be times, is redudant and should fix?
	// printf("Testing age %.16f in getlike\n", times);
	testAge = times;

	// printf("get like test age of %.16f\n", testAge);

	if (testAge < 0.0)
	{
		return 1000000000.0;
	}

	// The following should be thought about because it would be inconvienent to redo for each age
	// May want to make userReads a global, but think after this is implemented and working
	// TO DO TO DO TO DO!!!!
	int readStart = parameters[0], i;
	double eh0 = 3e-8, age_like = 0.0, ageIncr;
	//for(int j = 0; j < numquery; j++)
	//{
	//	brent_invector[j] = nodeages[treeAssign[j]][assignments[j]] + bls[treeAssign[j]][assignments[j]]/2;
	//}

	//age_like = minimize_brent(brent_invector, numquery, getlike_gamma_root_in_trifucation_sample_age_brent, 10000, &testAge);

	// Compute per-read likelihoods in parallel, store in fixed-indexed array to avoid
	// floating-point non-associativity from OpenMP reduction ordering.
	double L2_values[numquery - readStart];

	#pragma omp parallel
	{
		onDindic = 1; // each worker thread must set its own thread-local copy
		#pragma omp for schedule(dynamic, 1)
		for (i = readStart; i < numquery; i++)
		{
			// if(usedReads[i] != 85)
			//	continue;
			double p[3], invector[3], lowbound[3], upbound[3];
			int nfun = 0;
			p[0] = usedReads[i];
			p[1] = treeAssign[usedReads[i]];
			p[2] = assignments[usedReads[i]];

			invector[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]] / 2.0;
			lowbound[1] = eh0;
			upbound[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]] - eh0;

			L2_values[i - readStart] = Brent1D(invector, lowbound, upbound, 1, getlike_gamma_root_in_trifurcation_testAge, p, 3);

			// printf("Read %d\tTime: %.16f\tLikelihood: %.16f\n", usedReads[i], testAge, L2_values[i - readStart]);

			// printf("\tRead %d assigned to tree %d and node %d with likelihood %.16f with root placement of %.16f\n", usedReads[i], treeAssign[usedReads[i]], assignments[usedReads[i]], L2_values[i - readStart], nodeages[treeAssign[usedReads[i]]][assignments[usedReads[i]]]+invector[1]);
		}
	}

	// Sum serially in fixed order — deterministic regardless of thread count
	for (i = 0; i < (int)(numquery - readStart); i++)
	{
		age_like += L2_values[i];
		// printf("\t\tassignment %d like contribution %.16f\n", assignments[usedReads[readStart + i]], L2_values[i]);
	}
	// printf("Likelihood of age %.16f: %.16f\n", testAge, age_like);

	// printf("\tlikelihood sum of %.16f\n", age_like);

	return (age_like);
}

double getlike_ages_brent(double times, double parameters[7])
{
	// age for optimization should be times, is redudant and should fix?
	// printf("Testing age %.16f in getlike\n", times);
	testAge = times;

	// printf("get like test age of %.16f\n", testAge);

	if (testAge < 0.0)
	{
		return 1000000000.0;
	}

	// The following should be thought about because it would be inconvienent to redo for each age
	// May want to make userReads a global, but think after this is implemented and working
	// TO DO TO DO TO DO!!!!
	//int readStart = parameters[0], i, nfun;
	//double invector[3], lowbound[3], upbound[3], eh0 = 3e-8, age_like = 0.0, ageIncr, L2;
	double age_like = 0.0;
	//double p[3];
	double brent_invector[numquery];

	for(int j = 0; j < numquery; j++)
	{
		brent_invector[j] = nodeages[treeAssign[j]][assignments[j]] + bls[treeAssign[j]][assignments[j]]/2;
	}

	age_like = minimize_brent(brent_invector, numquery, getlike_gamma_root_in_trifucation_sample_age_brent, 10000, &testAge);

	//onDindic = 1;

	//// Drop any reads that cannot be of test age
	//for (i = readStart; i < numquery; i++)
	//{
	//	// if(usedReads[i] != 85)
	//	//	continue;
	//	p[0] = usedReads[i];
	//	p[1] = treeAssign[usedReads[i]];
	//	p[2] = assignments[usedReads[i]];

	//	nfun = 0;
	//	invector[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]] / 2.0;
	//	lowbound[1] = eh0;
	//	upbound[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]] - eh0;

	//	L2 = GoldenSection(invector, lowbound, upbound, 1, getlike_gamma_root_in_trifurcation_testAge, p, 3);

	//	// printf("Read %d\tTime: %.16f\tLikelihood: %.16f\n", usedReads[i], testAge, L2);

	//	// printf("\tRead %d assigned to tree %d and node %d with likelihood %.16f with root placement of %.16f\n", usedReads[i], treeAssign[usedReads[i]], assignments[usedReads[i]], L2, nodeages[treeAssign[usedReads[i]]][assignments[usedReads[i]]]+invector[1]);

	//	age_like += L2; // sum of log likelihoods
	//					// printf("\t\tassignment %d like contribution %.16f\n", assignments[usedReads[i]], L2);
	//}
	// printf("Likelihood of age %.16f: %.16f\n", testAge, age_like);

	// printf("\tlikelihood sum of %.16f\n", age_like);

	return (age_like);
}

// Gets Fisher information bounds for confidence interval
// finite method for second derivative
//(f(x+h) - 2f(x) + f(x-h))/h^2
//  optAge is x, f(x) is optLik
//  We are giving -lik, so remember (and test) to convert back again. Test if this is really needed
double confidenceIntervalFisher(double **par, double optAge, double optLik, int readStart)
{
	// Step size for centered finite-difference second derivative.
	// Optimal h ~ eps^(1/4) * |x| ~ 4e-4 * optAge; floor prevents
	// catastrophic cancellation when optAge is very small.
	double h = fmax(4e-4 * optAge, 1e-6), f_hx = 0.0, fx_h = 0.0, nSecondDeriv;
	int i, k, v, nfun;
	double p[3];
	double invector[3], lowbound[3], upbound[3], eh0 = 3e-8;

	// Drop any reads that cannot be of test age
	for (i = readStart; i < numquery; i++)
	{
		p[0] = usedReads[i];
		p[1] = treeAssign[usedReads[i]];
		p[2] = assignments[usedReads[i]];

		nfun = 0;
		invector[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]] / 2.0;
		lowbound[1] = eh0;
		upbound[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]] - eh0;

		// Does there need to be a test if this passes some internal node that reads need to be dropped in?
		testAge = optAge + h;
		f_hx += Brent1D(invector, lowbound, upbound, 1, getlike_gamma_root_in_trifurcation_testAge, p, 3);

		invector[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]] / 2.0;
		lowbound[1] = eh0;
		upbound[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]] - eh0;

		testAge = optAge - h;
		fx_h += Brent1D(invector, lowbound, upbound, 1, getlike_gamma_root_in_trifurcation_testAge, p, 3);
	}

	// printf("f(x-h): %.16f\nf(x+h):%.16f\nf(x) - f(x-h): %.16f\nf(x) - f(x+h): %.16f\nSec Deriv: %.16f\n", fx_h, f_hx, optLik - fx_h, optLik - f_hx, (2 * optLik - fx_h - f_hx)/(pow(h,2.0)));

	// equivalent to f_hx - 2*f_x + fx_h if getlike_gamma didn't return -Lik
	// Potentially implement as f_xh - 2 * f_x + fx_h so that -1 already incorporated to take the second derivate for
	// return((2 * optLik - fx_h - f_hx)/(pow(h,2.0)));

	//-1 * second derivative
	nSecondDeriv = (f_hx - 2 * optLik + fx_h) / (pow(h, 2.0));

	// printf("%.16f\n", nSecondDeriv);

	return (1.96 / sqrt(nSecondDeriv));
}

void maximize_like_jointly_for_all2D(double **par, int allTrees)
{
	printf("Starting maximize jointly for all\n");
	int i, k, v, nfun, readStart = 0, nodePointer = 0, twice = 0, oldReadStart, oldNodePointer;
	double p[3], L1, L2, secD, confI;
	double invector[3], lowbound[3], upbound[3], eh0 = 3e-8, nextNodeAge, est_age, est_age_lik, oldNodeAge;
	double incr = totMaxAge / 1000, oldIn1;
	usedReads = (int *)malloc(sizeof(int) * numquery);
	// printf("Max age of %.16f and incr of %.16f\n", totMaxAge, incr);

	orderReads();

	// This nextNodeAge represents the max without dropping reads
	nextNodeAge = nodeages[treeAssign[usedReads[readStart]]][assignments[usedReads[readStart]]] + bls[treeAssign[usedReads[readStart]]][assignments[usedReads[readStart]]];

	onDindic = 1;

	if (numquery > 1)
	{
		printf("Starting rough estimation for maximum bound\n");

		//p[0] = readStart;
		//p[1] = nodePointer;
		//p[2] = 0;
		//invector[1] = nextNodeAge / 2;
		//lowbound[1] = eh0;
		//upbound[1] = nextNodeAge - eh0;
		//nfun = 0;
		//onDindic = 0;

		//L1 = Brent1D(invector, lowbound, upbound, 1, getlike_ages, p, 3);

		//do{

		//	oldNodePointer = nodePointer;
		//	oldReadStart = readStart;
		//	oldNodeAge = nextNodeAge;
		//	oldIn1 = invector[1];
		//	//testAge = nextNodeAge;
		//	nextNodeAge = dropReads(&readStart, incr);
		//	p[0] = readStart;
		//	p[1] = nodePointer;
		//	p[2] = 0;
		//	invector[1] = nextNodeAge / 2;
		//	lowbound[1] = eh0;
		//	upbound[1] = nextNodeAge - eh0;

		//	L1 = getlike_ages(oldIn1, p);
		//	L2 = Brent1D(invector, lowbound, upbound, 1, getlike_ages, p, 3);

		//	printf("L1(%.16f): %.16f\tL2(%.16f): %.16f\tLLR: %.16f\t%d of %d\n", oldIn1, L1, invector[1], L2, 2 * (L1 - L2), readStart, numquery);

		//	//if (2 * (L1 - L2) < 1.0){
		//	//	printf("Old age is within 1.0 of new age, stopping\n");
		//	//	break;
		//	//}
		//}while(numquery - readStart> 1);

		// double time2 = (double) clock()/CLOCKS_PER_SEC;
		//  Rough optimization to find upper bound and reads to drop without running too much optimization
		//  To change: Instread of making this per assigned node, do it via time slices to make a little faster
		//  		Also maybe a GoldenSection search that is a little less stringent?
		//  		Am just trying to do a rough optimization - maybe if time is still long despite time slice change
		//  		Also maybe incorporate more options for the confidence interval, not just 95%
		//  		Would need to be able to calculate the z-score from the user given value
		//  		confI = 1.96 / sqrt(-secD);
		do
		{
			// printf("Testing max age of %.16f\n", nextNodeAge);
			// add fillers
			p[0] = readStart;
			p[1] = nodePointer;
			p[2] = 0;

			// Should consider the best way to pick these values...
			L1 = getlike_ages_rough(nextNodeAge - nextNodeAge / 100.0, p);
			L2 = getlike_ages_rough(nextNodeAge - nextNodeAge / 10.0, p);
			// printf("Results of L1 %.16f and L2 %.16f with max at %.16f with readStart at %d\n", L1, L2, nextNodeAge, readStart);
			//  Likelihood near bound is better, so drop and test next age range
			if (L1 <= L2)
			{
				// TO DO: scale time increase by maxAge
				oldNodePointer = nodePointer;
				oldReadStart = readStart;
				oldNodeAge = nextNodeAge;
				testAge = nextNodeAge;
				nextNodeAge = dropReads(&readStart, incr);
				// printf("%d of %d\n", readStart, numquery);
				if (readStart >= numquery - 1) // no reads left
				{
					// maybe exit if the likelihoods are still very different...
					if (L2 - L1 < 1.0)
					{
						printf("Warning: Using last round of dropped reads\n");
						break;
					}
					else
					{
						printf("Error: Could not optimize for bounds.\n");
						exit(0);
					}
				}
				twice = 0;
			}
			if (L1 > L2)
			{
				twice++;
				if (twice > 1)
					break;
				oldNodePointer = nodePointer;
				oldReadStart = readStart;
				oldNodeAge = nextNodeAge;
				testAge = nextNodeAge;
				nextNodeAge = dropReads(&readStart, incr);
				if (readStart == numquery) // no reads left
				{
					// I think this is ok as it means the last set of reads had the right trend
					break;
					// printf("Error: Could not find max age for opt!\n");
					// exit(0);
				}
			}
		} while (twice < 2); // DOUBLE CHECK THIS/Think of better way

		// Current stop gap... not great
		readStart = oldReadStart;
		nodePointer = oldNodePointer;
		nextNodeAge = oldNodeAge;

		printf("Maximum bound found, now finding optimum age.\n%lu of %lu reads left with maximum age of %.16f\n", numquery - readStart, numquery, nextNodeAge);
		// printf("The elapsed time for rough estimation is %.16f seconds\n", ( ((double) clock()) / CLOCKS_PER_SEC) - time2);
	}
	else
	{
		printf("Not enough reads for bound optimization, using all\n");
	}

	// printf("Remained merged assignments are ");
	// for (int z = readStart; z < numquery; z++)
	//{
	//	printf("%d ", assignments[usedReads[z]]);
	// }
	// printf("\n");

	fprintf(outfile, "reads_used=%lu\n", numquery - readStart);

	// Some bounds or fillers added
	p[0] = readStart;
	p[1] = nodePointer;
	p[2] = 0;
	invector[1] = nextNodeAge / 2;
	lowbound[1] = eh0;
	upbound[1] = nextNodeAge - eh0;
	nfun = 0;
	onDindic = 0; // Because of 2nd layer of optimization

	// double time3 = (double) clock()/CLOCKS_PER_SEC;

	// Decide how the inputs may need to change at some point I guess
	est_age_lik = Brent1D(invector, lowbound, upbound, 1, getlike_ages, p, 3);
	// printf("The elapsed time for age estimation is %.16f seconds\n", ( ((double) clock()) / CLOCKS_PER_SEC) - time3);

	est_age = invector[1];

	if (nextNodeAge - est_age < nextNodeAge / 100.0)
	{
		printf("Warning! Age estimate was near boundary of %.16f\n", nextNodeAge);
	}

	printf("Opt found, calculating confidence intervals\n");

	// Use fisher information to get rough confidence interval
	// Maybe give option for bootstrap confidence interval
	// double time4 = (double) clock()/CLOCKS_PER_SEC;
	confI = confidenceIntervalFisher(par, est_age, est_age_lik, readStart);
	// printf("The elapsed time for confidenceIntervalFisher is %.16f seconds\n", ( ((double) clock()) / CLOCKS_PER_SEC) - time4);

	// Also maybe incorporate more options for the confidence interval, not just 95%
	// Would need to be able to calculate the z-score from the user given value
	// confI = 1.96 / sqrt(-secD);

	printf("Estimated age is %.16f with likelihood %.16f and 95%% confidence interval [%.16f,%.16f]\n", est_age, est_age_lik, est_age - confI, est_age + confI);
	fprintf(outfile, "estimated_age=%.16f\n", est_age);
	fprintf(outfile, "likelihood=%.16f\n", est_age_lik);
	fprintf(outfile, "CI_lower=%.16f\n", est_age - confI);
	fprintf(outfile, "CI_upper=%.16f\n", est_age + confI);
	// printf("%.16f,%.16f\n", est_age, est_age_lik);

	// confidenceSearch(bounds, chiValue, maxAge, est_age, -est_age_lik, p);
	// printf("Estimated age is %.16f with likelihood %.16f and %.2f%% confidence interval [%.16f,%.16f]\n", est_age, est_age_lik, chiValue, bounds[0], bounds[2]);
	// printf("Confidence interval is [%.16f, %.16f] with likelihoods %.16f and %.16f\n", bounds[0], bounds[2], bounds[1], bounds[3]);

	if (allTrees)
	{
		for (unsigned long int i = 0; i < numTrees; i++)
		{
			// Some bounds or fillers added
			p[0] = readStart;
			p[1] = i;
			p[2] = 0;
			invector[1] = nextNodeAge / 2;
			lowbound[1] = eh0;
			upbound[1] = nextNodeAge - eh0;
			nfun = 0;
			onDindic = 0; // Because of 2nd layer of optimization

			// double time3 = (double) clock()/CLOCKS_PER_SEC;

			// Decide how the inputs may need to change at some point I guess
			est_age_lik = Brent1D(invector, lowbound, upbound, 1, getlike_ages_tree, p, 3);
			// printf("The elapsed time for age estimation is %.16f seconds\n", ( ((double) clock()) / CLOCKS_PER_SEC) - time3);

			est_age = invector[1];

			printf("Estimated age is %.16f with likelihood %.16f for tree %lu\n", est_age, est_age_lik, i);
		}
	}
}

int reassign_up(int move_node, int treeNum)
{
	int i;
	int assignmentChanged;
	double move_node_age = nodeages[treeNum][move_node] + bls[treeNum][move_node];

	// printf("Changing assignment of at least node %d in tree %d\n", move_node, treeNum);
	// printf("Original assignments: ");
	// for(i = 0; i < numquery; i++)
	//{
	//	printf("%d ", assignments[i]);
	// }
	// printf("\n");

	// assignmentChanged = getGFLPar(move_node, treeNum);

	// if(assignmentChanged == treeRoots[treeNum])
	//{
	//	return(-1);
	// }

	// printf("New assignments: ");
	for (i = 0; i < numquery; i++)
	{
		if (move_node_age >= assignAges[i])
		{
			assignments[i] = getGFLPar(assignments[i], treeAssign[i]);
			assignAges[i] = nodeages[treeAssign[i]][assignments[i]] + bls[treeAssign[i]][assignments[i]];
			if (assignments[i] == treeRoots[treeAssign[i]])
			{
				return (-1);
			}
		}
		// printf("%d ", assignments[i]);
	}
	// printf("\n");

	// indicate success
	return (0);
}

void storeOldOrder(int *oldAssign)
{
	for (unsigned long int i = 0; i < numquery; i++)
	{
		oldAssign[i] = assignments[i];
	}
}

void maximize_like_jointly_for_all2D_reassign(double **par, int allTrees)
{
	printf("Starting maximize jointly for all, reassign instead of dropping\n");
	int i, k, v, nfun, twice = 0, *oldAssign, reassign_success;
	double p[3], L1, L2, secD, confI;
	double lowbound[3], upbound[3], eh0 = 3e-8, nextNodeAge;
	double *invect_temp, *invectorL1, *invectorL2;
	double brent_invector[numquery+1], L1_test, L2_test;
	//printf("Max age of %.16f and incr of %.16f\n", totMaxAge, incr);

	invectorL1 = (double *)malloc(sizeof(double) * 3);
	invectorL2 = (double *)malloc(sizeof(double) * 3);
	oldAssign = (int *)malloc(sizeof(int) * numquery);
	usedReads = (int *)malloc(sizeof(int) * numquery);

	// these stay the same no matter
	p[0] = 0;
	p[1] = 0;
	p[2] = 0;
	lowbound[1] = eh0;
	nfun = 0;
	onDindic = 0; // Because of 2nd layer of optimization

	orderReads();
	// This nextNodeAge represents the max without dropping reads
	nextNodeAge = assignAges[usedReads[0]];

	 //printf("Assignments before age optimization:\n");
	 //for (unsigned long int i = 0; i < numquery; i++)
	 //{
	 //	printf("\t%d\t%d\t%d\n", usedReads[i], treeAssign[usedReads[i]], assignments[usedReads[i]]);
	 //}
	 //printf("\n");

	fprintf(outfile, "reads_used=%lu\n", numquery);

	// Some bounds or fillers added
	invectorL1[1] = nextNodeAge / 2;
	upbound[1] = nextNodeAge - eh0;
	// printf("Next node age is %.16f\n", nextNodeAge);
	L1 = Brent1D(invectorL1, lowbound, upbound, 1, getlike_ages, p, 3);
	//printf("L1 (%.16f, %.16f)\n", L1, invectorL1[1]);

	//brent_invector[numquery] = nextNodeAge/2;
	//for(int j = 0; j < numquery; j++)
	//{
	//	brent_invector[j] = nodeages[treeAssign[j]][assignments[j]] + bls[treeAssign[j]][assignments[j]]/2;
	//}
	//curAgeBound = nextNodeAge - eh0;
	//L1_test = minimize_brent(brent_invector, numquery+1, getlike_gamma_root_in_trifucation_sample_age_brent, 10000);
	//L1_test = GoldenSection(invectorL1, lowbound, upbound, 1, getlike_ages_brent, p, 3);
	//printf("\tvs L1_test (%.16f, %.16f)\n", L1_test, brent_invector[numquery]);

	//exit(0);

	//printf("L1 (%.16f, %.16f) vs L1_test (%.16f, %.16f)\n", L1, invectorL1[1], L1_test, brent_invector[numquery]);

	//exit(0);

	storeOldOrder(oldAssign);
	reassign_success = reassign_up(assignments[usedReads[0]], treeAssign[usedReads[0]]);

	if (numquery > 1)
	{
		while (1)
		{

			orderReads();
			// This nextNodeAge represents the max without dropping reads
			nextNodeAge = assignAges[usedReads[0]];
			// printf("Next node age is %.16f\n", nextNodeAge);
			// Some bounds or fillers added
			invectorL2[1] = nextNodeAge / 2;
			upbound[1] = nextNodeAge - eh0;
			curAgeBound = nextNodeAge - eh0;

			L2 = Brent1D(invectorL2, lowbound, upbound, 1, getlike_ages, p, 3);

			// printf("L1 (%.16f, %.16f) vs L2 (%.16f, %.16f)\n", L1, invectorL1[1], L2, invectorL2[1]);

			//brent_invector[numquery] = nextNodeAge/2;

			//for(int j = 0; j < numquery; j++)
			//{
			//	brent_invector[j] = nodeages[treeAssign[j]][assignments[j]] + bls[treeAssign[j]][assignments[j]]/2;
			//}

			//L2_test = minimize_brent(brent_invector, numquery+1, getlike_gamma_root_in_trifucation_sample_age_brent, 1000000);
			//printf("L2 (%.16f, %.16f) vs L2_test (%.16f, %.16f)\n", L2, invectorL2[1], L2_test, brent_invector[numquery]);

			// We are working with assuming L1 is better, store assignments/likelihoods when L2 is better
			if (L2 <= L1)
			{
				storeOldOrder(oldAssign);
				L1 = L2;
				invect_temp = invectorL1;
				invectorL1 = invectorL2;
				invectorL2 = invect_temp;
				invect_temp = NULL;
			}
			else
			{
				free(assignments);
				assignments = oldAssign;
				oldAssign = NULL;
				break;
			}

			reassign_success = reassign_up(assignments[usedReads[0]], treeAssign[usedReads[0]]);
			if (L1 > L2 || reassign_success == treeRoots[treeAssign[usedReads[0]]] || reassign_success == -1)
			{
				// Can no longer push reads up the tree if the new assignment is the root or if it failed
				free(assignments);
				assignments = oldAssign;
				oldAssign = NULL;
				if (reassign_success == -1 || reassign_success == treeRoots[treeAssign[usedReads[0]]])
				{
					printf("Opt stopped due to root of tree %d\n", treeAssign[usedReads[0]]);
				}
				break;
			}
		}
	}

	// double time3 = (double) clock()/CLOCKS_PER_SEC;

	if (nextNodeAge - invectorL1[1] < nextNodeAge / 100.0)
	{
		printf("Warning! Age estimate was near boundary of %.16f\n", nextNodeAge);
	}

	printf("Opt found, calculating confidence intervals\n");

	// Use fisher information to get rough confidence interval
	// Maybe give option for bootstrap confidence interval
	// double time4 = (double) clock()/CLOCKS_PER_SEC;
	confI = confidenceIntervalFisher(par, invectorL1[1], L1, 0);
	// printf("The elapsed time for confidenceIntervalFisher is %.16f seconds\n", ( ((double) clock()) / CLOCKS_PER_SEC) - time4);

	// Also maybe incorporate more options for the confidence interval, not just 95%
	// Would need to be able to calculate the z-score from the user given value
	// confI = 1.96 / sqrt(-secD);

	//printf("Assignments used in age estimation:\n");
	//for (unsigned long int i = 0; i < numquery; i++)
	//{
	//	printf("\t%d\t%d\t%d\n", usedReads[i], treeAssign[usedReads[i]], assignments[usedReads[i]]);
	//}
	//printf("\n");

	printf("Estimated age is %.16f with likelihood %.16f and 95%% confidence interval [%.16f,%.16f]\n", invectorL1[1], L1, invectorL1[1] - confI, invectorL2[1] + confI);
	fprintf(outfile, "estimated_age=%.16f\n", invectorL1[1]);
	fprintf(outfile, "likelihood=%.16f\n", L1);
	fprintf(outfile, "CI_lower=%.16f\n", invectorL1[1] - confI);
	fprintf(outfile, "CI_upper=%.16f\n", invectorL2[1] + confI);
	// printf("%.16f,%.16f\n", est_age, est_age_lik);

	// confidenceSearch(bounds, chiValue, maxAge, est_age, -est_age_lik, p);
	// printf("Estimated age is %.16f with likelihood %.16f and %.2f%% confidence interval [%.16f,%.16f]\n", est_age, est_age_lik, chiValue, bounds[0], bounds[2]);
	// printf("Confidence interval is [%.16f, %.16f] with likelihoods %.16f and %.16f\n", bounds[0], bounds[2], bounds[1], bounds[3]);

	// TO DO: REIMPLEMENT PER TREE VERSION
	// if(allTrees)
	//{
	//	for(int i = 0; i < numTrees; i++)
	//	{
	//		//Some bounds or fillers added
	//		p[0] = 0;
	//		p[1] = i;
	//		p[2] = 0;
	//		invector[1] = nextNodeAge/2;
	//		lowbound[1] = eh0;
	//		upbound[1] = nextNodeAge - eh0;
	//		nfun=0;
	//		onDindic=0;	// Because of 2nd layer of optimization

	//		//double time3 = (double) clock()/CLOCKS_PER_SEC;

	//		// Decide how the inputs may need to change at some point I guess
	//		est_age_lik = GoldenSection(invector,lowbound, upbound, 1, getlike_ages_tree, p, 3);
	//		//printf("The elapsed time for age estimation is %.16f seconds\n", ( ((double) clock()) / CLOCKS_PER_SEC) - time3);

	//		est_age = invector[1];

	//		printf("Estimated age is %.16f with likelihood %.16f for tree %d\n", est_age, est_age_lik, i);
	//	}
	//}

	free(invectorL1);
	free(invectorL2);
	free(oldAssign);
}

void maximize_like_jointly_for_all_noDrop2D(double **par, int allTrees)
{
	printf("Starting maximize jointly for all\n");
	int i, k, v, nfun, readStart = 0, nodePointer = 0, twice = 0, oldReadStart, oldNodePointer;
	double p[3], L1, L2, secD, confI;
	double invector[3], lowbound[3], upbound[3], eh0 = 3e-8, nextNodeAge, est_age, est_age_lik;
	usedReads = (int *)malloc(sizeof(int) * numquery);
	// double incr = totMaxAge/1000;

	// printf("Max age of %.16f and incr of %.16f\n", totMaxAge, incr);

	orderReads();

	// This nextNodeAge represents the max without dropping reads
	nextNodeAge = nodeages[treeAssign[usedReads[readStart]]][assignments[usedReads[readStart]]] + bls[treeAssign[usedReads[readStart]]][assignments[usedReads[readStart]]];

	printf("Max bound is %.16f\n", nextNodeAge);

	onDindic = 1;

	fprintf(outfile, "reads_used=%lu\n", numquery);

	// Some bounds or fillers added
	p[0] = readStart;
	p[1] = nodePointer;
	p[2] = 0;
	invector[1] = nextNodeAge / 2;
	lowbound[1] = eh0;
	upbound[1] = nextNodeAge - eh0;
	nfun = 0;
	onDindic = 0; // Because of 2nd layer of optimization

	// double time3 = (double) clock()/CLOCKS_PER_SEC;

	// Decide how the inputs may need to change at some point I guess
	est_age_lik = Brent1D(invector, lowbound, upbound, 1, getlike_ages, p, 3);
	// printf("The elapsed time for age estimation is %.16f seconds\n", ( ((double) clock()) / CLOCKS_PER_SEC) - time3);

	est_age = invector[1];

	if (nextNodeAge - est_age < nextNodeAge / 100.0)
	{
		printf("Warning! Age estimate was near boundary of %.16f\n", nextNodeAge);
	}

	printf("Opt found, calculating confidence intervals\n");

	// Use fisher information to get rough confidence interval
	// Maybe give option for bootstrap confidence interval
	// double time4 = (double) clock()/CLOCKS_PER_SEC;
	confI = confidenceIntervalFisher(par, est_age, est_age_lik, readStart);
	// printf("The elapsed time for confidenceIntervalFisher is %.16f seconds\n", ( ((double) clock()) / CLOCKS_PER_SEC) - time4);

	// Also maybe incorporate more options for the confidence interval, not just 95%
	// Would need to be able to calculate the z-score from the user given value
	// confI = 1.96 / sqrt(-secD);

	printf("Estimated age is %.16f with likelihood %.16f and 95%% confidence interval [%.16f,%.16f]\n", est_age, est_age_lik, est_age - confI, est_age + confI);
	fprintf(outfile, "estimated_age=%.16f\n", est_age);
	fprintf(outfile, "likelihood=%.16f\n", est_age_lik);
	fprintf(outfile, "CI_lower=%.16f\n", est_age - confI);
	fprintf(outfile, "CI_upper=%.16f\n", est_age + confI);

	// confidenceSearch(bounds, chiValue, maxAge, est_age, -est_age_lik, p);
	// printf("Estimated age is %.16f with likelihood %.16f and %.2f%% confidence interval [%.16f,%.16f]\n", est_age, est_age_lik, chiValue, bounds[0], bounds[2]);
	// printf("Confidence interval is [%.16f, %.16f] with likelihoods %.16f and %.16f\n", bounds[0], bounds[2], bounds[1], bounds[3]);

	if (allTrees)
	{
		for (unsigned long int i = 0; i < numTrees; i++)
		{
			// Some bounds or fillers added
			p[0] = readStart;
			p[1] = i;
			p[2] = 0;
			invector[1] = nextNodeAge / 2;
			lowbound[1] = eh0;
			upbound[1] = nextNodeAge - eh0;
			nfun = 0;
			onDindic = 0; // Because of 2nd layer of optimization

			// double time3 = (double) clock()/CLOCKS_PER_SEC;

			// Decide how the inputs may need to change at some point I guess
			est_age_lik = Brent1D(invector, lowbound, upbound, 1, getlike_ages_tree, p, 3);
			// printf("The elapsed time for age estimation is %.16f seconds\n", ( ((double) clock()) / CLOCKS_PER_SEC) - time3);

			est_age = invector[1];

			printf("Estimated age is %.16f with likelihood %.16f for tree %lu\n", est_age, est_age_lik, i);
		}
	}
}

void age_like_distribution_jointly_for_all2D_upperLimit(double **par, double topAge)
{
	printf("Creating age likelihood distribution for all\n");
	int i, k, v, nfun, readStart = 0;
	int numDrops = 0;
	double p[3], L1, L2, ageIncr, nextNodeAge, tempAge, tempAge2;
	double invector[3], lowbound[3], upbound[3], eh0 = 1e-8, age_like;
	usedReads = (int *)malloc(sizeof(int) * numquery);

	orderReads();

	// An important question here is why was the testAge for below so off?
	// nodePointer = assignments[usedReads[readStart]];
	// nextNodeAge = nodeages[nodeOrder[nodePointer]] + bls[nodeOrder[nodePointer]];
	nextNodeAge = nodeages[treeAssign[usedReads[readStart]]][assignments[usedReads[readStart]]] + bls[treeAssign[usedReads[readStart]]][assignments[usedReads[readStart]]];

	// printf("First node age as old %.16f and by assigning nodePointer %.16f\n")

	// printf("First node age %.16f\n", nextNodeAge);
	// printf("First node age %.16f from assignment %d with node age %.16f and branch length %.16f that has order %d\n", nextNodeAge, nodePointer, nodeages[nodeOrder[nodePointer]], bls[nodeOrder[nodePointer]], nodeOrder[nodePointer]);

	ageIncr = topAge / 10000.0; // make this an option later

	testAge = eh0;

	onDindic = 1;

	// assignmentMode of 0 means single assignment given
	age_like = 0.0;
	while (testAge <= totMaxAge && readStart < numquery && testAge <= topAge)
	{
		// Trimmed because you can't  actually compare the likelihoods, just the trends.
		if (nextNodeAge <= testAge)
		{
			// This is probably unneeded due to the while loop condition
			if (nextNodeAge >= totMaxAge)
			{
				break;
			}
			nextNodeAge = dropReads(&readStart, 0.01); // Not sure if this is correct, may want to come back/make two dropReads versions
			ageIncr = topAge / 100.0;
			numDrops++;

			testAge = eh0;
		}
		printf("Testing age %.16f with dropped rounds %d with nextNodeAge %.16f and read start at %d:\n", testAge, numDrops, nextNodeAge, readStart);
		for (i = readStart; i < numquery; i++)
		{
			p[0] = usedReads[i];
			p[1] = treeAssign[usedReads[i]];
			p[2] = assignments[usedReads[i]];

			nfun = 0;
			invector[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]] / 2.0;
			lowbound[1] = eh0;
			upbound[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]] - eh0;

			L2 = Brent1D(invector, lowbound, upbound, 1, getlike_gamma_root_in_trifurcation_testAge, p, 3);

			// printf("\tRead %d assigned to %d with likelihood %lf with root placement of %.16f\n", usedReads[i], assignments[usedReads[i]], L2, nodeages[assignments[usedReads[i]]]+invector[1]);

			age_like += L2; // sum of log likelihoods
		}
		printf("\tLikelihood: %.16f\n", age_like);
		age_like = 0.0;
		testAge += ageIncr;
	}
}

void age_like_distribution_jointly_for_all2D(double **par)
{
	printf("Creating age likelihood distribution for all\n");
	int i, k, v, nfun, readStart = 0;
	int numDrops = 0;
	double p[3], L1, L2, ageIncr, nextNodeAge, tempAge, tempAge2;
	double invector[3], lowbound[3], upbound[3], eh0 = 1e-8, age_like;
	usedReads = (int *)malloc(sizeof(int) * numquery);

	orderReads();

	// An important question here is why was the testAge for below so off?
	// nodePointer = assignments[usedReads[readStart]];
	// nextNodeAge = nodeages[nodeOrder[nodePointer]] + bls[nodeOrder[nodePointer]];
	nextNodeAge = nodeages[treeAssign[usedReads[readStart]]][assignments[usedReads[readStart]]] + bls[treeAssign[usedReads[readStart]]][assignments[usedReads[readStart]]];

	// printf("First node age as old %.16f and by assigning nodePointer %.16f\n")

	// printf("First node age %.16f\n", nextNodeAge);
	// printf("First node age %.16f from assignment %d with node age %.16f and branch length %.16f that has order %d\n", nextNodeAge, nodePointer, nodeages[nodeOrder[nodePointer]], bls[nodeOrder[nodePointer]], nodeOrder[nodePointer]);

	ageIncr = nextNodeAge / 100.0; // make this an option later

	testAge = eh0;

	onDindic = 1;

	// assignmentMode of 0 means single assignment given
	age_like = 0.0;
	while (testAge <= totMaxAge && readStart < numquery)
	{
		// Trimmed because you can't  actually compare the likelihoods, just the trends.
		if (nextNodeAge <= testAge)
		{
			// This is probably unneeded due to the while loop condition
			if (nextNodeAge >= totMaxAge)
			{
				break;
			}
			nextNodeAge = dropReads(&readStart, 0.01); // Not sure if this is correct, may want to come back/make two dropReads versions
			ageIncr = nextNodeAge / 100.0;
			numDrops++;

			testAge = eh0;
		}
		printf("Testing age %.16f with dropped rounds %d with nextNodeAge %.16f and read start at %d:\n", testAge, numDrops, nextNodeAge, readStart);
		for (i = readStart; i < numquery; i++)
		{
			p[0] = usedReads[i];
			p[1] = treeAssign[usedReads[i]];
			p[2] = assignments[usedReads[i]];

			nfun = 0;
			invector[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]] / 2.0;
			lowbound[1] = eh0;
			upbound[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]] - eh0;

			L2 = Brent1D(invector, lowbound, upbound, 1, getlike_gamma_root_in_trifurcation_testAge, p, 3);

			// printf("\tRead %d assigned to %d with likelihood %lf with root placement of %.16f\n", usedReads[i], assignments[usedReads[i]], L2, nodeages[assignments[usedReads[i]]]+invector[1]);

			age_like += L2; // sum of log likelihoods
		}
		printf("\tLikelihood: %.16f\n", age_like);
		age_like = 0.0;
		testAge += ageIncr;
	}
}

void read_branch_like_dist(double **par)
{
	printf("Getting read likelihood distribution over assigned branch\n");
	int i, k, v, nfun;
	double p[3], L1, L2, increment;
	double invector[3], lowbound[3], upbound[3], eh0 = 3e-8;

	onDindic = 1;

	// assignmentMode of 0 means single assignment given
	for (i = 0; i < numquery; i++)
	{
		p[0] = i;
		p[1] = treeAssign[i];
		p[2] = assignments[i];

		// printf("Sequence %d assigned to tree %d and node %d\n", i, (int)p[1], (int)p[2]);
		rooted = eh0;
		increment = bls[treeAssign[i]][assignments[i]] / 100.0;

		while (rooted < bls[treeAssign[i]][assignments[i]] - eh0)
		{

			invector[1] = 0.5;
			lowbound[1] = eh0;
			upbound[1] = 1.0 - eh0;

			L2 = Brent1D(invector, lowbound, upbound, 1, getlike_gamma_root_in_trifurcation_testRoot, p, 3);

			printf("%d,%d,%d,%.16f,%.16f\n", i, treeAssign[i], assignments[i], rooted, L2);
			rooted += increment;
		}
	}
}

void maximize_like_jointly_for_all2D_and_LR()
{
}

void maximize_like_jointly_for_all_gamma_plus_zero()
{
}

void postrior_assignment_for_all_gamma_plus_zero()
{
}

// Free globals
void freeData()
{
	// Things to free that are only single *
	free(assignments);
	free(startpos);
	free(assignAges);
	free(treeAssign);
	free(usedReads);

	free(numbases);
	free(treeRoots);
	free(maxAges);

	// More complex frees
	for (unsigned long int i = 0; i < numquery; i++)
	{
		free(QUERYDATA[i]);
		for (unsigned long int j = 0; j < readlength[i]; j++)
		{
			free(readlike[i][j]);
		}
		free(readlike[i]);
	}

	for (unsigned long int i = 0; i < numTrees; i++)
	{
		free(pi[i]);
		free(par[i]);
		for (unsigned long int j = 0; j < numseq[i]; j++)
		{
			free(DATA[i][j]);
		}
		free(DATA[i]);
		free(bls[i]);
		free(nodeages[i]);
		free(FRACLIKE[i]);
		free(statevector[i]);
		free(nodeOrder[i]);
		free(RRVAL[i]);
		free(RRVEC[i]);
		free(LRVEC[i]);
	}

	free(numseq);
	free(QUERYDATA);
	free(readlike);
	free(readlength);
	free(pi);
	free(par);
	free(DATA);
	free(bls);
	free(nodeages);
	free(FRACLIKE);
	free(statevector);
	free(nodeOrder);
	free(RRVAL);
	free(RRVEC);
	free(LRVEC);
}

static void print_usage(const char *prog)
{
	printf("Usage: %s -a ASSIGN_FILE -q QUERY_FILE -l LIK_DIR -p PARAM_DIR \\\n", prog);
	printf("              -n NUM_TREES -m MODE -e ERROR_FILE \\\n");
	printf("              -A ALL_TREES -c COMPARE_AGE -r MERGE_READS -k REASSIGN \\\n");
	printf("              -o OUT_PREFIX [-R] [-C COVERAGE] [-t THREADS]\n");
	printf("  -a/--assign-file:  Path to sample assignment file.\n");
	printf("  -q/--query-file:   Path to sample alignment file.\n");
	printf("  -l/--lik-dir:      Path to likelihood directory.\n");
	printf("  -p/--param-dir:    Path to parameter/tree directory.\n");
	printf("  -n/--num-trees:    Number of trees.\n");
	printf("  -m/--mode:         Estimation mode (0-9):\n");
	printf("                       0=maximize_like_seperately_for_all2D_Print\n");
	printf("                       1=maximize_like_seperately_for_all2D\n");
	printf("                       2=likelihoodratiotest_for_all\n");
	printf("                       3=maximize_like_jointly_for_all2D\n");
	printf("                       4=age_like_distribution_jointly_for_all2D\n");
	printf("                       5=age_like_distribution_jointly_for_all2D_upperLimit\n");
	printf("                       6=maximize_like_jointly_for_all_noDrop2D\n");
	printf("                       7=read_branch_like_dist\n");
	printf("                       8=readContour\n");
	printf("                       9=maximize_like_jointly_for_all2D_reassign\n");
	printf("  -e/--error-file:   Path to error profile file.\n");
	printf("  -A/--all-trees:    0/1 - estimate age for each tree individually or jointly.\n");
	printf("  -c/--compare-age:  Age for Likelihood Ratio Test (0 if not testing or use modern).\n");
	printf("  -r/--merge-reads:  0/1 - Merge reads assigned to the same edge.\n");
	printf("  -k/--reassign:     0=keep inputted, 1=full (bestAssignment), 2=tronko.\n");
	printf("  -o/--out-prefix:   Output file prefix (required). Creates <prefix>.summary.\n");
	printf("  -R/--write-reads:  (optional) Write merged read sequences to <prefix>.reads.\n");
	printf("  -C/--coverage:     (optional) Minimum coverage to keep merged read. Default '0.05f'.\n");
	printf("                     Examples: '0.05f' (5%% fraction), '50b' (50 base pairs).\n");
	printf("  -t/--threads:      (optional) Number of OpenMP threads. Default 1.\n");
}

int main(int argc, char *argv[])
{
	// time for whole program:
	// double time1 = (double) clock()/CLOCKS_PER_SEC;
	// srand(time(NULL) ^ getpid()); // Combine time and process ID for uniqueness
	srand(0);

	double L, compareAge;

	char assignfile[500], fraclikefile[500], querydatafile[500], referencedatafile[500], errorfile[500], strTree[100], tempFileName[500];
	int mode, allTrees, reassign_mode;
	char out_prefix[500] = {0};
	char num_trees_str[64] = {0};

	static struct option long_options[] = {
		{"assign-file",  required_argument, 0, 'a'},
		{"query-file",   required_argument, 0, 'q'},
		{"lik-dir",      required_argument, 0, 'l'},
		{"param-dir",    required_argument, 0, 'p'},
		{"num-trees",    required_argument, 0, 'n'},
		{"mode",         required_argument, 0, 'm'},
		{"error-file",   required_argument, 0, 'e'},
		{"out-prefix",   required_argument, 0, 'o'},
		{"all-trees",    required_argument, 0, 'A'},
		{"compare-age",  required_argument, 0, 'c'},
		{"merge-reads",  required_argument, 0, 'r'},
		{"reassign",     required_argument, 0, 'k'},
		{"coverage",     required_argument, 0, 'C'},
		{"threads",      required_argument, 0, 't'},
		{"write-reads",  no_argument,       0, 'R'},
		{0, 0, 0, 0}
	};
	const char *short_opts = "a:q:l:p:n:m:e:o:A:c:r:k:C:t:R";

	int opt, option_index = 0;
	int seen_a=0, seen_q=0, seen_l=0, seen_p=0, seen_n=0,
	    seen_m=0, seen_e=0, seen_o=0, seen_A=0, seen_c=0, seen_r=0, seen_k=0;
	int num_threads = 1;
	int write_reads_flag = 0;

	while ((opt = getopt_long(argc, argv, short_opts, long_options, &option_index)) != -1) {
		switch (opt) {
			case 'a': sprintf(assignfile, "%s", optarg);        seen_a=1; break;
			case 'q': sprintf(querydatafile, "%s", optarg);     seen_q=1; break;
			case 'l': sprintf(fraclikefile, "%s", optarg);      seen_l=1; break;
			case 'p': sprintf(referencedatafile, "%s", optarg); seen_p=1; break;
			case 'n': sprintf(num_trees_str, "%s", optarg);     seen_n=1; break;
			case 'm': mode = atoi(optarg);                       seen_m=1; break;
			case 'e': sprintf(errorfile, "%s", optarg);         seen_e=1; break;
			case 'o': sprintf(out_prefix, "%s", optarg);        seen_o=1; break;
		case 'R': write_reads_flag = 1;                      break;
			case 'A': allTrees = atoi(optarg);                  seen_A=1; break;
			case 'c': compareAge = atof(optarg);                seen_c=1; break;
			case 'r': toMerge = atoi(optarg);                   seen_r=1; break;
			case 'k':
			reassign_mode = atoi(optarg); seen_k=1;
			if (reassign_mode < 0 || reassign_mode > 2) {
				fprintf(stderr, "Error: --reassign must be 0, 1, or 2.\n");
				exit(1);
			}
			break;
			case 'C': {
				char *threshold_str = optarg;
				int len = strlen(threshold_str);
				char last_char = 'f';
				if (len > 0 && (threshold_str[len-1] == 'f' || threshold_str[len-1] == 'b')) {
					last_char = threshold_str[len-1];
					threshold_str[len-1] = '\0';
				}
				char *endptr_thresh;
				merge_coverage_threshold = strtod(threshold_str, &endptr_thresh);
				if (threshold_str == endptr_thresh || *endptr_thresh != '\0' || merge_coverage_threshold < 0) {
					fprintf(stderr, "Error: Invalid merge coverage threshold format: %s\n", optarg);
					fprintf(stderr, "Use format like '0.05f' or '50b'.\n");
					exit(1);
				}
				if (last_char == 'b') {
					merge_coverage_mode = MERGE_MODE_BP;
					if (merge_coverage_threshold != floor(merge_coverage_threshold)) {
						fprintf(stderr, "Warning: Base pair coverage threshold provided as non-integer. Using floor value %d.\n", (int)floor(merge_coverage_threshold));
						merge_coverage_threshold = floor(merge_coverage_threshold);
					}
				} else {
					merge_coverage_mode = MERGE_MODE_FRACTION;
					if (merge_coverage_threshold > 1.0) {
						fprintf(stderr, "Error: Fractional merge coverage threshold cannot be greater than 1.0.\n");
						exit(1);
					}
				}
				break;
			}
			case 't':
				num_threads = atoi(optarg);
				if (num_threads < 1) {
					fprintf(stderr, "Error: --threads must be >= 1.\n");
					exit(1);
				}
				break;
			default:
				print_usage(argv[0]); exit(1);
		}
	}

	if (!seen_a || !seen_q || !seen_l || !seen_p || !seen_n ||
	    !seen_m || !seen_e || !seen_o || !seen_A || !seen_c || !seen_r || !seen_k) {
		fprintf(stderr, "Error: missing required flag(s)\n");
		if (!seen_a) fprintf(stderr, "  missing: -a/--assign-file\n");
		if (!seen_q) fprintf(stderr, "  missing: -q/--query-file\n");
		if (!seen_l) fprintf(stderr, "  missing: -l/--lik-dir\n");
		if (!seen_p) fprintf(stderr, "  missing: -p/--param-dir\n");
		if (!seen_n) fprintf(stderr, "  missing: -n/--num-trees\n");
		if (!seen_m) fprintf(stderr, "  missing: -m/--mode\n");
		if (!seen_e) fprintf(stderr, "  missing: -e/--error-file\n");
		if (!seen_o) fprintf(stderr, "  missing: -o/--out-prefix\n");
		if (!seen_A) fprintf(stderr, "  missing: -A/--all-trees\n");
		if (!seen_c) fprintf(stderr, "  missing: -c/--compare-age\n");
		if (!seen_r) fprintf(stderr, "  missing: -r/--merge-reads\n");
		if (!seen_k) fprintf(stderr, "  missing: -k/--reassign\n");
		print_usage(argv[0]); exit(1);
	}

	{
		char summary_path[512], reads_path[512];
		snprintf(summary_path, sizeof(summary_path), "%s.summary", out_prefix);
		outfile = fopen(summary_path, "w");
		if (!outfile) { perror("Cannot open summary file"); exit(1); }

		if (write_reads_flag) {
			snprintf(reads_path, sizeof(reads_path), "%s.reads", out_prefix);
			readsfile = fopen(reads_path, "w");
			if (!readsfile) { perror("Cannot open reads file"); exit(1); }
		}
	}

#ifdef _OPENMP
	omp_set_num_threads(num_threads);
#endif
	printf("num_threads: %d\n", num_threads);

	// New checks for reading in files, following guidelines from get_frac_like.c
	// Get number of trees and check its a valid number
	// Reference: https://stackoverflow.com/questions/26080829/detecting-strtol-failure

	const char *nptr = num_trees_str;
	char *endptr = NULL;
	numTrees = 0;

	// reset errno to 0 before call for validation of number
	errno = 0;
	numTrees = strtol(nptr, &endptr, 10);

	/* test return to number and errno values,  */
	if (nptr == endptr || (errno == ERANGE && numTrees == LONG_MIN) || (errno == ERANGE && numTrees == LONG_MAX) || errno == EINVAL || (errno != 0 && numTrees == 0) || (errno == 0 && nptr && *endptr != 0))
	{
		printf("Error in reading number of trees. Please double check. If you think this is an error, please contact maya_lemmon-kishi@berkeley.edu\n");
		exit(0);
	}

	// printf("Number of trees being read %lu\n", numTrees);

	fprintf(outfile, "assign_file=%s\n", assignfile);
	fprintf(outfile, "query_file=%s\n", querydatafile);
	fprintf(outfile, "lik_dir=%s\n", fraclikefile);
	fprintf(outfile, "param_dir=%s\n", referencedatafile);
	fprintf(outfile, "num_trees=%lu\n", numTrees);
	fprintf(outfile, "mode=%d\n", mode);
	fprintf(outfile, "error_file=%s\n", errorfile);
	fprintf(outfile, "all_trees=%d\n", allTrees);
	fprintf(outfile, "compare_age=%f\n", compareAge);
	fprintf(outfile, "merge_reads=%d\n", toMerge);
	fprintf(outfile, "reassign_mode=%d\n", reassign_mode);
	fprintf(outfile, "coverage_threshold=%f\n", merge_coverage_threshold);
	fprintf(outfile, "coverage_mode=%s\n", merge_coverage_mode == MERGE_MODE_BP ? "bp" : "fraction");
	fprintf(outfile, "num_threads=%d\n", num_threads);

	numseq = (unsigned long int *)malloc(sizeof(unsigned long int) * numTrees);

	// Check if input directories exists
	// sprintf(likeDir, "%s", argv[3]);
	struct stat s;
	int err = stat(fraclikefile, &s);
	if (-1 == err)
	{
		if (ENOENT == errno)
		{
			printf("Likelihood directory does not exist. Please double check input\n");
			exit(0);
		}
		else
		{
			perror("Error in stat");
			exit(1);
		}
	}
	else
	{
		if (!S_ISDIR(s.st_mode))
		{
			printf("Likelihood path is not a directory. Please double check input\n");
			exit(0);
		}
	}

	// sprintf(treeDir, "%s", argv[4]);
	err = stat(referencedatafile, &s);
	if (-1 == err)
	{
		if (ENOENT == errno)
		{
			printf("Tree directory does not exist. Please double check input\n");
			exit(0);
		}
		else
		{
			perror("Error in stat");
			exit(1);
		}
	}
	else
	{
		if (!S_ISDIR(s.st_mode))
		{
			printf("Tree path is not a directory. Please double check input\n");
			exit(0);
		}
	}

	read_data(assignfile, fraclikefile, querydatafile, referencedatafile, errorfile, mode, reassign_mode);

	fprintf(outfile, "input_reads=%lu\n", numquery);

	// exit(0);

	if (mode == 0)
		maximize_like_seperately_for_all2D_Print(par);
	else if (mode == 1)
		maximize_like_seperately_for_all2D(par);
	else if (mode == 2)
		likelihoodratiotest_for_all(par, compareAge);
	else if (mode == 3)
		maximize_like_jointly_for_all2D(par, allTrees);
	else if (mode == 4)
		age_like_distribution_jointly_for_all2D(par);
	else if (mode == 5)
		age_like_distribution_jointly_for_all2D_upperLimit(par, 0.005);
	else if (mode == 6)
		maximize_like_jointly_for_all_noDrop2D(par, allTrees);
	else if (mode == 7)
		read_branch_like_dist(par);
	else if (mode == 8)
		readContour(par);
	else if (mode == 9)
		maximize_like_jointly_for_all2D_reassign(par, allTrees);
	else
		printf("Invalid mode. Valid modes are 0-9.\n");

	if (outfile)   fclose(outfile);
	if (readsfile) fclose(readsfile);

	freeNRinits(2);
	freetreememmory();
	freeData();

	// printf("The elapsed time for whole program is %.16f seconds\n", ( ((double) clock()) / CLOCKS_PER_SEC) - time1);
}
