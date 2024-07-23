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

//For testing, remove!
#include <time.h>

#define MINBL 0.0000001
#define MAXBL 5.0
#define NUMCAT 4
#define INF DBL_MAX
#define VERBOSE 0
#define testmax(a,b) \
	({ __typeof__ (a) _a = (a); \
	 __typeof__ (b) _b = (b); \
	 _a > _b ? _a : _b; })

#define INITIAL_BUFFER_SIZE 1000  // Initial buffer size, adjust as needed

//assignAges holds the max possible age for the assignment node of a read! So nodeage + bls of assignment node for that read
//double LRVEC[4][4], RRVEC[4][4], RRVAL[4], PMAT[3][NUMCAT][4][4];
double **LRVEC, **RRVEC, **RRVAL, PMAT[3][NUMCAT][4][4];	//I think
double **statevector, *****FRACLIKE, **nodeages, **bls, ***readlike, testAge, **pi, **par, *maxAges, totMaxAge, *assignAges, errorTest;
int numbase, numquery, queryagesknown, ***DATA, **QUERYDATA, *assignments, *readlength, *startpos, **nodeOrder, *usedReads, *treeAssign, *usedTrees;
int onDindic = 0;//hack to avoid passing this indicator around
int *numseq;
long numTrees;

int tip,comma=0; /*globals used to read in the tree. Old code - don't ask.*/


FILE *infile, *outfile;

//Adding tree information so that we can do node uncertainity!!
struct node {
	int up[2];	// for storing children node number, assuming binary. -1 means no children (current node is leaf)
	int down;	// for storing parent node number
	double bl;	//Scaled branch length to parent
};

// initializing tree into global
struct node *tree;


/*subfunction needed by ‘getclade*/
void linknodes(int i,int j,int nodee) /*linking i down to nodee and j down to nodee*/
{
	tree[nodee].up[0]=j;
	tree[nodee].up[1]=i;
	tree[i].down=nodee;
	tree[j].down=nodee;
}

/*subfunction needed by ‘getclade*/
/** What does this subfunction do exactly? **/
// is tip used to track number of leaves observed while reading in the file?
// 	if so, should there be a check that the number of tips cannot be larger than numleaves?
int specsearch(int numleaves)
{
	char ch;
	int i=1;
	ch = fgetc(infile);
	if ((ch!=')')&&(ch!='(')&&(ch!=',')&&(ch!=' ')&&(ch!='\t')&&(ch!='\n')&&(ch!=EOF)){
		ungetc(ch, infile);
		fscanf(infile,"%d",&tip);
		while ((ch=fgetc(infile))!=':'&&(i<10))
			i++;
		tree[tip+numleaves-2].up[0]=-1;
		tree[tip+numleaves-2].up[1]=-1;
		while ((ch=(fgetc(infile)))==' ');
		ungetc(ch,infile);
		fscanf(infile,"%lf",&tree[tip+numleaves-2].bl);
		/*printf("\nbranchlength of node %i =%f",tip+numleaves-1,tree[tip+numleaves-2].bl);*/
		return 1;
	}
	else {
		ungetc(ch, infile);
		return 0;
	}
}

/*subfunction needed by ‘getclade*/
/** Gets node number... but how? **/
int getnodenumb()
{
	char c;
	int i,j=0;
	fpos_t position;
	i=0;
	fgetpos(infile, &position);
	do{
		c=fgetc(infile);
		if (c==',')
			i++;
		if (c=='(')
			j=j-1;
		if (c==')')
			j++;
	} while ((j<0)&&(c!=EOF));
	fsetpos(infile,&position);
	return (i+comma+1);
}

/*some old code for reading a Newick tree*/
int getclade(int numleaves)
{
	int n1, n2, n3;
	char ch;

	do{
		if (specsearch(numleaves)==1){
			/*tip++;*/
			return tip+(numleaves-1);
		}
		ch = fgetc(infile);
		if (ch==','){
			comma++;
		}
		if (ch==')'){
			if ((ch=fgetc(infile))!=':'){
				ungetc(ch,infile);
			}
			else {
				do{
					ch=(fgetc(infile));
				}while((ch=='\n')||(ch==' '));
				ungetc(ch,infile);
				fscanf(infile,"%lf",&tree[n3-1].bl);
			}
			// returns root node of the (sub)tree
			return n3;
		}
		if (ch=='('){
			n3=getnodenumb();
			// Recursive function down each subtree
			n1=getclade(numleaves);
			n2=getclade(numleaves);
			linknodes(n1-1,n2-1,n3-1);
		}
	} while (ch!=';');

	return -1;
}

void allocatetreememmory(int numleaves)
{
	int i;

	// Number of nodes based off full binary tree
	tree=malloc((numleaves*2-1)*(sizeof(struct node)));
}

void freetreememmory(void)
{
	free(tree);
}

/*old code for printing a tree*/
void printtree(int numleaves, int root)

{
	int i;

	printf("\nPRINTING TREE\n");
	for (i=0; i<2*numleaves-1; i++)
	{
		if (tree[i].up[0] != -1)
			printf("Node %i: up: (%i, %i) down: %i",i,tree[i].up[0],tree[i].up[1],tree[i].down);
		else
			printf(" Node %i (leaf): up: (%i, %i) down: %i",i,tree[i].up[0],tree[i].up[1],tree[i].down);
		if (i != root)
			printf(" (bl: %f)\n",tree[i].bl);
		else printf(" (root)\n");
	}
}

/*Code getting node to leaf length assuming ultrametric (all paths the same length) tree*/
double getMaxAge(int curNode)
{
	if(tree[curNode].up[0] == -1)
	{
		return(tree[curNode].bl);
	}
	else
	{
		return(tree[curNode].bl + getMaxAge(tree[curNode].up[0]));
	}
}

//get gfl node number of children when given ratePlacer assignment node encoding
void getGFLChildren(int curNode, int childNodes[3], int treeNum)
{
	// leaf, no children of node to check assignment
	if(curNode < numseq[treeNum])
	{
		childNodes[1] = -1;
		childNodes[2] = -1;
	}
	else	//non-leaf
	{
		// convert curNode to gfl node
		curNode = curNode - numseq[treeNum];

		if(tree[tree[curNode].up[0]].up[0] == -1)
		{
			//Children are leaf
			childNodes[1] = tree[curNode].up[0] - numseq[treeNum] + 1;
		}
		else
		{
			//children are internal
			childNodes[1] = tree[curNode].up[0] + numseq[treeNum];
		}

		if(tree[tree[curNode].up[1]].up[0] == -1)
		{
			//Children are leaf
			childNodes[2] = tree[curNode].up[1] - numseq[treeNum] + 1;
		}
		else
		{
			//children are internal
			childNodes[2] = tree[curNode].up[1] + numseq[treeNum];
		}
	}	
}

//get gfl nodes of parent and sibling when given ratePlacer assignment node
//ASSUMES curNode IS A LEAF NODE
void getGFLParSib(int curNode, int parSib[3], int treeNum)
{
	printf("%d converted to ", curNode);
	//convert leaf curNode to gfl node
	curNode = curNode + numseq[treeNum] - 1;
	printf("%d\n", curNode);

	//get parent
	parSib[1] = tree[curNode].down + numseq[treeNum];

	printf("Parent is %d, converted to %d\n", tree[curNode].down, parSib[1]);

	//get sibling
	int par = tree[curNode].down;

	printf("Children of parent are %d and %d\n", tree[par].up[0], tree[par].up[1]);

	if(tree[par].up[0] == curNode)
	{
		//parSib[2] = tree[par].up[1] - numseq + 1;	
		if(tree[tree[par].up[1]].up[0] == -1)
		{
			//Children are leaf
			parSib[2] = tree[par].up[1] - numseq[treeNum] + 1;
		}
		else
		{
			//children are internal
			parSib[2] = tree[par].up[1] + numseq[treeNum];
		}
		
	}
	else
	{
		//parSib[2] = tree[par].up[0] - numseq + 1;
		if(tree[tree[par].up[0]].up[0] == -1)
		{
			//Children are leaf
			parSib[2] = tree[par].up[0] - numseq[treeNum] + 1;
		}
		else
		{
			//children are internal
			parSib[2] = tree[par].up[0] + numseq[treeNum];
		}

	}

	printf("Testing also %d %d\n", parSib[1], parSib[2]);
}

int getGFLPar(int curNode, int treeNum)
{
	if(curNode < numseq[treeNum])
	{
		//leaf node
		curNode = curNode + numseq[treeNum] - 1;	
	}
	else
	{
		//internal node
		curNode = curNode - numseq[treeNum];
	}

	//Parent is always internal node
	return(tree[curNode].down + numseq[treeNum]);
}

//******************FUNCTIONS TO DEAL WITH THE GAMMA DISTRIBUION********************************************
double LnGamma (double alpha)
{
	/* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.
	   Stirling's formula is used for the central polynomial part of the procedure.
	   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
	   Communications of the Association for Computing Machinery, 9:684
	   */
	double x=alpha, f=0, z;

	if (x<7) {
		f=1;  z=x-1;
		while (++z<7)  f*=z;
		x=z;   f=-log(f);
	}
	z = 1/(x*x);
	return  f + (x-0.5)*log(x) - x + .918938533204673
		+ (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
				+.083333333333333)/x;
}

double IncompleteGamma (double x, double alpha, double ln_gamma_alpha)
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
	double p=alpha, g=ln_gamma_alpha;
	double accurate=1e-20, overflow=1e30;
	double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];

	if (x==0) return (0);
	if (x<0 || p<=0) return (-1);

	factor=exp(p*log(x)-x-g);   
	if (x>1 && x>=p) goto l30;
	/* (1) series expansion */
	gin=1;  term=1;  rn=p;
l20:
	rn++;
	term*=x/rn;   gin+=term;

	if (term > accurate) goto l20;
	gin*=factor/p;
	goto l50;
l30:
	/* (2) continued fraction */
	a=1-p;   b=a+x+1;  term=0;
	pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
	gin=pn[2]/pn[3];
l32:
	a++;  b+=2;  term++;   an=a*term;
	for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
	if (pn[5] == 0) goto l35;
	rn=pn[4]/pn[5];   dif=fabs(gin-rn);
	if (dif>accurate) goto l34;
	if (dif<=accurate*rn) goto l42;
l34:
	gin=rn;
l35:
	for (i=0; i<4; i++) pn[i]=pn[i+2];
	if (fabs(pn[4]) < overflow) goto l32;
	for (i=0; i<4; i++) pn[i]/=overflow;
	goto l32;
l42:
	gin=1-factor*gin;

l50:
	/*printf("Incompletegamma got %f %f %f and returned %f\n",x,  alpha,  ln_gamma_alpha,gin);*/
	//printf("");
	return (gin);
}

double PointNormal (double prob)
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
	double a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245;
	double a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495;
	double b2=.531103462366, b3=.103537752850, b4=.0038560700634;
	double y, z=0, p=prob, p1;

	p1 = (p<0.5 ? p : 1-p);
	if (p1<1e-20) return (-9999);

	y = sqrt (log(1/(p1*p1)));
	z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
	return (p<0.5 ? -z : z);
}



double PointChi2 (double prob, double v)
{
	/* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
	   returns -1 if in error.   0.000002<prob<0.999998
	   RATNEST FORTRAN by
	   Best DJ & Roberts DE (1975) The percentage points of the
	   Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
	   Converted into C by Ziheng Yang, Oct. 1993.
	   */
	double e=.5e-6, aa=.6931471805, p=prob, g;
	double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;

	if (p<.000002 || p>.999998 || v<=0) return (-1);

	g = LnGamma (v/2);
	xx=v/2;   c=xx-1;
	if (v >= -1.24*log(p)) goto l1;

	ch=pow((p*xx*exp(g+xx*aa)), 1/xx);
	if (ch-e<0) return (ch);
	goto l4;
l1:
	if (v>.32) goto l3;
	ch=0.4;   a=log(1-p);
l2:
	q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
	t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
	ch-=(1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
	if (fabs(q/ch-1)-.01 <= 0) goto l4;
	else                       goto l2;

l3:
	x=PointNormal (p);
	p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0);
	if (ch>2.2*v+6)  ch=-2*(log(1-p)-c*log(.5*ch)+g);
l4:
	q=ch;   p1=.5*ch;
	if ((t=IncompleteGamma (p1, xx, g))<0) {
		printf ("\nerr IncompleteGamma");
		return (-1);
	}
	p2=p-t;
	t=p2*exp(xx*aa+g+p1-c*log(ch));
	b=t/ch;  a=0.5*t-b*c;

	s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
	s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
	s3=(210+a*(462+a*(707+932*a)))/2520;
	s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
	s5=(84+264*a+c*(175+606*a))/2520;
	s6=(120+c*(346+127*c))/5040;
	ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
	if (fabs(q/ch-1) > e) goto l4;

	return (ch);
}



#define PointGamma(prob,alpha,beta) PointChi2(prob,2.0*(alpha))/(2.0*(beta))


double CDFfunGamma(double x, double par[2])

{
	return IncompleteGamma(par[1]*x,par[0],LnGamma(par[0]));
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

#define BASE        2    /* base of floating point arithmetic */
#define DIGITS     40    /* no. of digits to the base BASE in the fraction */
#define MAXITER    30    /* max. no. of iterations to converge */

#define pos(i,j,n)      ((i)*(n)+(j))

/*A is the matrix, rr = root real(nx1), ri = root imaginary(nx1), vr = real part of eigenvector(nxn), w = working space (size of 2n), set job to 1 (decides if both eigenvectors and eigenvalues should be calculated*/

int eigen(int job, double A[], int n, double rr[], double ri[],
		double vr[], double vi[], double w[]);
void balance(double mat[], int n, int *low, int *hi, double scale[]);
void unbalance(int n, double vr[], double vi[], int low, int hi,
		double scale[]);
int realeig(int job, double mat[], int n,int low, int hi, double valr[],
		double vali[], double vr[], double vi[]);
void elemhess(int job, double mat[], int n, int low, int hi,
		double vr[], double vi[], int work[]);

int eigen(int job, double A[], int n, double rr[], double ri[],
		double vr[], double vi[], double work[])
{
	/*  double work[n*2]: working space
	*/
	int low,hi,i,j,k, it, istate=0;
	double tiny=sqrt(pow((double)BASE,(double)(1-DIGITS))), t;

	balance(A,n,&low,&hi,work);
	elemhess(job,A,n,low,hi,vr,vi, (int*)(work+n));
	if (-1 == realeig(job,A,n,low,hi,rr,ri,vr,vi)) return (-1);
	if (job) unbalance(n,vr,vi,low,hi,work);

	/* sort, added by Z. Yang */
	for (i=0; i<n; i++) {
		for (j=i+1,it=i,t=rr[i]; j<n; j++)
			if (t<rr[j]) { t=rr[j]; it=j; }
		rr[it]=rr[i];   rr[i]=t;
		t=ri[it];       ri[it]=ri[i];  ri[i]=t;
		for (k=0; k<n; k++) {
			t=vr[k*n+it];  vr[k*n+it]=vr[k*n+i];  vr[k*n+i]=t;
			t=vi[k*n+it];  vi[k*n+it]=vi[k*n+i];  vi[k*n+i]=t;
		}
		if (fabs(ri[i])>tiny) istate=1;
	}

	return (istate) ;
}

/* complex funcctions
*/

complex compl (double re,double im)
{
	complex r;

	r.re = re;
	r.im = im;
	return(r);
}

/*complex conj (complex a)
  {
  a.im = -a.im;
  */

#define csize(a) (fabs(a.re)+fabs(a.im))

complex cplus (complex a, complex b)
{
	complex c;
	c.re = a.re+b.re;
	c.im = a.im+b.im;
	return (c);
}

complex cminus (complex a, complex b)
{
	complex c;
	c.re = a.re-b.re;
	c.im = a.im-b.im;
	return (c);
}

complex cby (complex a, complex b)
{
	complex c;
	c.re = a.re*b.re-a.im*b.im ;
	c.im = a.re*b.im+a.im*b.re ;
	return (c);
}

complex cdiv (complex a,complex b)
{
	double ratio, den;
	complex c;

	if (fabs(b.re) <= fabs(b.im)) {
		ratio = b.re / b.im;
		den = b.im * (1 + ratio * ratio);
		c.re = (a.re * ratio + a.im) / den;
		c.im = (a.im * ratio - a.re) / den;
	}
	else {
		ratio = b.im / b.re;
		den = b.re * (1 + ratio * ratio);
		c.re = (a.re + a.im * ratio) / den;
		c.im = (a.im - a.re * ratio) / den;
	}
	return(c);
}

/*complex cexp (complex a)
  {
  complex c;
  c.re = exp(a.re);
  if (fabs(a.im)==0) c.im = 0;
  else  { c.im = c.re*sin(a.im); c.re*=cos(a.im); }
  return (c);
  }*/

complex cfactor (complex x, double a)
{
	complex c;
	c.re = a*x.re;
	c.im = a*x.im;
	return (c);
}

int cxtoy (complex x[], complex y[], int n)
{
	int i;
	FOR (i,n) y[i]=x[i];
	return (0);
}

int cmatby (complex a[], complex b[], complex c[], int n,int m,int k)
	/* a[n*m], b[m*k], c[n*k]  ......  c = a*b
	*/
{
	int i,j,i1;
	complex t;

	FOR (i,n)  FOR(j,k) {
		for (i1=0,t=compl(0,0); i1<m; i1++)
			t = cplus (t, cby(a[i*m+i1],b[i1*k+j]));
		c[i*k+j] = t;
	}
	return (0);
}

int cmatout (FILE * fout, complex x[], int n, int m)
{
	int i,j;
	for (i=0,FPN(fout); i<n; i++,FPN(fout))
		FOR(j,m) fprintf(fout, "%7.3f%7.3f  ", x[i*m+j].re, x[i*m+j].im);
	return (0);
}

int cmatinv( complex x[], int n, int m, double space[])
{
	/* x[n*m]  ... m>=n
	*/
	int i,j,k, *irow=(int*) space;
	double xmaxsize, ee=1e-20;
	complex xmax, t,t1;

	FOR(i,n)  {
		xmaxsize = 0.;
		for (j=i; j<n; j++) {
			if ( xmaxsize < csize (x[j*m+i]))  {
				xmaxsize = csize (x[j*m+i]);
				xmax = x[j*m+i];
				irow[i] = j;
			}
		}
		if (xmaxsize < ee)   {
			printf("\nDet goes to zero at %8d!\t\n", i+1);
			return(-1);
		}
		if (irow[i] != i) {
			FOR(j,m) {
				t = x[i*m+j];
				x[i*m+j] = x[irow[i]*m+j];
				x[ irow[i]*m+j] = t;
			}
		}
		t = cdiv (compl(1,0), x[i*m+i]);
		FOR(j,n) {
			if (j == i) continue;
			t1 = cby (t,x[j*m+i]);
			FOR(k,m)  x[j*m+k] = cminus (x[j*m+k], cby(t1,x[i*m+k]));
			x[j*m+i] = cfactor (t1, -1);
		}
		FOR(j,m)   x[i*m+j] = cby (x[i*m+j], t);
		x[i*m+i] = t;
	}
	for (i=n-1; i>=0; i--) {
		if (irow[i] == i) continue;
		FOR(j,n)  {
			t = x[j*m+i];
			x[j*m+i] = x[j*m+irow[i]];
			x[ j*m+irow[i]] = t;
		}
	}
	return (0);
}


void balance(double mat[], int n,int *low, int *hi, double scale[])
{
	/* Balance a matrix for calculation of eigenvalues and eigenvectors
	*/
	double c,f,g,r,s;
	int i,j,k,l,done;
	/* search for rows isolating an eigenvalue and push them down */
	for (k = n - 1; k >= 0; k--) {
		for (j = k; j >= 0; j--) {
			for (i = 0; i <= k; i++) {
				if (i != j && fabs(mat[pos(j,i,n)]) != 0) break;
			}

			if (i > k) {
				scale[k] = j;

				if (j != k) {
					for (i = 0; i <= k; i++) {
						c = mat[pos(i,j,n)];
						mat[pos(i,j,n)] = mat[pos(i,k,n)];
						mat[pos(i,k,n)] = c;
					}

					for (i = 0; i < n; i++) {
						c = mat[pos(j,i,n)];
						mat[pos(j,i,n)] = mat[pos(k,i,n)];
						mat[pos(k,i,n)] = c;
					}
				}
				break;
			}
		}
		if (j < 0) break;
	}

	/* search for columns isolating an eigenvalue and push them left */

	for (l = 0; l <= k; l++) {
		for (j = l; j <= k; j++) {
			for (i = l; i <= k; i++) {
				if (i != j && fabs(mat[pos(i,j,n)]) != 0) break;
			}
			if (i > k) {
				scale[l] = j;
				if (j != l) {
					for (i = 0; i <= k; i++) {
						c = mat[pos(i,j,n)];
						mat[pos(i,j,n)] = mat[pos(i,l,n)];
						mat[pos(i,l,n)] = c;
					}

					for (i = l; i < n; i++) {
						c = mat[pos(j,i,n)];
						mat[pos(j,i,n)] = mat[pos(l,i,n)];
						mat[pos(l,i,n)] = c;
					}
				}

				break;
			}
		}

		if (j > k) break;
	}

	*hi = k;
	*low = l;

	/* balance the submatrix in rows l through k */

	for (i = l; i <= k; i++) {
		scale[i] = 1;
	}

	do {
		for (done = 1,i = l; i <= k; i++) {
			for (c = 0,r = 0,j = l; j <= k; j++) {
				if (j != i) {
					c += fabs(mat[pos(j,i,n)]);
					r += fabs(mat[pos(i,j,n)]);
				}
			}

			if (c != 0 && r != 0) {
				g = r / BASE;
				f = 1;
				s = c + r;

				while (c < g) {
					f *= BASE;
					c *= BASE * BASE;
				}

				g = r * BASE;

				while (c >= g) {
					f /= BASE;
					c /= BASE * BASE;
				}

				if ((c + r) / f < 0.95 * s) {
					done = 0;
					g = 1 / f;
					scale[i] *= f;

					for (j = l; j < n; j++) {
						mat[pos(i,j,n)] *= g;
					}

					for (j = 0; j <= k; j++) {
						mat[pos(j,i,n)] *= f;
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
void unbalance(int n,double vr[],double vi[], int low, int hi, double scale[])
{
	int i,j,k;
	double tmp;

	for (i = low; i <= hi; i++) {
		for (j = 0; j < n; j++) {
			vr[pos(i,j,n)] *= scale[i];
			vi[pos(i,j,n)] *= scale[i];
		}
	}

	for (i = low - 1; i >= 0; i--) {
		if ((k = (int)scale[i]) != i) {
			for (j = 0; j < n; j++) {
				tmp = vr[pos(i,j,n)];
				vr[pos(i,j,n)] = vr[pos(k,j,n)];
				vr[pos(k,j,n)] = tmp;

				tmp = vi[pos(i,j,n)];
				vi[pos(i,j,n)] = vi[pos(k,j,n)];
				vi[pos(k,j,n)] = tmp;
			}
		}
	}

	for (i = hi + 1; i < n; i++) {
		if ((k = (int)scale[i]) != i) {
			for (j = 0; j < n; j++) {
				tmp = vr[pos(i,j,n)];
				vr[pos(i,j,n)] = vr[pos(k,j,n)];
				vr[pos(k,j,n)] = tmp;

				tmp = vi[pos(i,j,n)];
				vi[pos(i,j,n)] = vi[pos(k,j,n)];
				vi[pos(k,j,n)] = tmp;
			}
		}
	}
}

/*
 * Reduce the submatrix in rows and columns low through hi of real matrix mat to
 * Hessenberg form by elementary similarity transformations
 */
void elemhess(int job,double mat[],int n,int low,int hi, double vr[],
		double vi[], int work[])
{
	/* work[n] */
	int i,j,m;
	double x,y;

	for (m = low + 1; m < hi; m++) {
		for (x = 0,i = m,j = m; j <= hi; j++) {
			if (fabs(mat[pos(j,m-1,n)]) > fabs(x)) {
				x = mat[pos(j,m-1,n)];
				i = j;
			}
		}

		if ((work[m] = i) != m) {
			for (j = m - 1; j < n; j++) {
				y = mat[pos(i,j,n)];
				mat[pos(i,j,n)] = mat[pos(m,j,n)];
				mat[pos(m,j,n)] = y;
			}

			for (j = 0; j <= hi; j++) {
				y = mat[pos(j,i,n)];
				mat[pos(j,i,n)] = mat[pos(j,m,n)];
				mat[pos(j,m,n)] = y;
			}
		}

		if (x != 0) {
			for (i = m + 1; i <= hi; i++) {
				if ((y = mat[pos(i,m-1,n)]) != 0) {
					y = mat[pos(i,m-1,n)] = y / x;

					for (j = m; j < n; j++) {
						mat[pos(i,j,n)] -= y * mat[pos(m,j,n)];
					}

					for (j = 0; j <= hi; j++) {
						mat[pos(j,m,n)] += y * mat[pos(j,i,n)];
					}
				}
			}
		}
	}
	if (job) {
		for (i=0; i<n; i++) {
			for (j=0; j<n; j++) {
				vr[pos(i,j,n)] = 0.0; vi[pos(i,j,n)] = 0.0;
			}
			vr[pos(i,i,n)] = 1.0;
		}

		for (m = hi - 1; m > low; m--) {
			for (i = m + 1; i <= hi; i++) {
				vr[pos(i,m,n)] = mat[pos(i,m-1,n)];
			}

			if ((i = work[m]) != m) {
				for (j = m; j <= hi; j++) {
					vr[pos(m,j,n)] = vr[pos(i,j,n)];
					vr[pos(i,j,n)] = 0.0;
				}
				vr[pos(i,m,n)] = 1.0;
			}
		}
	}
}

/*
 * Calculate eigenvalues and eigenvectors of a real upper Hessenberg matrix
 * Return 1 if converges successfully and 0 otherwise
 */

int realeig(int job,double mat[],int n,int low, int hi, double valr[],
		double vali[], double vr[],double vi[])
{
	complex v;
	double p=0,q=0,r=0,s=0,t,w,x,y,z=0,ra,sa,norm,eps;
	int niter,en,i,j,k,l,m;
	double precision  = pow((double)BASE,(double)(1-DIGITS));

	eps = precision;
	for (i=0; i<n; i++) {
		valr[i]=0.0;
		vali[i]=0.0;
	}
	/* store isolated roots and calculate norm */
	for (norm = 0,i = 0; i < n; i++) {
		for (j = max(0,i-1); j < n; j++) {
			norm += fabs(mat[pos(i,j,n)]);
		}
		if (i < low || i > hi) valr[i] = mat[pos(i,i,n)];
	}
	t = 0;
	en = hi;

	while (en >= low) {
		niter = 0;
		for (;;) {

			/* look for single small subdiagonal element */

			for (l = en; l > low; l--) {
				s = fabs(mat[pos(l-1,l-1,n)]) + fabs(mat[pos(l,l,n)]);
				if (s == 0) s = norm;
				if (fabs(mat[pos(l,l-1,n)]) <= eps * s) break;
			}

			/* form shift */

			x = mat[pos(en,en,n)];

			if (l == en) {             /* one root found */
				valr[en] = x + t;
				if (job) mat[pos(en,en,n)] = x + t;
				en--;
				break;
			}

			y = mat[pos(en-1,en-1,n)];
			w = mat[pos(en,en-1,n)] * mat[pos(en-1,en,n)];

			if (l == en - 1) {                /* two roots found */
				p = (y - x) / 2;
				q = p * p + w;
				z = sqrt(fabs(q));
				x += t;
				if (job) {
					mat[pos(en,en,n)] = x;
					mat[pos(en-1,en-1,n)] = y + t;
				}
				if (q < 0) {                /* complex pair */
					valr[en-1] = x+p;
					vali[en-1] = z;
					valr[en] = x+p;
					vali[en] = -z;
				}
				else {                      /* real pair */
					z = (p < 0) ? p - z : p + z;
					valr[en-1] = x + z;
					valr[en] = (z == 0) ? x + z : x - w / z;
					if (job) {
						x = mat[pos(en,en-1,n)];
						s = fabs(x) + fabs(z);
						p = x / s;
						q = z / s;
						r = sqrt(p*p+q*q);
						p /= r;
						q /= r;
						for (j = en - 1; j < n; j++) {
							z = mat[pos(en-1,j,n)];
							mat[pos(en-1,j,n)] = q * z + p *
								mat[pos(en,j,n)];
							mat[pos(en,j,n)] = q * mat[pos(en,j,n)] - p*z;
						}
						for (i = 0; i <= en; i++) {
							z = mat[pos(i,en-1,n)];
							mat[pos(i,en-1,n)] = q * z + p * mat[pos(i,en,n)];
							mat[pos(i,en,n)] = q * mat[pos(i,en,n)] - p*z;
						}
						for (i = low; i <= hi; i++) {
							z = vr[pos(i,en-1,n)];
							vr[pos(i,en-1,n)] = q*z + p*vr[pos(i,en,n)];
							vr[pos(i,en,n)] = q*vr[pos(i,en,n)] - p*z;
						}
					}
				}
				en -= 2;
				break;
			}
			if (niter == MAXITER) return(-1);
			if (niter != 0 && niter % 10 == 0) {
				t += x;
				for (i = low; i <= en; i++) mat[pos(i,i,n)] -= x;
				s = fabs(mat[pos(en,en-1,n)]) + fabs(mat[pos(en-1,en-2,n)]);
				x = y = 0.75 * s;
				w = -0.4375 * s * s;
			}
			niter++;
			/* look for two consecutive small subdiagonal elements */
			for (m = en - 2; m >= l; m--) {
				z = mat[pos(m,m,n)];
				r = x - z;
				s = y - z;
				p = (r * s - w) / mat[pos(m+1,m,n)] + mat[pos(m,m+1,n)];
				q = mat[pos(m+1,m+1,n)] - z - r - s;
				r = mat[pos(m+2,m+1,n)];
				s = fabs(p) + fabs(q) + fabs(r);
				p /= s;
				q /= s;
				r /= s;
				if (m == l || fabs(mat[pos(m,m-1,n)]) * (fabs(q)+fabs(r)) <=
						eps * (fabs(mat[pos(m-1,m-1,n)]) + fabs(z) +
							fabs(mat[pos(m+1,m+1,n)])) * fabs(p)) break;
			}
			for (i = m + 2; i <= en; i++) mat[pos(i,i-2,n)] = 0;
			for (i = m + 3; i <= en; i++) mat[pos(i,i-3,n)] = 0;
			/* double QR step involving rows l to en and columns m to en */
			for (k = m; k < en; k++) {
				if (k != m) {
					p = mat[pos(k,k-1,n)];
					q = mat[pos(k+1,k-1,n)];
					r = (k == en - 1) ? 0 : mat[pos(k+2,k-1,n)];
					if ((x = fabs(p) + fabs(q) + fabs(r)) == 0) continue;
					p /= x;
					q /= x;
					r /= x;
				}
				s = sqrt(p*p+q*q+r*r);
				if (p < 0) s = -s;
				if (k != m) {
					mat[pos(k,k-1,n)] = -s * x;
				}
				else if (l != m) {
					mat[pos(k,k-1,n)] = -mat[pos(k,k-1,n)];
				}
				p += s;
				x = p / s;
				y = q / s;
				z = r / s;
				q /= p;
				r /= p;
				/* row modification */
				for (j = k; j <= (!job ? en : n-1); j++){
					p = mat[pos(k,j,n)] + q * mat[pos(k+1,j,n)];
					if (k != en - 1) {
						p += r * mat[pos(k+2,j,n)];
						mat[pos(k+2,j,n)] -= p * z;
					}
					mat[pos(k+1,j,n)] -= p * y;
					mat[pos(k,j,n)] -= p * x;
				}
				j = min(en,k+3);
				/* column modification */
				for (i = (!job ? l : 0); i <= j; i++) {
					p = x * mat[pos(i,k,n)] + y * mat[pos(i,k+1,n)];
					if (k != en - 1) {
						p += z * mat[pos(i,k+2,n)];
						mat[pos(i,k+2,n)] -= p*r;
					}
					mat[pos(i,k+1,n)] -= p*q;
					mat[pos(i,k,n)] -= p;
				}
				if (job) {             /* accumulate transformations */
					for (i = low; i <= hi; i++) {
						p = x * vr[pos(i,k,n)] + y * vr[pos(i,k+1,n)];
						if (k != en - 1) {
							p += z * vr[pos(i,k+2,n)];
							vr[pos(i,k+2,n)] -= p*r;
						}
						vr[pos(i,k+1,n)] -= p*q;
						vr[pos(i,k,n)] -= p;
					}
				}
			}
		}
	}

	if (!job) return(0);
	if (norm != 0) {
		/* back substitute to find vectors of upper triangular form */
		for (en = n-1; en >= 0; en--) {
			p = valr[en];
			if ((q = vali[en]) < 0) {            /* complex vector */
				m = en - 1;
				if (fabs(mat[pos(en,en-1,n)]) > fabs(mat[pos(en-1,en,n)])) {
					mat[pos(en-1,en-1,n)] = q / mat[pos(en,en-1,n)];
					mat[pos(en-1,en,n)] = (p - mat[pos(en,en,n)]) /
						mat[pos(en,en-1,n)];
				}
				else {
					v = cdiv(compl(0.0,-mat[pos(en-1,en,n)]),
							compl(mat[pos(en-1,en-1,n)]-p,q));
					mat[pos(en-1,en-1,n)] = v.re;
					mat[pos(en-1,en,n)] = v.im;
				}
				mat[pos(en,en-1,n)] = 0;
				mat[pos(en,en,n)] = 1;
				for (i = en - 2; i >= 0; i--) {
					w = mat[pos(i,i,n)] - p;
					ra = 0;
					sa = mat[pos(i,en,n)];
					for (j = m; j < en; j++) {
						ra += mat[pos(i,j,n)] * mat[pos(j,en-1,n)];
						sa += mat[pos(i,j,n)] * mat[pos(j,en,n)];
					}
					if (vali[i] < 0) {
						z = w;
						r = ra;
						s = sa;
					}
					else {
						m = i;
						if (vali[i] == 0) {
							v = cdiv(compl(-ra,-sa),compl(w,q));
							mat[pos(i,en-1,n)] = v.re;
							mat[pos(i,en,n)] = v.im;
						}
						else {                      /* solve complex equations */
							x = mat[pos(i,i+1,n)];
							y = mat[pos(i+1,i,n)];
							v.re = (valr[i]- p)*(valr[i]-p) + vali[i]*vali[i] - q*q;
							v.im = (valr[i] - p)*2*q;
							if ((fabs(v.re) + fabs(v.im)) == 0) {
								v.re = eps * norm * (fabs(w) +
										fabs(q) + fabs(x) + fabs(y) + fabs(z));
							}
							v = cdiv(compl(x*r-z*ra+q*sa,x*s-z*sa-q*ra),v);
							mat[pos(i,en-1,n)] = v.re;
							mat[pos(i,en,n)] = v.im;
							if (fabs(x) > fabs(z) + fabs(q)) {
								mat[pos(i+1,en-1,n)] =
									(-ra - w * mat[pos(i,en-1,n)] +
									 q * mat[pos(i,en,n)]) / x;
								mat[pos(i+1,en,n)] = (-sa - w * mat[pos(i,en,n)] -
										q * mat[pos(i,en-1,n)]) / x;
							}
							else {
								v = cdiv(compl(-r-y*mat[pos(i,en-1,n)],
											-s-y*mat[pos(i,en,n)]),compl(z,q));
								mat[pos(i+1,en-1,n)] = v.re;
								mat[pos(i+1,en,n)] = v.im;
							}
						}
					}
				}
			}
			else if (q == 0) {                             /* real vector */
				m = en;
				mat[pos(en,en,n)] = 1;
				for (i = en - 1; i >= 0; i--) {
					w = mat[pos(i,i,n)] - p;
					r = mat[pos(i,en,n)];
					for (j = m; j < en; j++) {
						r += mat[pos(i,j,n)] * mat[pos(j,en,n)];
					}
					if (vali[i] < 0) {
						z = w;
						s = r;
					}
					else {
						m = i;
						if (vali[i] == 0) {
							if ((t = w) == 0) t = eps * norm;
							mat[pos(i,en,n)] = -r / t;
						}
						else {            /* solve real equations */
							x = mat[pos(i,i+1,n)];
							y = mat[pos(i+1,i,n)];
							q = (valr[i] - p) * (valr[i] - p) + vali[i]*vali[i];
							t = (x * s - z * r) / q;
							mat[pos(i,en,n)] = t;
							if (fabs(x) <= fabs(z)) {
								mat[pos(i+1,en,n)] = (-s - y * t) / z;
							}
							else {
								mat[pos(i+1,en,n)] = (-r - w * t) / x;
							}
						}
					}
				}
			}
		}
		/* vectors of isolated roots */
		for (i = 0; i < n; i++) {
			if (i < low || i > hi) {
				for (j = i; j < n; j++) {
					vr[pos(i,j,n)] = mat[pos(i,j,n)];
				}
			}
		}
		/* multiply by transformation matrix */

		for (j = n-1; j >= low; j--) {
			m = min(j,hi);
			for (i = low; i <= hi; i++) {
				for (z = 0,k = low; k <= m; k++) {
					z += vr[pos(i,k,n)] * mat[pos(k,j,n)];
				}
				vr[pos(i,j,n)] = z;
			}
		}
	}
	/* rearrange complex eigenvectors */
	for (j = 0; j < n; j++) {
		if (vali[j] != 0) {
			for (i = 0; i < n; i++) {
				vi[pos(i,j,n)] = vr[pos(i,j+1,n)];
				vr[pos(i,j+1,n)] = vr[pos(i,j,n)];
				vi[pos(i,j+1,n)] = -vi[pos(i,j,n)];
			}
			j++;
		}
	}
	return(0);
}

int matinv( double x[], int n, int m, double space[])
{
	/* x[n*m]  ... m>=n
	*/
	register int i,j,k;
	int *irow=(int*) space;
	double ee=1.0e-20, t,t1,xmax;
	double det=1.0;

	FOR (i,n)  {
		xmax = 0.;
		for (j=i; j<n; j++) {
			if (xmax < fabs(x[j*m+i]))  {
				xmax = fabs( x[j*m+i] );
				irow[i] = j;
			}
		}
		det *= xmax;
		if (xmax < ee)   {
			printf("\nDet becomes zero at %3d!\t\n", i+1);
			return(-1);
		}
		if (irow[i] != i) {
			FOR (j,m) {
				t = x[i*m+j];
				x[i*m+j] = x[irow[i] * m + j];
				x[ irow[i] * m + j] = t;
			}
		}
		t = 1./x[i*m+i];
		FOR (j,n) {
			if (j == i) continue;
			t1 = t*x[j*m+i];
			FOR(k,m)  x[j*m+k] -= t1*x[i*m+k];
			x[j*m+i] = -t1;
		}
		FOR(j,m)   x[i*m+j] *= t;
		x[i*m+i] = t;
	}                            /* i  */
	for (i=n-1; i>=0; i--) {
		if (irow[i] == i) continue;
		FOR(j,n)  {
			t = x[j*m+i];
			x[j*m+i] = x[ j*m + irow[i] ];
			x[ j*m + irow[i] ] = t;
		}
	}
	return (0);
}

void make_transition_prob_matrices(double t[3], int treeNum)
{
	int i, j, k, v, n;
	double EXPOS[4], T;

	for (v=0; v<3; v++)
	{
		//printf("branch %i\n",v);
		for (n=0; n<NUMCAT; n++)
		{
			//printf("Category %i\n",n);
			T = statevector[treeNum][n]*t[v];
			for (k=0; k<4; k++)
				EXPOS[k] = exp(T*RRVAL[treeNum][k]);
			for (i=0; i<4; i++)
			{
				for (j=0; j<4; j++)
				{
					PMAT[v][n][i][j] = 0.0;
					for (k=0; k<4; k++)
						PMAT[v][n][i][j] += RRVEC[treeNum][k * 4 + j]*LRVEC[treeNum][i * 4 + k]*EXPOS[k];
					//printf("%.16f:",PMAT[v][n][i][j]);
					if(PMAT[v][n][i][j] <= 0.0)
					{
						PMAT[v][n][i][j] = 0.00000001;
						// printf("%.16f:",PMAT[v][n][i][j]);
					}
					PMAT[v][n][i][j] = log(PMAT[v][n][i][j]);
					//printf("%.16f\t",PMAT[v][n][i][j]);
				} //printf("\n");
			}
		}
	}
}

double mydistance(double v1[], double v2[], int start1, int start2, int n)
{
	int i;
	double sum=0.0;

	for (i=0; i<n; i++)
		sum+=(v1[i+start1]-v2[i+start2])*(v1[i+start1]-v2[i+start2]);
	return sqrt(sum);
}

// log sum exp trick sum(exp(X)-maxX) + maxX
//Calculates the log sum exp of log like + log prob for the tree calculations
//This was implemented after scaling at base/pos was underflowing
// always used on dim 4 vectors
double logSumExp(double X[4])
{
	int i;
	double maxX = X[0];
	double sumX = 0.0;

	//Find max
	for(i = 1; i < 4; i++)
	{
		if(X[i] > maxX)
			maxX = X[i];
	}

	//sum(exp(x) - maxX)
	for(i = 0; i < 4; i++)
	{
		//printf("\t\t\t\t\tX[%d]: %.16f\n", i, X[i]);
		sumX += exp(X[i] - maxX);
	}

	return(log(sumX) + maxX);
}


//transtion probability matrix has already been diagonalized
//this function assumes that the only thing that changes between function calls is the branch lengths or the query sequence and its placement
//node1 and node2 are the nodes around the edge to which the sequence has been assigned
//seq is the identifer of the sequence
//pi are the nucletoide frequencies - this could be made a global to avoid passing them around
//Alignment data should be stored in DATA with -1 indicating missing data
//upper and lower bounds are parameters 7-10 when 0-counting
double getlike_gamma_root_in_trifurcation(double times[3], double parameters[7])
{
	int i, j,k, b, c, v, po, node, seq, treeNum; 
	double Like, t[3], A[4], B[4], C[4];
	//diagonalizaiton has previously been done   

	if (onDindic==1) {	// Learn more about what these do, seem to be hitting these values
		if (times[2] < parameters[3] || times[2] > parameters[4]){
			//printf("\t\tError 1: %.16f < %.16f or %.16f > %.16f\n", times[2], parameters[3], times[2], parameters[4]);
			return 1000000000.0;//If an actual likelihood is smaler than this we are screwed
		}
	}
	else {
		if ((times[2] < parameters[5] || times[2] > parameters[6]) || (times[1] < parameters[3] || times[1] > parameters[4])){
			//printf("\t\tError 2: (%.16f < %.16f or %.16f > %.16f) or (%.16f < %.16f or %.16f > %.16f)\n", times[2], parameters[5], times[2], parameters[6], times[1], parameters[3], times[1], parameters[4]);
			return 1000000000.0;//If an actual likelihood is smaler than this we are screwed
		}
	}	

	node = (int)parameters[2];
	seq = (int)parameters[0];
	treeNum = (int)parameters[1];
	//	printf("Analysing read %d\n", seq);
	//Why copy this? TO DO: Remove this copy and change below code to reflect the change
	//for (i=0; i<4; i++) 
	//	printf("pi[%d]: %.16lf\t", i, pi[treeNum][i]);
	//	pi[i]=parameters[i+3];
	
	//printf("\n");

	t[1] = times[2];  //length from node to position where query joins
	t[0] = (times[2]+nodeages[treeNum][node])*times[1]; //length from age of query node to position where query joins
	t[2] = bls[treeNum][node]-times[2]; //length from position where query joins to parent node
	
	// printf("node: %d\tseq: %d\tt0: %.16f\tt1: %.16f\tt2: %.16f\n", node, seq, t[0], t[1], t[2]);
	//printf("t0: %.16f\tt1: %.16f\tt2: %.16f\n", t[0], t[1], t[2]);

	make_transition_prob_matrices(t, treeNum);

	//printf("Errors setting of %d, sequence %d, assignment %d\n", errors, seq, node);

	Like=0.0;

	// printf("Sequence %d\n", seq);

	for (i=startpos[seq]; i<readlength[seq]+startpos[seq]; i++){ 
		po = i-startpos[seq];//i keeps track of the positon in the ref sequences while po is the positoon in the read
		b = QUERYDATA[seq][po];

		// printf("\tBase %d\n", b);

		if (b!=-1){
			for (j=0; j<NUMCAT; j++){
				// printf("\t\tj: %d\n", j);
				for (k=0; k<4; k++){
					//printf("\t\t\tk: %d\n", k);
					for (v=0; v<4; v++)
					{
						A[v] = pi[treeNum][v] + PMAT[0][j][v][k] + readlike[seq][po][v];
					//	//printf("\t\t\t\tpi[%d] is %.16f anf PMAT[0] is %.16f anf reaflike is %.16f\n", v, pi[treeNum][v], PMAT[0][j][v][k], readlike[seq][po][v]);
					}
					B[k] = logSumExp(A);
					//B[k] = pi[treeNum][b] + PMAT[0][j][b][k];
					// printf("\t\t\t\t\tB:%lf\tpi:%lf\tPMAT[0]:%lf\n", B[k], pi[treeNum][b], PMAT[0][j][b][k]);
					if (node>=numseq[treeNum]){//If not leaf node. Assumes the leaf nodes are numbered from 0 to numseq-1
						for (v=0; v<4; v++)
						{
							// printf("\t\t\t\t\tPMAT[1] is %.16f and FRACLIKE[%d] is %.16f\n", PMAT[1][j][k][v], v+4, FRACLIKE[treeNum][i][node][j][v+4]);
							A[v] = PMAT[1][j][k][v] + FRACLIKE[treeNum][i][node][j][v+4];
						}
						B[k] += logSumExp(A);
						// printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
						for (v=0; v<4; v++)
						{
							// printf("\t\t\t\t\tPMAT[2] is %.16f and FRACLIKE[%d] is %.16f\n", PMAT[2][j][k][v], v, FRACLIKE[treeNum][i][node][j][v]);
							A[v] = PMAT[2][j][k][v] + FRACLIKE[treeNum][i][node][j][v];
						}
						B[k] += logSumExp(A);
						// printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
					}
					else { //If leaf node
						if ((c=DATA[treeNum][node][i])>-1)
						{
							B[k] += PMAT[1][j][k][c];
						}
						for (v=0; v<4; v++)//add position here into fraclike
						{
							// printf("\t\t\t\t\tPMAT[2] is %.16f and FRACLIKE[%d] is %.16f\n", PMAT[2][j][k][v], v, FRACLIKE[treeNum][i][node][j][v]);
							A[v] = PMAT[2][j][k][v] + FRACLIKE[treeNum][i][node][j][v];
						}
						B[k] += logSumExp(A);
						// printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
					}
					// printf("\t\t\t\t\tB[%d]: %lf\n", k, B[k]);
				}
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

	//printf("\tLikelihood: %.16f\n", Like);
	//printf("\t\t%d,%f,%f,%f,%f\n", seq, (1.0-times[1])*(nodeages[treeNum][seq]+times[2]), times[1], times[2], Like);

	return -Like; //Notice: a scaling factor of NUMCAT^(number of sites) is missing
}

double getlike_gamma_root_in_trifurcation_Print_Lik(double times[3], double parameters[7])
{
	int i, j,k, b, c, v, po, node, seq, treeNum; 
	double Like, t[3], A[4], B[4], C[4];
	//diagonalizaiton has previously been done   

	if (onDindic==1) {	// Learn more about what these do, seem to be hitting these values
		if (times[2] < parameters[3] || times[2] > parameters[4]){
			//printf("\t\tError 1: %.16f < %.16f or %.16f > %.16f\n", times[2], parameters[3], times[2], parameters[4]);
			return 1000000000.0;//If an actual likelihood is smaler than this we are screwed
		}
	}
	else {
		if ((times[2] < parameters[5] || times[2] > parameters[6]) || (times[1] < parameters[3] || times[1] > parameters[4])){
			//printf("\t\tError 2: (%.16f < %.16f or %.16f > %.16f) or (%.16f < %.16f or %.16f > %.16f)\n", times[2], parameters[5], times[2], parameters[6], times[1], parameters[3], times[1], parameters[4]);
			return 1000000000.0;//If an actual likelihood is smaler than this we are screwed
		}
	}	

	node = (int)parameters[2];
	seq = (int)parameters[0];
	treeNum = (int)parameters[1];
	//	printf("Analysing read %d\n", seq);
	//Why copy this? TO DO: Remove this copy and change below code to reflect the change
	//for (i=0; i<4; i++) 
	//	printf("pi[%d]: %.16lf\t", i, pi[treeNum][i]);
	//	pi[i]=parameters[i+3];
	
	//printf("\n");

	t[1] = times[2];  //length from node to position where query joins
	t[0] = (times[2]+nodeages[treeNum][node])*times[1]; //length from age of query node to position where query joins
	t[2] = bls[treeNum][node]-times[2]; //length from position where query joins to parent node
	
	//printf("node: %d\tseq: %d\ttree: %d\tt0: %.16f\tt1: %.16f\tt2: %.16f\n", node, seq, treeNum, t[0], t[1], t[2]);
	//printf("t0: %.16f\tt1: %.16f\tt2: %.16f\n", t[0], t[1], t[2]);

	make_transition_prob_matrices(t, treeNum);

	//printf("Errors setting of %d, sequence %d, assignment %d\n", errors, seq, node);

	Like=0.0;

	// printf("Sequence %d\n", seq);

	for (i=startpos[seq]; i<readlength[seq]+startpos[seq]; i++){ 
		po = i-startpos[seq];//i keeps track of the positon in the ref sequences while po is the positoon in the read
		b = QUERYDATA[seq][po];

		//printf("\tBase %d of pos %d with start %d and readlength %d\n", b, po, startpos[seq], readlength[seq]);

		if (b!=-1){
			for (j=0; j<NUMCAT; j++){
				// printf("\t\tj: %d\n", j);
				for (k=0; k<4; k++){
					//printf("\t\t\tk: %d\n", k);
					for (v=0; v<4; v++)
					{
						A[v] = pi[treeNum][v] + PMAT[0][j][v][k] + readlike[seq][po][v];
						//printf("\t\t\t\tpi[%d] is %.16f anf PMAT[0] is %.16f anf reaflike is %.16f\n", v, pi[treeNum][v], PMAT[0][j][v][k], readlike[seq][po][v]);
					}
					B[k] = logSumExp(A);
					//B[k] = pi[treeNum][b] + PMAT[0][j][b][k];
					// printf("\t\t\t\t\tB:%lf\tpi:%lf\tPMAT[0]:%lf\n", B[k], pi[treeNum][b], PMAT[0][j][b][k]);
					if (node>=numseq[treeNum]){//If not leaf node. Assumes the leaf nodes are numbered from 0 to numseq-1
						for (v=0; v<4; v++)
						{
							// printf("\t\t\t\t\tPMAT[1] is %.16f and FRACLIKE[%d] is %.16f\n", PMAT[1][j][k][v], v+4, FRACLIKE[treeNum][i][node][j][v+4]);
							A[v] = PMAT[1][j][k][v] + FRACLIKE[treeNum][i][node][j][v+4];
						}
						B[k] += logSumExp(A);
						// printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
						for (v=0; v<4; v++)
						{
							// printf("\t\t\t\t\tPMAT[2] is %.16f and FRACLIKE[%d] is %.16f\n", PMAT[2][j][k][v], v, FRACLIKE[treeNum][i][node][j][v]);
							A[v] = PMAT[2][j][k][v] + FRACLIKE[treeNum][i][node][j][v];
						}
						B[k] += logSumExp(A);
						// printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
					}
					else { //If leaf node
						if ((c=DATA[treeNum][node][i])>-1)
						{
							B[k] += PMAT[1][j][k][c];
						}
						for (v=0; v<4; v++)//add position here into fraclike
						{
							// printf("\t\t\t\t\tPMAT[2] is %.16f and FRACLIKE[%d] is %.16f\n", PMAT[2][j][k][v], v, FRACLIKE[treeNum][i][node][j][v]);
							A[v] = PMAT[2][j][k][v] + FRACLIKE[treeNum][i][node][j][v];
						}
						B[k] += logSumExp(A);
						// printf("\t\t\t\t\tB:%lf\n", logSumExp(A));
					}
					// printf("\t\t\t\t\tB[%d]: %lf\n", k, B[k]);
				}
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

	//printf("\tLikelihood: %.16f\n", Like);
	printf("\t\t%d,%d,%d,%f,%f,%f,%f\n", seq, node, treeNum, (1.0-times[1])*(nodeages[treeNum][seq]+times[2]), times[1], times[2], Like);

	return -Like; //Notice: a scaling factor of NUMCAT^(number of sites) is missing
}

// Same as above, but used to print out site scores
double getlike_gamma_root_in_trifurcation_Print(double times[3], double parameters[3])
{
	int i, j,k, b, c, v, po, node, seq, treeNum; 
	double Like, t[3], A[4], B[4], C[4];

	//Don't need to do the error checking like normally

	node = (int)parameters[2];
	seq = (int)parameters[0];
	treeNum = (int)parameters[1];
	//for (i=0; i<4; i++) 
	//	pi[i]=parameters[i+3];

	t[1] = times[2];  //length from node to position where query joins
	t[0] = (times[2]+nodeages[treeNum][node])*times[1]; //length from age of query node to position where query joins
	t[2] = bls[treeNum][node]-times[2]; //length from position where query joins to parent node
	
	//printf("node: %d\tseq: %d\ttree: %d\tt0: %.16f\tt1: %.16f\tt2: %.16f\n", node, seq, treeNum, t[0], t[1], t[2]);

	make_transition_prob_matrices(t, treeNum);


	Like=0.0;


	for (i=startpos[seq]; i<readlength[seq]+startpos[seq]; i++){ 
		po = i-startpos[seq];
		b = QUERYDATA[seq][po];


		if (b!=-1){
			for (j=0; j<NUMCAT; j++){
				for (k=0; k<4; k++){
					for (v=0; v<4; v++)
					{
						A[v] = pi[treeNum][v] + PMAT[0][j][v][k] + readlike[seq][po][v];
					}
					B[k] = logSumExp(A);
					//B[k] = pi[treeNum][b] + PMAT[0][j][b][k];
					if (node>=numseq[treeNum]){//If not leaf node. Assumes the leaf nodes are numbered from 0 to numseq-1
						for (v=0; v<4; v++)
						{
							A[v] = PMAT[1][j][k][v] + FRACLIKE[treeNum][i][node][j][v+4];
						}
						B[k] += logSumExp(A);
						for (v=0; v<4; v++)
						{
							A[v] = PMAT[2][j][k][v] + FRACLIKE[treeNum][i][node][j][v];
						}
						B[k] += logSumExp(A);
					}
					else { //If leaf node
						if ((c=DATA[treeNum][node][i])>-1)
						{
							B[k] += PMAT[1][j][k][c];
						}
						for (v=0; v<4; v++)//add position here into fraclike
						{
							A[v] = PMAT[2][j][k][v] + FRACLIKE[treeNum][i][node][j][v];
						}
						B[k] += logSumExp(A);
					}
				}
				C[j] = logSumExp(B);

			}
			printf("%d, %.16f\n", i, logSumExp(C));
			Like += logSumExp(C);
		}
	}

	return -Like; //Notice: a scaling factor of NUMCAT^(number of sites) is missing
}

double getlike_gamma_root_in_trifurcation_L0(double times[3], double parameters[7])
{

	double T[3];

	T[2] = times[1];  
	T[1] = 1.0;
	return getlike_gamma_root_in_trifurcation(T, parameters);
}

double getlike_gamma_root_in_trifurcation_testAge(double rootPlace, double parameters[7])
{
	double T[3], returnValue;

	T[2] = rootPlace;
	T[1] = 1.0 - testAge/(nodeages[(int)parameters[1]][(int)parameters[2]] + T[2]);

	// Since its been switched to 2D optimization, allow for the bounds to catch <0 and >1
	if (T[1] < 0.0)
	{
		// This means that the testAge > node + T[2], requiring a "negative" branch. So changing to 0 does not make sense and should not be done
		//printf("\t\ttestAge of %.16f: T[1] < 0.0 of %.16f!!\n", testAge, T[1]); 
		return 1000000000.0;
	}
	else if (T[1] > 1.0)	//This must be between 0 and 1
	{
		// This only would happen if negative T[2] + node age occurs, which shouldn't? Again, makes more sense to penalize heavily with high likelihood instead of changing to 1
		//printf("\t\ttestAge of %.16f: T[1] > 1.0 of %.16f!!\n", testAge, T[1]);
		return 1000000000.0;
	}

	return getlike_gamma_root_in_trifurcation(T, parameters);
}


void inittransitionmatrix()

{
	int i, j;
	double sum, RIVAL[4], RIVEC[4][4],  A[4][4], workspace[8];

	//TO DO: Can probably make these a single malloc
	RRVAL = (double **)malloc(numTrees * sizeof(double *));
	RRVEC = (double **)malloc(numTrees * sizeof(double *));
	LRVEC = (double **)malloc(numTrees * sizeof(double *));

	for(int treeNum = 0; treeNum < numTrees; treeNum++)
	{
		if (usedTrees[treeNum] == 0)
		{
			// no reads were assigned to this tree, can skip
			continue;
		}
		RRVAL[treeNum] = (double *)malloc(4 * sizeof(double));
		RRVEC[treeNum] = (double *)malloc(4*4*sizeof(double));
			//(double **)malloc(4 * sizeof(double *));
		LRVEC[treeNum] = (double *)malloc(4*4*sizeof(double));
			//(double **)malloc(4 * sizeof(double *));
		//for(i = 0; i < 4; i++)
		//{
		//	RRVEC[treeNum][i] = (double *)malloc(4 * sizeof(double));
		//	LRVEC[treeNum][i] = (double *)malloc(4 * sizeof(double));
		//}

		A[0][1]=pi[treeNum][1]*par[treeNum][0];
		A[0][2]=pi[treeNum][2]*par[treeNum][1];
		A[0][3]=pi[treeNum][3]*par[treeNum][2];
		A[1][0]=pi[treeNum][0]*par[treeNum][0];
		A[1][2]=pi[treeNum][2]*par[treeNum][3];
		A[1][3]=pi[treeNum][3]*par[treeNum][4];
		A[2][0]=pi[treeNum][0]*par[treeNum][1];
		A[2][1]=pi[treeNum][1]*par[treeNum][3];
		A[2][3]=pi[treeNum][3]*par[treeNum][5];
		A[3][0]=pi[treeNum][0]*par[treeNum][2];
		A[3][1]=pi[treeNum][1]*par[treeNum][4];
		A[3][2]=pi[treeNum][2]*par[treeNum][5];


		//printf("GTR Free Parameter: %lf, %lf, %lf, %lf, %lf, %lf\n", par[treeNum][0], par[treeNum][1], par[treeNum][2], par[treeNum][3], par[treeNum][4],par[treeNum][5]);

		//printf("-q_ii\n");
		for (i=0; i<4; i++)
		{
			// Can't you just subtract directly from A[i][i]??
			A[i][i]=0.0;
			sum=0.0;
			for (j=0; j<4; j++)
				sum = sum + A[i][j];
			A[i][i] = -sum;
			//printf("%d\t%lf\n", i, A[i][i]);
		}

		//printf("Pi - 0: %lf\t1: %lf\t2: %lf\t3: %lf\n", pi[treeNum][0], pi[treeNum][1], pi[treeNum][2], pi[treeNum][3]);

		// Edit here to get Q (stored as A[][]) and pi/frequencies to calculate the estimated number of differences
		if (eigen(1, A[0], 4, RRVAL[treeNum], RIVAL, RRVEC[treeNum], RIVEC[0], workspace) != 0)
		{
			printf("Transitions matrix did not converge or contained non-real values!\n");
			exit(-1);
		}
		for (i=0; i<4; i++){
			for (j=0; j<4; j++){
				LRVEC[treeNum][i * 4 +j] = RRVEC[treeNum][i * 4 + j];
				//printf("LRVEC[%d][%d][%d]: %lf\n", treeNum, i, j, LRVEC[treeNum][i * 4 + j]);
			}
		}
		if (matinv(RRVEC[treeNum],4, 4, workspace) != 0)	//TO DO: MAKE SURE THE RRVEC INPUT IS CORRECT FOR THIS FUNCTION
			printf("Could not invert matrix!\nResults may not be reliable!\n");
	}
}


/*some old code for reading a fasta file and storing DNA as ints*/
// Probably needs same fix as get_frac_like.c's readseq
int readseq(int *nb, int treeNum)
{
	int i, j, num, k=0;
	char c;

	fscanf(infile,"%i %i",&num,nb);
	//printf("there are %i species and %i bases in tree %i\n",num,*nb, treeNum);
	if (num < 3) {printf("This is for more than two seqs only!\n"); exit(-1);}
	DATA[treeNum] = (int**)malloc(num*(sizeof(int*)));
	for (i=0; i<num; i++)
		DATA[treeNum][i] = (int*)malloc(*nb*(sizeof(int)));
	while ((c=(fgetc(infile)))!='\n');
	do
	{
		for (i=0; i<num; i++)
		{
			j=k;
			while ((c=tolower((fgetc(infile))))!='\n')
			{
				if ((c!=' ')&&(c!='\t')) {
					if (c=='a') DATA[treeNum][i][j] = 0;
					else if (c=='c') DATA[treeNum][i][j] = 1;
					else if (c=='g') DATA[treeNum][i][j]  = 2;
					else if (c=='t') DATA[treeNum][i][j]  = 3;
					else if (c=='n'||c=='-'||c=='~') DATA[treeNum][i][j]  = -1;
					else{
						printf("\nBAD BASE (%c) in species %i base %i",c,i+1,j+1);
						scanf("%i",&i);
						exit(-1);
					}
					j++;
					if (i==(num-1)) k++;
				}
			}		
		}
	} while (k<*nb);
	return num;
}

int isblankorreturn(char c)
{
	if (c<33) return 1;
	else return 0;
}

// Probably needs same fix as get_frac_like.c
int read_query_data(int num)
{
	int i, j, v, length, k=0, treeNum, assignment;
	char c;

	for (i=0; i<num; i++)
	{
		fscanf(infile,"%i %i ",&length, &startpos[i]);
		readlength[i]=length;
		QUERYDATA[i]=malloc(length*(sizeof(int)));
		for (j=0; j<length; j++){
			do {
				fscanf(infile,"%c",&c);
			} while (isblankorreturn(c)==1);
			c = tolower(c);
			if (c=='a') QUERYDATA[i][j] = 0;
			else if (c=='c') QUERYDATA[i][j] = 1;
			else if (c=='g') QUERYDATA[i][j]  = 2;
			else if (c=='t') QUERYDATA[i][j]  = 3;
			else if (c=='n'||c=='-'||c=='~') QUERYDATA[i][j]  = -1;
			else{
				printf("\nBAD BASE (%c) in query %i base %i",c,i+1,j+1);
				exit(-1);
			}
		}
		//printf("\n");
	}
	return num;
	//printf("Errors setting of %d, sequence %d, assignment %d\n", errors, seq, node);
} 



void get_fractionalike(int treeNum)
{
	int i, j, k, v, inin;
	double a;
	char c;

	//fractional likelihoods. dimension: [numbase][2*numseq-1][NUMCAT][4 (leaf nodes) or 8 (internal nodes)]
	FRACLIKE[treeNum] = (double ****)malloc(numbase*(sizeof(double***)));
	for (i=0; i<numbase; i++){
		FRACLIKE[treeNum][i] = (double ***) malloc((2 * numseq[treeNum] - 1) * (sizeof(double**)));
		for (j=0; j < 2 * numseq[treeNum] - 1; j++){
			FRACLIKE[treeNum][i][j] = (double **)malloc( NUMCAT * (sizeof(double*)));
			if (j<numseq[treeNum]){
				for (k=0;k<NUMCAT; k++)
					FRACLIKE[treeNum][i][j][k] = (double *)malloc(4 * (sizeof(double)));
			}
			else {
				for (k=0;k<NUMCAT; k++)
					FRACLIKE[treeNum][i][j][k] = (double *)malloc( 8  *(sizeof(double)));
			}
		}
	}
	//read in the fractional likelihoods
	for (i=0; i<NUMCAT; i++)
	{
		do {
			fscanf(infile,"%c",&c);
		} while (isblankorreturn(c)==1);
		if (c!='C') {
			printf("error reading fractional likelihoods (C%i: %c != C)\n",i,c);
			exit(-1);
		}
		fscanf(infile,"%i",&inin);
		if (inin!=i+1) {
			printf("error reading fractional likelihoods (c%i: %i != %i)\n",i,inin,i+1);
			exit(-1);
		}
		for (j=0; j<numbase; j++){
			do {
				fscanf(infile,"%c",&c);
			} while (isblankorreturn(c)==1);
			if (c!='S') {
				printf("error reading fractional likelihoods (s%i: '%c' != S)\n",j,c);
				exit(-1);
			}
			fscanf(infile,"%i",&inin);
			if (inin!=j+1) {
				printf("error reading fractional likelihoods(s%i: %i != %i)\n",j,inin,j);
				exit(-1);
			}
			do {
				fscanf(infile,"%c",&c);
			} while (isblankorreturn(c)==1);
			if (c!=':') {
				printf("error reading fractional likelihoods (s%i: '%c' != :)\n",j,c);
				exit(-1);
			}

			// leaf nodes first
			for (k=0; k < numseq[treeNum]; k++)
			{
				for (v=0; v<4; v++){
					fscanf(infile,"%lf",&a);
					FRACLIKE[treeNum][j][k][i][v]=a;
				}
			}

			// non-leaf nodes with 8
			for (k = numseq[treeNum]; k < 2 * numseq[treeNum] - 1; k++)
			{
				//For fractional likelihood
				for (v=0; v<4; v++){
					fscanf(infile,"%lf",&a);
					FRACLIKE[treeNum][j][k][i][v]=a;
				}

				//For conditional likelihood
				for (v=0; v<4; v++){
					fscanf(infile,"%lf",&a);
					FRACLIKE[treeNum][j][k][i][v+4]=a;
				}

			}

		}
	}
}

// How likely read is not error
// make default error 0.01
void make_readfraclike()
{
	int i, j, k, pos, read, parseNum;
	double e, ec, err0, err1, err2, err3;

	// TO DO: fix this with more reasonable general errors
	// Should probably make this an option for users, so need to update
	e = log(0.01);	//Error for parts of sequence not covered by errorProfile
	ec = log(1.0000 - 0.01);

	readlike = malloc(numquery*(sizeof(double**)));

	// All reads get general error
	for (i = 0; i < numquery; i++)
	{
		readlike[i] = malloc(readlength[i]*(sizeof(double*)));
		for (j=0; j<readlength[i]; j++){
			if (QUERYDATA[i][j] != -1)
			{
				readlike[i][j] = malloc(4*(sizeof(double)));
				for (k=0; k<4; k++)
				{
					if (k==QUERYDATA[i][j])
					{
						//readlike[i][j][k] = ec + errorTestLog;
						readlike[i][j][k] = ec;
					}
                        	        else
					{
						//readlike[i][j][k] = e + errorTestLog;
						readlike[i][j][k] = e;
					}
				}
			}
		}
	}

	// fix error based on error profile
	// Should be more flexible to custom error profiles even if not the most efficient
	//while (fscanf(infile,"%d %d %lf %lf %lf %lf", &read, &pos, &err0, &err1, &err2, &err3) == 6)
	//{
	//	if(QUERYDATA[read][pos] != -1)
	//	{
	//		readlike[read][pos][0] = err0;
	//		readlike[read][pos][1] = err1;
	//		readlike[read][pos][2] = err2;
	//		readlike[read][pos][3] = err3;
	//	}
	//	else
	//	{
	//		printf("Warning, error profile includes positions not in query alignment. Please review error profile, but ratePlacer is proceeding.\n");
	//	}
    	//}	

	////// Use for a test function in the future
	////for (i = 0; i < numquery; i++)
	////{
	////	for (j = 0; j < readlength[i]; j++)
	////	{
	////		if (QUERYDATA[i][j] != -1)
	////		{
	////			printf("Read %d pos %d QUERY %d a: %f c: %f g: %f t: %f\n", i, j, QUERYDATA[i][j], readlike[i][j][0], readlike[i][j][1], readlike[i][j][2], readlike[i][j][3]);
	////		}
	////		else
	////		{
	////			printf("Read %d pos %d QUERY %d no profile\n", i, j, QUERYDATA[i][j]);
	////		}
	////	}
	////}
	//
	//// Last line read in was not properly read in
	//if(fscanf(infile, "%d", &read) != EOF)
	//{
	//	printf("Error in reading in error profile\n");
	//	exit(0);
	//}
}

//numseq: number of sequences in reference data set
//numquery: number of sequences in the query data
//numbase: total length of alignment of reference sequences
//NUMCAT: number of categories in the gamma distribution 
void read_data(char *assignfile, char *fraclikefile, char *querydatafile, char *referencedatafile, char *errorfile)
{
	int i, j, k, v, nin, numcat;
	double a, b, checksum, MINLIKE=-INF;
	char tempFileName[500], strTree[500];
	
	usedTrees = (int *)calloc(numTrees, sizeof(int));
	if(usedTrees == NULL)
	{
		printf("ERROR: Malloc failure for usedTrees\n");
		exit(0);
	}

	//This may be a good one to make into single malloc at some point since its a known size
	pi = (double **)malloc(numTrees * sizeof(double*));
	par = (double **)malloc(numTrees * sizeof(double*));	

	DATA = (int***)malloc(numTrees*(sizeof(int**)));
	bls = (double **)malloc(numTrees * sizeof(double*));
	nodeages = (double **)malloc(numTrees * sizeof(double *));
	FRACLIKE = (double *****)malloc(numTrees *(sizeof(double****)));
	maxAges = (double *)malloc(numTrees * sizeof(double));
	statevector = (double **)malloc(numTrees * sizeof(double*));
	nodeOrder = (int **)malloc(numTrees * sizeof(int*));

	totMaxAge = -INF;

	//READING IN ASSIGNMENTS
	if (VERBOSE) printf("Reading in species assignments\n");
	if (NULL==(infile=fopen(assignfile,"r")))
	{
		puts ("Cannot open infile with assignments: ");
		exit(-1);
	}
	fscanf(infile,"%i\n",&numquery); 

	assignments = (int *)malloc(numquery*(sizeof(int)));
	treeAssign = (int *)malloc(numquery * (sizeof(int)));
	assignAges = (double *)malloc(numquery * (sizeof(double)));

	// Determine how nodes are stored and determine if script is needed to convert node names etc
	for (i=0; i<numquery; i++){
		//read in tree assignment first
		fscanf(infile,"%i",&v);
		treeAssign[i] = v;
		usedTrees[v] = 1;
		if(treeAssign[i] < 0 || treeAssign[i] >= numTrees)
		{
			printf("Error reading assignments for trees with %d\n", treeAssign[i]);
			exit(-1);
		}
		//node assignment
		fscanf(infile,"%i",&v); 
		assignments[i] = v-1;   //notice that we here convert from counting from 1 to counting from 0
	}

	fclose(infile);

	for(int treeNum = 0; treeNum < numTrees; treeNum++)
	{
		if(usedTrees[treeNum] == 0)
		{
			// no reads were assigned to this tree, don't need to read in or allocate more memory for this tree
			continue;
		}
		//READING IN OUTPUT FROM GET_FRAC_LIKE WITH RACTIONAL LIKELIHOODS, NODE AGES, AND MORE
		sprintf(strTree, "%01d", treeNum);
		tempFileName[0] = '\0';
		strcat(tempFileName, fraclikefile);
		strcat(tempFileName, "/");
		strcat(tempFileName, strTree);
		strcat(tempFileName, "_likelihood.txt");
		if (NULL==(infile=fopen(tempFileName,"r")))
		{
			printf("Cannot open infile with fractional likelihoods: %s\n", tempFileName);
			exit(-1);
		}
		//line 1
		// THERE IS AN ASSUMPTION THAT NUMBASE AND NUMCAT ALWAYS THE SAME, SHOULD MAKE A CHECK FOR THAT
		fscanf(infile,"%i %i %i\n",&numseq[treeNum], &numbase, &numcat);
		if (VERBOSE) printf("Reading in fractional likelihoods for %i sequences\n",numseq[treeNum]);
		if (NUMCAT != numcat){
			printf("Wrong number of categories for the discretization of the gamma distribution");
			exit(-1);
		}
		statevector[treeNum] = malloc(NUMCAT*(sizeof(double)));
		//line 2
		for (i=0;i<NUMCAT; i++)
			fscanf(infile,"%lf ",&statevector[treeNum][i]);

		nodeages[treeNum] = (double *)malloc((2 * numseq[treeNum] - 1) * (sizeof(double)));  //nodeages: ages of internal nodes in reference data
		bls[treeNum] = (double *)malloc((2 * numseq[treeNum] - 1) * (sizeof(double))); //bls: branch lengths associated with each node
		nodeOrder[treeNum] = (int *)malloc((2 * numseq[treeNum] - 1) * (sizeof(int))); //nodeOrder: order of nodes in likelihood file reflect order based on nodeage + bls

		// line 3 - numseq + 3, node ages
		for (i=0; i<2*numseq[treeNum]-1; i++){
			fscanf(infile,"%i",&j);
			if (j<0 || j>2*numseq[treeNum]-2) {
				printf("Error reading node ages");
				exit(-1);
			}
			if (fgetc(infile) != ':') {
				printf("Error reading node ages");
				exit(-1);}
			else { 
				fscanf(infile,"%lf",&a);
				fscanf(infile,"%lf",&b);
				nodeages[treeNum][j] = a;
				bls[treeNum][j] = b;
				nodeOrder[treeNum][i] = j;
				//printf("%d\t%.16f\t%.16f\n", j, a, b);
			}
		}

		fscanf(infile, "%lf", &maxAges[treeNum]);
		//printf("%.16f\n", maxAges[treeNum]);
		//if(treeNum > 0 && maxAges[treeNum] > totMaxAge)
		if(maxAges[treeNum] > totMaxAge)
		{
			totMaxAge = maxAges[treeNum];
		}

		get_fractionalike(treeNum); //reads in all the fractional likelihoods
		fclose(infile);

		//READING IN REFRENCE SEQUENCE DATA
		tempFileName[0] = '\0';
		strcat(tempFileName, referencedatafile);
		strcat(tempFileName, "/");
		strcat(tempFileName, strTree);
		strcat(tempFileName, "_reference.txt");
		if (VERBOSE) printf("Reading in reference sequences: ");
		if (NULL==(infile=fopen(tempFileName,"r")))
		{
			printf("Cannot open infile with reference sequence data");
			exit(-1);
		}

		// readseq returns num nodes and alters input with number of bases
		i = readseq(&j, treeNum);
		if (i!=numseq[treeNum]|| j!=numbase)
		{
			printf("Number of sequences and sequencelengths do no match in input files\n");
			exit(-1);
		}
		fclose(infile);

		pi[treeNum] = (double*)malloc(4 * sizeof(double));
		par[treeNum] = (double*)malloc(6 * sizeof(double));

		//READING IN GTR parameters DATA
		tempFileName[0] = '\0';
		strcat(tempFileName, referencedatafile);
		strcat(tempFileName, "/");
		strcat(tempFileName, strTree);
		strcat(tempFileName, "_parameter.txt");
		if (VERBOSE) printf("Reading in GTR+Gamma parameters\n");
		if (NULL==(infile=fopen(tempFileName,"r")))
		{
			printf("Cannot open infile with GTR+Gamma parameters");
			exit(-1);
		}
		checksum=0.0;
		// line 1
		fscanf(infile,"%lf",&a);

		// line 2: nucleotide frequencies
		for (i=0; i<4; i++){
			fscanf(infile,"%lf",&pi[treeNum][i]);
			checksum+=pi[treeNum][i];
		}

		if (checksum > 1.0 + 3e-7 || checksum < 1.0 - 3e-7){	//TO DO: bring back to -8 and figure out why tolerances are not being met with COI database
			printf("WARNING: Nucleotide frequencies not properly scaled or not read correctly (checksum: %.15f)\n",checksum);
			exit(-1);
		}

		for (i=0; i<6; i++)
		{
			fscanf(infile,"%lf",&par[treeNum][i]);
		}

		fclose(infile);

		// printf("Printing node ages for tree %d\n", treeNum);

		// for(int i = 0; i < numseq[treeNum]; i++)
		// {
		// 	printf("\tNode %d: %.16f\n", i, nodeages[treeNum][i]);
		// }
	}


	//READING IN QUERY DATA
	if (VERBOSE) printf("Reading in query data: ");
	if (NULL==(infile=fopen(querydatafile,"r")))
	{
		puts ("Cannot open infile with query data\n");
		exit(-1);
	}
	// Line 1, number of queries
	fscanf(infile,"%i",&numquery);
	printf("There are %i query sequences\n",numquery);
	QUERYDATA = malloc(numquery*(sizeof(int*)));
	startpos = malloc(numquery*(sizeof(int)));
	readlength = malloc(numquery*(sizeof(int)));

	// This expression also reads in data to QUERYDATA
	if (read_query_data(numquery) != numquery)
	{
		printf(" Different number of query sequences found in query data file and in assignment file\n");      
		exit(-1);
	}	

	////This code below needs to be modified if queryages can be specified
	//// Ask for future I guess
	//queryagesknown=0;

	//if (queryagesknown)
	//{
	//	nodeages = malloc(numquery*(sizeof(double)));
	//	for (i=0; i<numquery; i++)
	//		fscanf(infile,"%lf",&nodeages[i]);
	//}
	fclose(infile);

	
	//Testing if read assignments are valid
	//Here because we need to know number of references in the tree
	if (VERBOSE) printf("Testing validity of node assignments.\n");
	for (i=0; i<numquery; i++){
		if (assignments[i]<0 || assignments[i] > 2*numseq[treeAssign[i]]-1)
		{
			printf("Error reading assignments for read %d with assignment %d from tree %d with %d references\n", i, assignments[i], treeAssign[i], numseq[treeAssign[i]]);
			exit(-1);    
		}
	}



	//READING IN ERROR
	if (VERBOSE) printf("Reading in error profile\n");
	if (NULL==(infile=fopen(errorfile,"r")))
	{
		puts ("Cannot open infile with error profile: ");
		exit(-1);
	}

	make_readfraclike();

	fclose(infile);
	//exit(0);
}

// Recursive search of children subtrees for maximum likelihood assignment
void searchChildren(double p[3], int *L, double *L_lik, int root)
{
	//p[0]: readNum
	//p[1]: treeNum
	//p[2]: curNode
	int nfun, testNodes[3];
	double invector[3], lowbound[3], upbound[3], eh0=3e-8, testLik;

	testNodes[0] = (int)p[2];

	getGFLChildren((int)p[2], testNodes, (int)p[1]);

	//Search left subtree
	if(testNodes[1] < numseq[(int)p[1]])
	{
		//Left child is leaf
		p[2] = testNodes[1];

		invector[1] = 0.5;
		invector[2] = bls[(int)p[1]][testNodes[1]]/2.0;
		lowbound[1] = eh0;
		lowbound[2] = eh0;
		upbound[1] = 1.0-eh0;
		upbound[2] = bls[(int)p[1]][testNodes[1]]-eh0;
		nfun=0;

		// TO DO: Make a version of this that only tests like 100 iterations for speed and compare the two versions to see if they agree or not!!
		testLik = findmax_amoeba_limited(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation, p, 3);

		printf("\tL of node %d is %.16f with a: %.16f, root: %.16f given age: %.16f bls: %.16f, est age: %.16f\n", testNodes[1], testLik, invector[1], invector[2], nodeages[(int)p[1]][testNodes[1]], bls[(int)p[1]][testNodes[1]], (1.0-invector[1])*(nodeages[(int)p[1]][testNodes[1]]+invector[2]));

		if(testLik < *L_lik)
		{
			*L = testNodes[1];
			*L_lik = testLik;
		}
	}
	else
	{
		//Left child is internal
		p[2] = (double)testNodes[1];
		searchChildren(p, L, L_lik, root);
	}
	
	//Search right subtree
	if(testNodes[2] < numseq[(int)p[1]])
	{
		//Right child is leaf
		p[2] = testNodes[2];

		invector[1] = 0.5;
		invector[2] = bls[(int)p[1]][testNodes[2]]/2.0;
		lowbound[1] = eh0;
		lowbound[2] = eh0;
		upbound[1] = 1.0-eh0;
		upbound[2] = bls[(int)p[1]][testNodes[2]]-eh0;
		nfun=0;

		// TO DO: Make a version of this that only tests like 100 iterations for speed and compare the two versions to see if they agree or not!!
		testLik = findmax_amoeba_limited(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation, p, 3);

		printf("\tL of node %d is %.16f with a: %.16f, root: %.16f given age: %.16f bls: %.16f, est age: %.16f\n", testNodes[2], testLik, invector[1], invector[2], nodeages[(int)p[1]][testNodes[2]], bls[(int)p[1]][testNodes[2]], (1.0-invector[1])*(nodeages[(int)p[1]][testNodes[2]]+invector[2]));

		if(testLik < *L_lik)
		{
			*L = testNodes[2];
			*L_lik = testLik;
		}
	}
	else
	{
		//Right child is internal
		//Right child is internal
		p[2] = (double)testNodes[2];
		searchChildren(p, L, L_lik, root);
	}

	//Test curNode
	//Skip root
	if (testNodes[0] != root)
	{
		p[2] = testNodes[0];

		invector[1] = 0.5;
		invector[2] = bls[(int)p[1]][testNodes[0]]/2.0;
		lowbound[1] = eh0;
		lowbound[2] = eh0;
		upbound[1] = 1.0-eh0;
		upbound[2] = bls[(int)p[1]][testNodes[0]]-eh0;
		nfun=0;

		// TO DO: Make a version of this that only tests like 100 iterations for speed and compare the two versions to see if they agree or not!!
		testLik = findmax_amoeba_limited(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation, p, 3);

		printf("\tL of node %d is %.16f with a: %.16f, root: %.16f given age: %.16f bls: %.16f, est age: %.16f\n", testNodes[0], testLik, invector[1], invector[2], nodeages[(int)p[1]][testNodes[0]], bls[(int)p[1]][testNodes[0]], (1.0-invector[1])*(nodeages[(int)p[1]][testNodes[0]]+invector[2]));

		if(testLik < *L_lik)
		{
			*L = testNodes[0];
			*L_lik = testLik;
		}
	}
}

// Search to root for maximum likelihood assignment
// Tests current node too
void searchLineage(double p[3], int *L, double *L_lik, int root)
{
	int nfun, treeNum = p[1], testNode = p[2];
	double invector[3], lowbound[3], upbound[3], eh0=3e-8, testLik;

	while(testNode != root)
	{
		invector[1] = 0.5;
		invector[2] = bls[treeNum][testNode]/2.0;
		lowbound[1] = eh0;
		lowbound[2] = eh0;
		upbound[1] = 1.0-eh0;
		upbound[2] = bls[treeNum][testNode]-eh0;
		nfun=0;

		// TO DO: Make a version of this that only tests like 100 iterations for speed and compare the two versions to see if they agree or not!!
		testLik = findmax_amoeba_limited(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation, p, 3);

		printf("\tL of node %d is %.16f with a: %.16f, root: %.16f given age: %.16f bls: %.16f, est age: %.16f\n", testNode, testLik, invector[1], invector[2], nodeages[treeNum][testNode], bls[treeNum][testNode], (1.0-invector[1])*(nodeages[treeNum][testNode]+invector[2]));

		if(testLik < *L_lik)
		{
			*L = testNode;
			*L_lik = testLik;
		}

		testNode = getGFLPar(testNode, treeNum);
		printf("Testing next node %d when max is %d\n", testNode, *L);
	}

}

//This function will test for the best assignment of reads since tronko returns node that is maximal, but not the edge
void bestAssignment(int root, int treeNum)
{
	int i, L1;
	double p[3], L1_lik;

	printf("Testing assignments\n");

	for (i = 0; i < numquery; i++)
	{
		//Slow, but better than keeping all trees into memory?
		if(treeAssign[i] != treeNum)
		{
			continue;
		}

		printf("Sequence %d of tree %d and node %d\n", i, treeNum, assignments[i]);

		p[0] = i;	//read num does not change
		p[1] = treeAssign[i];	//tree num does not change
		p[2] = assignments[i]; 
		L1 = assignments[i];
		L1_lik = INFINITY;

		//Tronko assignment is the root
		if(p[2] == root)
		{
			printf("Root\n");
			searchChildren(p, &L1, &L1_lik, root);

			assignments[i] = L1;
			assignAges[i] =  nodeages[treeNum][L1] + bls[treeNum][L1];
		}
		//Tronko assignment is a leaf
		else if(p[2] < numseq[treeNum])
		{
			printf("Leaf\n");
			searchLineage(p, &L1, &L1_lik, root);

			assignments[i] = L1;
			assignAges[i] =  nodeages[treeNum][L1] + bls[treeNum][L1];
		}
		else
		{
			printf("Internal\n");
			searchChildren(p, &L1, &L1_lik, root);
			p[2] = assignments[i]; //So don't double up on nodes
			searchLineage(p, &L1, &L1_lik, root);

			assignments[i] = L1;
			assignAges[i] =  nodeages[treeNum][L1] + bls[treeNum][L1];
		}
	}
}

void maximize_like_seperately_for_all2D_Print(double **par)
{
	printf("Starting maximize seperately for all\n");
	int i, k, v, nfun;
	double p[3];
	double L1, invector[3],lowbound[3], upbound[3], eh0=3e-8;

	onDindic=0;


	// assignmentMode of 0 means single assignment given
	for (i=0; i<numquery; i++)
	{
		p[0] = i;
		p[1] = treeAssign[i];
		p[2] = assignments[i];

		invector[1] = 0.5;
		invector[2] = bls[treeAssign[i]][assignments[i]]/2.0;
		lowbound[1] = eh0;
		lowbound[2] = eh0;
		upbound[1] = 1.0-eh0;
		upbound[2] = bls[treeAssign[i]][assignments[i]]-eh0;
		nfun=0;

		printf("Sequence %d\n", i);
		L1 = findmax_amoeba(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation, p, 3);	
		printf("\tParameter estimates %.16f %.16f: %.16f\n", invector[1],invector[2],L1);
		printf("\tAssignment %d of tree %d age: %.16f\n", assignments[i], treeAssign[i], nodeages[treeAssign[i]][assignments[i]]);
		printf("\tsequence age: %.16f\n",(1.0-invector[1])*(nodeages[treeAssign[i]][assignments[i]]+invector[2]));
		printf("Site scores:\n");
		getlike_gamma_root_in_trifurcation_Print(invector,p);
	}
}

void maximize_like_seperately_for_all2D(double **par)
{
	printf("Starting maximize seperately for all\n");
	int i, k, v, nfun;
	double p[3];
	double L1, invector[3],lowbound[3], upbound[3], eh0=3e-8;

	onDindic=0;


	// assignmentMode of 0 means single assignment given
	for (i=0; i<numquery; i++)
	{
		p[0] = i;
		p[1] = treeAssign[i];
		p[2] = assignments[i];

		invector[1] = 0.5;
		invector[2] = bls[treeAssign[i]][assignments[i]]/2.0;
		lowbound[1] = eh0;
		lowbound[2] = eh0;
		upbound[1] = 1.0-eh0;
		upbound[2] = bls[treeAssign[i]][assignments[i]]-eh0;
		nfun=0;

		L1 = findmax_amoeba(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation, p, 3);
		printf("Sequence %d\n", i);
		printf("\tParameter estimates %.16f %.16f: %.16f\n", invector[1],invector[2],L1);
		printf("\tAssignment %d age: %.16f\n", assignments[i], nodeages[treeAssign[i]][assignments[i]]);
		printf("\tsequence age: %.16f\n",(1.0-invector[1])*(nodeages[treeAssign[i]][assignments[i]]+invector[2]));
	}
}

void likelihoodratiotest_for_all(double **par)
{
	printf("Starting maximize seperately with likelihood ratio test for all\n");
	int i, k, v, nfun;
	double p[3], L1, L2;
	double invector[3],lowbound[3], upbound[3], eh0=3e-8;

	// assignmentMode of 0 means single assignment given
	for (i=0; i<numquery; i++)
	{
		p[0] = i;
		p[1] = treeAssign[i];
		p[2] = assignments[i];

		invector[1] = 0.5;
		invector[2] = bls[treeAssign[i]][assignments[i]]/2.0;
		lowbound[1] = eh0;
		lowbound[2] = eh0;
		upbound[1] = 1.0-eh0;
		upbound[2] = bls[treeAssign[i]][assignments[i]]-eh0;
		//printf("Checking bls %.16f and node age %.16f\n", upbound[2], upbound[1]);
		nfun=0;
		onDindic=0;
		L1 = findmax_amoeba(invector,lowbound, upbound, 2, getlike_gamma_root_in_trifurcation, p, 3);

		printf("Sequence %d\n", i);
		printf("\tParameter estimates %.16f %.16f (%.16f): %.16f\n", invector[1],invector[2], bls[treeAssign[i]][assignments[i]],L1);
		printf("\tAssignment %d age of tree %d: %.16f\n", assignments[i], treeAssign[i], nodeages[treeAssign[i]][assignments[i]]);
		printf("\tsequence age: %.16f\n",(1.0-invector[1])*(nodeages[treeAssign[i]][assignments[i]]+invector[2]));

		invector[1] = bls[treeAssign[i]][assignments[i]]/2.0;
		lowbound[1] = eh0;
		upbound[1] = bls[treeAssign[i]][assignments[i]]-eh0;
		invector[2]=-1;
		nfun=0;
		onDindic=1;

		// NEEDS TO BE SWITCHED TO GOLDEN SECTION
		L2 = findmax_amoeba(invector,lowbound, upbound, 1, getlike_gamma_root_in_trifurcation_L0, p, 3); 
		
		printf("\tParameter estimates llr %.16f: %.16f\n", invector[1],L2);
		printf("\tLikelihood ratio statistic: %.16f\n",2.0*(L1-L2));
	}
}

// Will "drop reads" by changing readStart based on timeInc and test age
// New update: Will not be based on the node ages themselves as we're now making
// this for across trees
double dropReads(int *readStart, double timeInc)
{
	int readDropped = 0;
	double test = testAge;

	do{
		test += timeInc;
		//TO DO: Try to make assignAges work here since that would be faster!! 
		while(*readStart < numquery-1 && nodeages[treeAssign[usedReads[*readStart]]][assignments[usedReads[*readStart]]] + bls[treeAssign[usedReads[*readStart]]][assignments[usedReads[*readStart]]] <= test) //similar -2 as above
		//while(*readStart < numquery && assignAges[assignments[usedReads[*readStart]]] <= test)
		{
			*readStart = *readStart + 1; 
			readDropped = 1;
			//printf("Read dropped!\n");
		}
		//if(*readStart >= numquery-1)
		//{
		//	readDropped = 1;
		//}
		//printf("\tTesting %f\n", test);
		//test = nodeages[nodeOrder[*nodePointer+1]] + bls[nodeOrder[*nodePointer+1]];
	}while(*readStart < numquery-1 && readDropped == 0);


	return nodeages[treeAssign[usedReads[*readStart]]][assignments[usedReads[*readStart]]] + bls[treeAssign[usedReads[*readStart]]][assignments[usedReads[*readStart]]];


	//int readDropped = 0;
	////double test = testAge;
	//double test = assignAges[*readStart];

	//do{
	//	test += timeInc;
	//	//printf("\tTesting: %lf\n", test);
	//	while(*readStart < numquery && assignAges[usedReads[*readStart]] <= test)
	//	{
	//		*readStart = *readStart + 1;
	//		readDropped = 1;
	//	}
	//}while(readDropped == 0 && *readStart < numquery);

	////printf("New test age of %lf\n", test);

	//if(readDropped == 0)
	//{
	//	return(-1);
	//}
	//else
	//{
	//	//return(test);
	//	return(assignAges[*readStart]);
	//}
	
	//int readDropped = 0;
	////int checkNumNodesDropped = 0;
	//double test =  nodeages[nodeOrder[*nodePointer]] + bls[nodeOrder[*nodePointer]] + timeInc;

	////printf("Test age with node %d is %.16f\n", *nodePointer, test);
	//do{
	//	//Remove some nodes so that the age changes	
	//	while(*nodePointer < 2 * numseq[treeNum] - 1 && nodeages[nodeOrder[*nodePointer]] + bls[nodeOrder[*nodePointer]] <= test)	//makes sure at least one node remains...
	//	{
	//		*nodePointer = *nodePointer + 1;
	//		//checkNumNodesDropped++;
	//	}

	//	test = nodeages[nodeOrder[*nodePointer]] + bls[nodeOrder[*nodePointer]];

	//	//printf("Test age with node %d is %.16f\n", *nodePointer, test);

	//	while(*readStart < numquery && nodeages[assignments[usedReads[*readStart]]] + bls[assignments[usedReads[*readStart]]] <= test)
	//	{
	//		*readStart = *readStart + 1; 
	//		readDropped = 1;
	//		//printf("Read dropped!\n");
	//	}

	//	if(*readStart == numquery)
	//	{
	//		printf("Unable to find maximum bound while keeping any reads.\nExiting\n");
	//		exit(0);
	//	}	
	//}while(readDropped == 0 && *nodePointer < 2 * numseq[treeNum] - 1 && *readStart < numquery);
	//// || nodeages[nodeOrder[*nodePointer]] + bls[nodeOrder[*nodePointer]] > test);

	//return nodeages[nodeOrder[*nodePointer]] + bls[nodeOrder[*nodePointer]];
}

//merge functionality, using in-place merge sort to reduce memory load
// Will base merge on nodeAge[i]+nodebl[i] to get age in which reads assigned to i will need to cut in ratePlacer
// for reference: https://www.geeksforgeeks.org/merge-sort/ and https://www.geeksforgeeks.org/in-place-merge-sort/
void merge(int left, int right, int mid)
{
	// If already merged, are there cases where this may be true but not sorted?
	if (assignAges[usedReads[mid]] <= assignAges[usedReads[mid+1]])
	{
		return;
	}

	int midstart = mid + 1;

	while(left <= mid && midstart <= right)
	{
		// if elements are in right place (ie left < right)
		if(assignAges[usedReads[left]] <= assignAges[usedReads[midstart]])
		{
			left++;
		}
		else
		{
			int index = midstart, tempOrder = usedReads[midstart];
			//double tempAge = assignAges[midstart];

			// shift all elements between left and midstart to right by 1
			while (index != left)
			{
				usedReads[index] = usedReads[index - 1];
				//assignAges[index] = assignAges[index - 1];
				index--;
			}
			usedReads[left] = tempOrder;
			//assignAges[left] = tempAge;

			// Update
			left++;
			mid++;
			midstart++;	// may be redundant, but removes some extra calc?
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
	usedReads = (int *)malloc(sizeof(int) * numquery);

	//TO DO: Determine if this is the most efficient way
	//or if I should have merge sort work on all data associated with
	//the reads instead
	for(int i = 0; i < numquery; i++)
	{
		usedReads[i] = i;
	}

	mergeSort(0, numquery-1);

	//for(int i = 0; i < numquery; i++)
	//{
	//	printf("%d is age %.16f and %.16f\n", i, assignAges[usedReads[i]], nodeages[treeAssign[usedReads[i]]][assignments[usedReads[i]]] + bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]]);
	//}
}

//// Initialization where reads are ordered by their assignment node order
//void orderReads()
//{
//	int arrayPointer = 0;
//
//	usedReads = (int *)malloc(sizeof(int) * numquery);
//	
//	//Should be numseq[treeNum],nodeOrder[treeNum] and running over trees, but currently testing on single tree
//	for (int i = 0; i < 2 * numseq[0] - 1; i++)
//	{
//		for (int j = 0; j < numquery; j++)
//		{
//			// if the assignment node matches
//			if (assignments[j] == nodeOrder[0][i])
//			{
//				usedReads[arrayPointer] = j;
//				arrayPointer++;
//			}
//		}
//	}
//}

//getlike_ages but for a single tree
double getlike_ages_tree(double times, double parameters[7])
{
	//age for optimization should be times, is redudant and should fix?
	testAge = times;

	if(testAge < 0.0)
	{
		return 1000000000.0;
	}

	// The following should be thought about because it would be inconvienent to redo for each age
	// May want to make userReads a global, but think after this is implemented and working
	// TO DO TO DO TO DO!!!!
	int readStart = parameters[0], tree = parameters[1], i, nfun;
	double invector[3],lowbound[3], upbound[3], eh0=3e-8, age_like = 0.0, ageIncr, L2;
	double p[3];

	onDindic = 1;

	// Drop any reads that cannot be of test age
	for (i = readStart; i < numquery; i++)
	{
		if(treeAssign[usedReads[i]] != tree)
		{
			continue;
		}

		p[0] = usedReads[i];
		p[1] = treeAssign[usedReads[i]];
		p[2] = assignments[usedReads[i]];	

		nfun=0;
		invector[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]]/2.0;
		lowbound[1] = eh0;
		upbound[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]]-eh0;

		L2 =  GoldenSection(invector, lowbound, upbound, 1, getlike_gamma_root_in_trifurcation_testAge, p, 3);

		age_like += L2; //sum of log likelihoods
	}

	return(age_like);
}

double getlike_ages(double times, double parameters[7])
{
	//age for optimization should be times, is redudant and should fix?
	//printf("Testing age %.16f in getlike\n", times);
	testAge = times;

	if(testAge < 0.0)
	{
		return 1000000000.0;
	}

	// The following should be thought about because it would be inconvienent to redo for each age
	// May want to make userReads a global, but think after this is implemented and working
	// TO DO TO DO TO DO!!!!
	int readStart = parameters[0], i, nfun;
	double invector[3],lowbound[3], upbound[3], eh0=3e-8, age_like = 0.0, ageIncr, L2;
	double p[3];

	onDindic = 1;

	// Drop any reads that cannot be of test age
	for (i = readStart; i < numquery; i++)
	{
		//if(usedReads[i] != 85)
		//	continue;
		p[0] = usedReads[i];
		p[1] = treeAssign[usedReads[i]];
		p[2] = assignments[usedReads[i]];	

		nfun=0;
		invector[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]]/2.0;
		lowbound[1] = eh0;
		upbound[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]]-eh0;

		L2 =  GoldenSection(invector, lowbound, upbound, 1, getlike_gamma_root_in_trifurcation_testAge, p, 3);

		//printf("Read %d\tTime: %.16f\tLikelihood: %.16f\n", usedReads[i], testAge, L2);

		//printf("\tRead %d assigned to tree %d and node %d with likelihood %.16f with root placement of %.16f\n", usedReads[i], treeAssign[usedReads[i]], assignments[usedReads[i]], L2, nodeages[treeAssign[usedReads[i]]][assignments[usedReads[i]]]+invector[1]);


		age_like += L2; //sum of log likelihoods
	}
	//printf("Likelihood of age %.16f: %.16f\n", testAge, age_like);

	return(age_like);
}

//Gets Fisher information bounds for confidence interval
//finite method for second derivative
//(f(x+h) - 2f(x) + f(x-h))/h^2
// optAge is x, f(x) is optLik
// We are giving -lik, so remember (and test) to convert back again. Test if this is really needed
double confidenceIntervalFisher(double **par, double optAge, double optLik, int readStart)
{
	//Need to figure out good h!
	double h = optAge/100000, f_hx = 0.0, fx_h = 0.0, nSecondDeriv;
	int i, k, v, nfun;
	double p[3];
	double invector[3],lowbound[3], upbound[3], eh0=3e-8;

	// Drop any reads that cannot be of test age
	for (i = readStart; i < numquery; i++)
	{
		p[0] = usedReads[i];
		p[1] = treeAssign[usedReads[i]];
		p[2] = assignments[usedReads[i]];

		nfun=0;
		invector[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]]/2.0;
		lowbound[1] = eh0;
		upbound[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]]-eh0;

		//Does there need to be a test if this passes some internal node that reads need to be dropped in?
		testAge = optAge + h;
		f_hx += GoldenSection(invector,lowbound, upbound, 1, getlike_gamma_root_in_trifurcation_testAge, p, 3);

		invector[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]]/2.0;
		lowbound[1] = eh0;
		upbound[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]]-eh0;

		testAge = optAge - h;
		fx_h +=	GoldenSection(invector,lowbound, upbound, 1, getlike_gamma_root_in_trifurcation_testAge, p, 3); 
	}

	//printf("f(x-h): %.16f\nf(x+h):%.16f\nf(x) - f(x-h): %.16f\nf(x) - f(x+h): %.16f\nSec Deriv: %.16f\n", fx_h, f_hx, optLik - fx_h, optLik - f_hx, (2 * optLik - fx_h - f_hx)/(pow(h,2.0)));

	//equivalent to f_hx - 2*f_x + fx_h if getlike_gamma didn't return -Lik
	//Potentially implement as f_xh - 2 * f_x + fx_h so that -1 already incorporated to take the second derivate for
	//return((2 * optLik - fx_h - f_hx)/(pow(h,2.0)));
	
	//-1 * second derivative
	nSecondDeriv = (f_hx - 2 * optLik + fx_h)/(pow(h,2.0));

	//printf("%.16f\n", nSecondDeriv);

	return(1.96 / sqrt(nSecondDeriv));
}

void maximize_like_jointly_for_all2D(double **par, int allTrees)
{
	printf("Starting maximize jointly for all\n");
	int i, k, v, nfun, readStart = 0, nodePointer = 0, twice = 0, oldReadStart, oldNodePointer;
	double p[3], L1, L2, secD, confI;
	double invector[3],lowbound[3], upbound[3], eh0=3e-8, nextNodeAge, est_age, est_age_lik;
	double incr = totMaxAge/1000;

	//printf("Max age of %.16f and incr of %.16f\n", totMaxAge, incr);

	orderReads();

	// This nextNodeAge represents the max without dropping reads
	nextNodeAge = nodeages[treeAssign[usedReads[readStart]]][assignments[usedReads[readStart]]] + bls[treeAssign[usedReads[readStart]]][assignments[usedReads[readStart]]];

	onDindic = 1; 

	if (numquery > 1)
	{
		printf("Starting rough estimation for maximum bound\n");
		//double time2 = (double) clock()/CLOCKS_PER_SEC;
		// Rough optimization to find upper bound and reads to drop without running too much optimization
		// To change: Instread of making this per assigned node, do it via time slices to make a little faster
		// 		Also maybe a GoldenSection search that is a little less stringent?
		// 		Am just trying to do a rough optimization - maybe if time is still long despite time slice change
		do{
			//printf("Testing max age of %.16f\n", nextNodeAge);
			//add fillers
			p[0] = readStart;
			p[1] = nodePointer;
			p[2] = 0;

			//Should consider the best way to pick these values...
			L1 = getlike_ages(nextNodeAge - nextNodeAge/100.0, p); 
			L2 = getlike_ages(nextNodeAge - nextNodeAge/10.0, p);
			//printf("Results of L1 %.16f and L2 %.16f with max at %.16f with readStart at %d\n", L1, L2, nextNodeAge, readStart);
			// Likelihood near bound is better, so drop and test next age range
			if(L1 <= L2)
			{
				// TO DO: scale time increase by maxAge
				oldNodePointer = nodePointer;
				oldReadStart = readStart;
				testAge = nextNodeAge;
				nextNodeAge = dropReads(&readStart, incr);
				//printf("%d of %d\n", readStart, numquery);
				if(readStart >= numquery - 1)	//no reads left 
				{
					//maybe exit if the likelihoods are still very different...
					if(L2 - L1 < 1.0)
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
			if(L1 > L2)
			{
				twice++;
				oldNodePointer = nodePointer;
				oldReadStart = readStart;
				testAge = nextNodeAge;
				nextNodeAge = dropReads(&readStart, incr);
				if(readStart == numquery)	//no reads left 
				{
					//I think this is ok as it means the last set of reads had the right trend
					break;
					//printf("Error: Could not find max age for opt!\n");
					//exit(0);
				}

			}
		}while(twice < 2);	//DOUBLE CHECK THIS/Think of better way

		//Current stop gap... not great
		readStart = oldReadStart;
		nodePointer = oldNodePointer;

		printf("Maximum bound found, now finding optimum age\n");
		//printf("The elapsed time for rough estimation is %.16f seconds\n", ( ((double) clock()) / CLOCKS_PER_SEC) - time2);
	}
	else
	{
		printf("Not enough reads for bound optimization, using all\n");
	}

	//Some bounds or fillers added
	p[0] = readStart;
	p[1] = nodePointer;
	p[2] = 0;
	invector[1] = nextNodeAge/2;
	lowbound[1] = eh0;
	upbound[1] = nextNodeAge - eh0;
	nfun=0;
	onDindic=0;	// Because of 2nd layer of optimization

	//double time3 = (double) clock()/CLOCKS_PER_SEC;

	// Decide how the inputs may need to change at some point I guess
	est_age_lik = GoldenSection(invector,lowbound, upbound, 1, getlike_ages, p, 3);
	//printf("The elapsed time for age estimation is %.16f seconds\n", ( ((double) clock()) / CLOCKS_PER_SEC) - time3);

	est_age = invector[1];
	
	if(nextNodeAge - est_age < nextNodeAge/100.0)
	{
		printf("Warning! Age estimate was near boundary of %.16f\n", nextNodeAge);
	}

	printf("Opt found, calculating confidence intervals\n");

	//Use fisher information to get rough confidence interval
	//Maybe give option for bootstrap confidence interval
	//double time4 = (double) clock()/CLOCKS_PER_SEC;
	confI = confidenceIntervalFisher(par, est_age, est_age_lik, readStart);
	//printf("The elapsed time for confidenceIntervalFisher is %.16f seconds\n", ( ((double) clock()) / CLOCKS_PER_SEC) - time4);

	//Also maybe incorporate more options for the confidence interval, not just 95%
	//Would need to be able to calculate the z-score from the user given value
	//confI = 1.96 / sqrt(-secD);

	printf("Estimated age is %.16f with likelihood %.16f and 95%% confidence interval [%.16f,%.16f]\n", est_age, est_age_lik, est_age - confI, est_age + confI);

	//confidenceSearch(bounds, chiValue, maxAge, est_age, -est_age_lik, p);	
	//printf("Estimated age is %.16f with likelihood %.16f and %.2f%% confidence interval [%.16f,%.16f]\n", est_age, est_age_lik, chiValue, bounds[0], bounds[2]);
	//printf("Confidence interval is [%.16f, %.16f] with likelihoods %.16f and %.16f\n", bounds[0], bounds[2], bounds[1], bounds[3]);
	
	if(allTrees)
	{
		for(int i = 0; i < numTrees; i++)
		{
			//Some bounds or fillers added
			p[0] = readStart;
			p[1] = i;
			p[2] = 0;
			invector[1] = nextNodeAge/2;
			lowbound[1] = eh0;
			upbound[1] = nextNodeAge - eh0;
			nfun=0;
			onDindic=0;	// Because of 2nd layer of optimization

			//double time3 = (double) clock()/CLOCKS_PER_SEC;

			// Decide how the inputs may need to change at some point I guess
			est_age_lik = GoldenSection(invector,lowbound, upbound, 1, getlike_ages_tree, p, 3);
			//printf("The elapsed time for age estimation is %.16f seconds\n", ( ((double) clock()) / CLOCKS_PER_SEC) - time3);

			est_age = invector[1];

			printf("Estimated age is %.16f with likelihood %.16f for tree %d\n", est_age, est_age_lik, i);
		}
	}
}

//Something in here is causing a seg fault :(
//NEED TO REWRITE AT SOME POINT
void age_like_distribution_jointly_for_all2D(double **par)
{
	printf("Creating age likelihood distribution for all\n");
	int i, k, v, nfun, readStart = 0;
	int numDrops = 0;
	double p[3], L1, L2, ageIncr, nextNodeAge, tempAge, tempAge2;
	double invector[3],lowbound[3], upbound[3], eh0=1e-8, age_like;

	orderReads();

	//An important question here is why was the testAge for below so off?
	//nodePointer = assignments[usedReads[readStart]];
	//nextNodeAge = nodeages[nodeOrder[nodePointer]] + bls[nodeOrder[nodePointer]];
	nextNodeAge = nodeages[treeAssign[usedReads[i]]][assignments[usedReads[readStart]]] + bls[treeAssign[usedReads[i]]][assignments[usedReads[readStart]]];

	//printf("First node age as old %.16f and by assigning nodePointer %.16f\n")
	
	//printf("First node age %.16f\n", nextNodeAge);
	//printf("First node age %.16f from assignment %d with node age %.16f and branch length %.16f that has order %d\n", nextNodeAge, nodePointer, nodeages[nodeOrder[nodePointer]], bls[nodeOrder[nodePointer]], nodeOrder[nodePointer]);

	ageIncr = nextNodeAge/10000.0;		// make this an option later

	testAge = eh0;

	onDindic=1;

	// assignmentMode of 0 means single assignment given
	age_like = 0.0;
	while (testAge <= totMaxAge)
	{
		// Trimmed because you can't  actually compare the likelihoods, just the trends.	
		if (nextNodeAge <= testAge)
		{	
			// This is probably unneeded due to the while loop condition
			if(nextNodeAge >= totMaxAge)
			{
				break;
			}
			nextNodeAge = dropReads(&readStart, 0.01);	//Not sure if this is correct, may want to come back/make two dropReads versions
			ageIncr = nextNodeAge/10000.0;
			numDrops++;

			testAge = eh0;
		}
		printf("Testing age %.16f with dropped rounds %d with nextNodeAge %.16f and read start at %d:\n", testAge, numDrops, nextNodeAge, readStart);
		for (i=readStart; i<numquery; i++)
		{
			p[0] = usedReads[i];
			p[1] = treeAssign[usedReads[i]];
			p[2] = assignments[usedReads[i]];

			nfun=0;
			invector[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]]/2.0;
			lowbound[1] = eh0;
			upbound[1] = bls[treeAssign[usedReads[i]]][assignments[usedReads[i]]]-eh0;

			L2 =  GoldenSection(invector,lowbound, upbound, 1, getlike_gamma_root_in_trifurcation_testAge, p, 3); 

			//printf("\tRead %d assigned to %d with likelihood %lf with root placement of %.16f\n", usedReads[i], assignments[usedReads[i]], L2, nodeages[assignments[usedReads[i]]]+invector[1]);

			age_like += L2; //sum of log likelihoods
		}
		printf("\tLikelihood: %.16f\n", age_like);
		age_like = 0.0;
		testAge += ageIncr;
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


int main(int argc, char *argv[])
{
	//time for whole program:
	//double time1 = (double) clock()/CLOCKS_PER_SEC;

	double L;

	char assignfile[500], fraclikefile[500], querydatafile[500], referencedatafile[500], errorfile[500], strTree[100], tempFileName[500];
	int mode, allTrees;

	//To do: make a more user friendly command/flag interface
	if (argc != 10)
	{
		//printf("Specify name of five infiles: assignmentfile,fractionallikehoodfile, querydatafile, referencedatafile, and GTR+Gamma parameterfile, assingment mode, and a likelihood mode option\nMaximum name length: 30 characters\n");
		printf("Specify path of sample assignment file, sample alignment file, likelihood directory, parameter/tree directory, number of trees, assignment mode, errorprofile, output file, and 0/1 for each tree age estimate\n");
		exit(-1);
	}

	sprintf(assignfile, "%s", argv[1]);
	sprintf(querydatafile, "%s", argv[2]);
	sprintf(fraclikefile, "%s", argv[3]);	//likelihoods
	sprintf(referencedatafile, "%s", argv[4]);	//parameters and reference data
	//sprintf(GTRAparfile, "%s", argv[5]);
	const char *nptr = argv[5];
	mode=atoi(argv[6]);
	sprintf(errorfile, "%s", argv[7]);
	allTrees = atoi(argv[9]);
	//errorTest = atof(argv[10]);
	
	//printf("Mode: %d\n", mode);

	//New checks for reading in files, following guidelines from get_frac_like.c
	//Get number of trees and check its a valid number
    	//Reference: https://stackoverflow.com/questions/26080829/detecting-strtol-failure
    	
    	char *endptr = NULL;  
    	numTrees = 0;

    	//reset errno to 0 before call for validation of number
    	errno = 0;
    	numTrees = strtol (nptr, &endptr, 10);

    	/* test return to number and errno values,  */
    	if (nptr == endptr || (errno == ERANGE && numTrees == LONG_MIN) || (errno == ERANGE && numTrees == LONG_MAX) || errno == EINVAL || (errno != 0 && numTrees == 0) || (errno == 0 && nptr && *endptr != 0))
    	{
    	    printf ("Error in reading number of trees. Please double check. If you think this is an error, please contact maya_lemmon-kishi@berkeley.edu\n");
    	    exit(0);
    	} 
    	
    	//printf("Number of trees being read %lu\n", numTrees);

	numseq = (int*)malloc(sizeof(int) * numTrees);

	//Check if input directories exists
   	//sprintf(likeDir, "%s", argv[3]);
	struct stat s;
	int err = stat(fraclikefile, &s);
	if(-1 == err) {
	    if(ENOENT == errno) {
	        printf("Likelihood directory does not exist. Please double check input\n");
		exit(0);
	    } else {
	        perror("Error in stat");
	        exit(1);
	    }
	} else {
	    if(!S_ISDIR(s.st_mode)) {
	        printf("Likelihood path is not a directory. Please double check input\n");
		exit(0);
	    }
	}

	//sprintf(treeDir, "%s", argv[4]);
	err = stat(referencedatafile, &s);
	if(-1 == err) {
	    if(ENOENT == errno) {
	        printf("Tree directory does not exist. Please double check input\n");
		exit(0);
	    } else {
	        perror("Error in stat");
	        exit(1);
	    }
	} else {
	    if(!S_ISDIR(s.st_mode)) {
	        printf("Tree path is not a directory. Please double check input\n");
		exit(0);
	    }
	}
	

	
	read_data(assignfile, fraclikefile, querydatafile, referencedatafile, errorfile);

	doNRinits(2);
	inittransitionmatrix();
	//unrolled
	for(int i = 0; i < numTrees; i++)
	{
		if(usedTrees[i] == 0)
		{
			// no reads were assigned to this tree, skip
			continue;
		}

		pi[i][0] = log(pi[i][0]);
		pi[i][1] = log(pi[i][1]);
		pi[i][2] = log(pi[i][2]);
		pi[i][3] = log(pi[i][3]);
	}

	for(int treeNum = 0; treeNum < numTrees; treeNum++)
	{
		if(usedTrees[treeNum] == 0)
		{
			// no reads were assigned to this tree, skip
			continue;
		}
		//printf("Testing end of file read\n");
		sprintf(strTree, "%01d", treeNum);
		tempFileName[0] = '\0';
		strcat(tempFileName, referencedatafile);
		strcat(tempFileName, "/");
		strcat(tempFileName, strTree);
		strcat(tempFileName, "_reference.txt");

		if (NULL==(infile=fopen(tempFileName,"r")))
		{
			printf("Cannot open infile with tree data\n");
			exit(-1);
		}
		// Calculate the file size
		fseek(infile, 0, SEEK_END);
		long file_size = ftell(infile);

		// Move the file pointer to the beginning of the last line
		fseek(infile, 0, SEEK_SET);

		// Read the file character by character
		char ch;
		int line_length = 0;
		int last_line_length = 0;

		while ((ch = fgetc(infile)) != EOF) {
			line_length++;

			if (ch == '\n') {
				last_line_length = line_length;
				line_length = 0;
			}
		}

		fseek(infile, -last_line_length, SEEK_END);

		//Globals from Rasmus code, may want to fix...
		tip,comma=0;

		allocatetreememmory(numseq[treeNum]);
		int root = getclade(numseq[treeNum]) - 1 + numseq[treeNum];	//converts to ratePlacer node

		fclose(infile);

		////printtree(numseq, root);
		if(mode != 0 && mode != 6)
		{
			bestAssignment(root, treeNum);
		}

		// No need to keep in memory
		freetreememmory();
	}

	if (mode == 0 || mode == 5)
	{
		// 0 for no edge reassignment
		// 5 for edge reassignment
		maximize_like_seperately_for_all2D_Print(par);
	}
	else if (mode==1 || mode == 6)
		maximize_like_seperately_for_all2D(par);
	else if (mode==2)
		likelihoodratiotest_for_all(par);
	else if (mode == 3)
		maximize_like_jointly_for_all2D(par, allTrees);
	else if (mode == 4)
		printf("Likelihood distribution for all currently broken\n");
		//age_like_distribution_jointly_for_all2D(par);
	// Should we make a similar function but for each read??
	else
		printf("Please specify run mode.\n 1 for age of each read, 2 for LLR of each read, 3 for sample age estimation, and 4 for likelihood surface\n");

	freeNRinits(2);

	//printf("The elapsed time for whole program is %.16f seconds\n", ( ((double) clock()) / CLOCKS_PER_SEC) - time1);
}



