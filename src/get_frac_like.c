// Created by Rasmus Nielsen
// get_frac_like.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "tools.h"
#include <errno.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>

#define STATESPACE 4 //for the four nucleotides
#define NUMCAT 4/*number of categories in the discretization of the gamma for the nucleotide substituion model*/

FILE *infile, *infile2, *outfile;

double UFC;	// What is this?
// PMAT[i][j] is transition probability from i to j
double LRVEC[STATESPACE][STATESPACE], RRVEC[STATESPACE][STATESPACE], RRVAL[STATESPACE], PMAT1[STATESPACE][STATESPACE], PMAT2[STATESPACE][STATESPACE];
int tip,comma=0; /*globals used to read in the tree. Old code - don't ask.*/
double *statevector;
int DEBUG = 0;
int DEBUGGAMMA = 0;
int CHECKACCURACY = 0;
int CALCFULLLIKE = 0;


struct node {
    int up[2];	// for storing children node number, assuming binary. -1 means no children (current node is leaf)
    int down;	// for storing parent node number
    int *seq;	//sequence stored at the node
    double bl;	//Scaled branch length to parent
    double blstore;	//Branch length from newick file to parent
    double *like;	//size of 4, for each possible base, the conditional likelihood P(node = i|subtree)  
    double *condlike;	//size of 8, for each possible base, stores the fractional likelihood of the node (likelihood of tree - node's subree)
    double *posterior;	//size of 4, for each possible base, stores the fractional likelihood * V(bl_np)
    double SCALEFACTOR1;
    double SCALEFACTOR2;
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
 *             by Wilkinson and Reinsch, 1971
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
/*
complex conj (complex a)
{
    a.im = -a.im;
    return(a);
}*/

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
        }fprintf(outfile,"C%i\n",i+1);
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



void inittransitionmatrix(double pi[4])
{
    int i, j;
    double checksum, sum, RIVAL[STATESPACE], RIVEC[STATESPACE][STATESPACE],  A[STATESPACE][STATESPACE], workspace[2*STATESPACE];
    double par[6];

   /* for (i=0; i<4; i++){
         checksum=0.0;
        for (j=0; j<4; j++){
            fscanf(infile2,"%lf",&A[i][j]);
            checksum += A[i][j];
            }
        if (checksum != 0.0){
            printf("Rate matrix not properly scaled or not read correctly\n");
            printf("Checksum row %i: %lf\n",i,checksum);
            exit(-1);
            }
        }*/
    checksum=0.0;
    for (i=0; i<4; i++){
        fscanf(infile2,"%lf",&pi[i]);
        checksum+=pi[i];
        }
    if (checksum != 1.0){               // Should add the error term here maybe
        printf("WARNING: Nucleotide frequencies not properly scaled or not read correctly (checksum: %.15f)\n",checksum);
        //exit(-1);
	}
    for (i=0; i<6; i++)
        fscanf(infile2,"%lf",&par[i]);

    //printf("GTR Parameters: %lf, %lf, %lf, %lf, %lf, and %lf\n", par[0], par[1], par[2], par[3], par[4], par[5]);

    A[0][1]=pi[1]*par[0];
    A[0][2]=pi[2]*par[1];
    A[0][3]=pi[3]*par[2];
    A[1][0]=pi[0]*par[0];
    A[1][2]=pi[2]*par[3];
    A[1][3]=pi[3]*par[4];
    A[2][0]=pi[0]*par[1];
    A[2][1]=pi[1]*par[3];
    A[2][3]=pi[3]*par[5]; //unscaled rate of GT = 1.0
    A[3][0]=pi[0]*par[2];
    A[3][1]=pi[1]*par[4];
    A[3][2]=pi[2]*par[5]; //unscaled rate of GT = 1.0

    //printf("Q matrix diagonals: ");
    for (i=0; i<4; i++)
        {
        A[i][i]=0.0;
        sum=0.0;
        for (j=0; j<4; j++)
            sum = sum + A[i][j];
        A[i][i] = -sum;
	//printf("%lf, ", A[i][i]);
        }
    //printf("\n");


    if (eigen(1, A[0], STATESPACE, RRVAL, RIVAL, RRVEC[0], RIVEC[0], workspace) != 0)
        {
        printf("Transitions matrix did not converge or contained non-real values!\n");
        exit(-1);
        }
    for (i=0; i<STATESPACE; i++)
        for (j=0; j<STATESPACE; j++)
            LRVEC[i][j] = RRVEC[i][j];
    if (matinv(RRVEC[0],STATESPACE, STATESPACE, workspace) != 0)
        printf("Could not invert matrix!\nResults may not be reliable!\n");
}



void maketransitionmatrix(int matnum, double t)
{
    int i, j, k;
    double EXPOS[STATESPACE];
    
    for (k=0; k<STATESPACE; k++)
        EXPOS[k] = exp(t*RRVAL[k]);
    for (i=0; i<STATESPACE; i++)
    {
        for (j=0; j<STATESPACE; j++)
        {
            if (matnum==0) PMAT1[i][j] = 0.0;
            else PMAT2[i][j] = 0.0;
            for (k=0; k<STATESPACE; k++)
            	{
                if (matnum==0) PMAT1[i][j] =  PMAT1[i][j] + RRVEC[k][j]*LRVEC[i][k]*EXPOS[k];
                else PMAT2[i][j] =  PMAT2[i][j] + RRVEC[k][j]*LRVEC[i][k]*EXPOS[k];
        	}
    	}
	}
/*	   if (matnum==1){
printf("PMAT2 (t=%lf)\n",t);
for (i=0; i<STATESPACE; i++){
     for (j=0; j<STATESPACE; j++)
     printf("%lf ",PMAT2[i][j]);
    printf("\n");} 
}*/
}

/*this calculates conditional likelihoods assuming leaf nodes have already been initialized */
void makecon(int node, int site)

{
    int i, j, D1, D2;
    double L, max;
   
    // Get node number for each child
    D1 = tree[node].up[0];
    D2 = tree[node].up[1];

    max=0.0;

    // Recursion to children nodes since we need likelihood information from each
    if (tree[D1].up[0] > -1) 
    	makecon(D1, site);
    if (tree[D2].up[0] > -1) 
    	makecon(D2, site);

    // Generate transition matrix based on branch length to child
    maketransitionmatrix(0, tree[D1].bl);
    maketransitionmatrix(1, tree[D2].bl);

    // For each possible nucleotide
    for (i=0; i<STATESPACE; i++){
	// For child 1, store likelihood in L
    	if (tree[D1].up[0] > -1){
	    // When there is grandchildren, sum over child's conditional likelihoods
	    L=0;
	    for (j=0; j<STATESPACE; j++)
           	L += PMAT1[i][j] * tree[D1].like[j];
	}
	else if (tree[D1].seq[site] > -1)
	    // When there is no grandchildren (ie children are leaf) and have sequence information
	    L = PMAT1[i][tree[D1].seq[site]];
        else
	    // When there is no grandchildren (ie children are leaf) and base is unknown (n, -, ~)
	    L = 1.0;

	// For child 2, store likelihood in tree[node].like[i]
    	if (tree[D2].up[0] > -1){
	    // When there is grandchildren, sum over child's conditional likelihoods
	    tree[node].like[i] = 0.0;
	    for (j=0; j<STATESPACE; j++)
		tree[node].like[i] += PMAT2[i][j] * tree[D2].like[j];
    	}
    	else if (tree[D2].seq[site] > -1)
	    // When there is no grandchildren (ie children are leaf) and have sequence information
	    tree[node].like[i] = PMAT2[i][tree[D2].seq[site]];
        else
	    // When there is no grandchildren (ie children are leaf) and base is unknown (n, -, ~)
	    tree[node].like[i]  = 1.0;

	// Gets conditional likelihood of i by multipling likelihoods of each subtree
    	tree[node].like[i] = tree[node].like[i] * L;

	// For scaling purposes
	if (tree[node].like[i]  > max)
	    max = tree[node].like[i];
    }

    // Scale conditional
    for (i=0; i<STATESPACE; i++)
        tree[node].like[i]=tree[node].like[i]/max;

    max = log(max);

    // Calculate scaling factor and save
    tree[node].SCALEFACTOR1 = max + tree[D1].SCALEFACTOR1 + tree[D2].SCALEFACTOR1;
    UFC = UFC + max;
}

//this function calculates the likelihood and stores the fractional likelihoods up to a scaling factor in each node
// ?? stores??
double tlike(int root, int site, double pi[4])

{ /*assumes uniform prior distribution on the state space*/ 
    int i, j, s;
    double Lsum=0.0;
    
    UFC=0.0;
    if (tree[root].up[0]==-1) {
	printf("Data contain only one species\n");
	exit(-1);
    }	

    // Caculates the conditional likelihood of each node above root
    makecon(root, site);

    // Sums conditional likelihood to get likelihood of tree
    for (i=0; i<STATESPACE; i++)
        Lsum += pi[i]*tree[root].like[i];
  /*  printf("Root likelihoods:");
    for (i=0; i<STATESPACE; i++)
        printf("%lf ",tree[root].like[i]);
    printf("\n");*/
    return log(Lsum)+UFC; 
}


void scale_branch_lengths(int numleaves, double lambda)
{
    int i;
    for (i=0; i<2*numleaves-1; i++)
        tree[i].bl=tree[i].blstore*lambda;
}

// Root is of given node, I don't think of the actual root?
void calculate_full_likelihood(int root, int numbase, int numleaves, double pi[4])
{ /*assumes uniform prior distribution on the state space*/
    
    int i, j, v;
    double L, CL, Lsum=0.0;
    
    if (tree[root].up[0]==-1) {
	    printf("Data contain only one species or have other error\n");
	    exit(-1);
    }

    for (i=0; i<numbase; i++){
        CL=0.0;
        for (j=0; j<NUMCAT; j++){
            UFC=0.0;
            scale_branch_lengths(numleaves, statevector[j]); //this is in general a slow way of doing this
            makecon(root, i);
            L=0;
            for (v=0; v<STATESPACE; v++)
                L += pi[v]*tree[root].like[v];
            CL += L*exp(UFC);//then the underflow protection doesn't work - so this is only inteneded when there is a small number of species
	    if (DEBUG)
		printf("Unscaled likelihood site %i, cat %i: %.15f\n",i,j,log(L*exp(UFC)));
	}
 	if (DEBUG)
	    printf("Likelihood site %i: %.15f\n",i,log(CL/NUMCAT));
	Lsum += log(CL/NUMCAT);
    }
    printf("Full Like: %.15f\n",Lsum);
}

//this algorithm runs from the root to calculate the conditional likelihood given all data towards the parental node for each node
void makeposterior(int node, int site)
{
    int i,j, parent, otherb;
    double max, bl, sum, templike[STATESPACE];
    
    max = 0.0;
    parent = tree[node].down;
    bl = tree[node].bl;

    maketransitionmatrix(0, bl);

    // Get the node number of sibling
    if ((otherb = tree[parent].up[0])==node)
	otherb = tree[parent].up[1];

    maketransitionmatrix(1, tree[otherb].bl);

    for (i=0; i<STATESPACE; i++){
        templike[i]=0;

        if (tree[otherb].up[0] > -1){
	    // If sibling is not leaf
	    for (j=0; j<STATESPACE; j++)
		templike[i]=templike[i]+tree[otherb].like[j]*PMAT2[i][j];
        }
        else if (tree[otherb].seq[site] > -1)
	    // If sibling is leaf and base is known (a,c,g,t)
            templike[i]=PMAT2[i][tree[otherb].seq[site]];//edited here as well
        else
	    // When sibling is leaf and base is unknown (n, -, ~)
	    templike[i] = 1.0;

	// Multiply sibling's subtree * parent's posterior (ie rest of tree to parent's node)
        templike[i]=templike[i]*tree[parent].posterior[i];

        //I think this is it
	// Store conditional likelihood (likelihood of tree - node's subtree)
        tree[node].condlike[i] = log(templike[i])+tree[otherb].SCALEFACTOR1+tree[parent].SCALEFACTOR2;
    }

    if (tree[otherb].up[0]>-1)
	tree[node].SCALEFACTOR2 = tree[otherb].SCALEFACTOR1;
    else
	tree[node].SCALEFACTOR2=0.0;

    // Calculate posterior of node by taking conditional likelihood and * with branch length to parent
    for (i=0; i<STATESPACE; i++){
        tree[node].posterior[i]=0.0;

        for (j=0; j<STATESPACE; j++){
            tree[node].posterior[i] = tree[node].posterior[i] + PMAT1[i][j]*templike[j];
        }

        if (tree[node].posterior[i]  > max)
            max = tree[node].posterior[i];
    }

    for (i=0; i<STATESPACE; i++)
        tree[node].posterior[i]=tree[node].posterior[i]/max; 

    tree[node].SCALEFACTOR2 += log(max)+tree[parent].SCALEFACTOR2;

    // If there are children
    if (tree[node].up[0] > -1)
    {
	//Recurse this function to children
    	makeposterior(tree[node].up[0], site);
    	makeposterior(tree[node].up[1], site);
    }
}

//finds the fractional likelihoods for all nodes 
double get_likes(int numleaves, int root, int site, double pi[4])
{
    int i, j;
    double L;

    L = tlike(root, site, pi); /*need to call likelihood first*/ 

    if (DEBUG) printf("Likelihood for site %i: %.15f\n",site,L);

    for (i=0; i<STATESPACE; i++)
        tree[root].posterior[i]=1.0;

    tree[root].SCALEFACTOR2=0.0;
    tree[root].SCALEFACTOR1=UFC;

    for (i=numleaves-1; i<2*numleaves-1;i++)
        tree[i].SCALEFACTOR1=0.0;

    makeposterior(tree[root].up[0], site);
    makeposterior(tree[root].up[1], site);

    // Conditional + fractional likelihoods for internal
    for (j=0; j<numleaves-1; j++){
        for (i=0; i<4; i++){
            tree[j].condlike[i+4] = log(tree[j].like[i])+tree[j].SCALEFACTOR1;
           // printf("node %i, site: %i, base %i: %lf\n",j, site,i,tree[j].like[i+4]);
        }
    }

    // Conditional + fractional likelihoods for leafs
    for (j=numleaves-1; j<2*numleaves-1; j++){
        for (i=0; i<4; i++){
            if (tree[j].seq[site]==i)
		tree[j].condlike[i+4] = 0.0;
            else
		tree[j].condlike[i+4] = -100000000.0;
        }
    }

    return L;
}




void allocatetreememmory(int numleaves)
{
    int i;
   
    // Number of nodes based off full binary tree
    tree=malloc((numleaves*2-1)*(sizeof(struct node)));

    for (i=0; i<(numleaves*2-1); i++){
            tree[i].condlike = malloc(8*(sizeof(double)));
            tree[i].like = malloc(4*(sizeof(double)));
            tree[i].posterior = malloc(4*(sizeof(double)));
    }
}


/*some old code for reading a fasta file and storing DNA as ints*/
/** TO DO: Need to fix so don't need to awk input file **/
int readseq(int *numbase)
{
    int i, j, numleaves, test, k=0;
    char c;
   
    // Read in first line that contains number of reference species and length of reference sequences
    fscanf(infile,"%i %i",&numleaves,numbase);
    //printf("there are %i species and %i bases\n",numleaves,*numbase);

    if (numleaves < 3)
    {
	    printf("This is for more than two seq.s only!\n");
	    exit(-1);
    }

    // Allocation of tree based on number of leaf nodes
    allocatetreememmory(numleaves);

    for (i=numleaves-1; i<2*numleaves-1; i++)
        tree[i].seq=malloc((*numbase)*(sizeof(int)));
    while ((c=(fgetc(infile)))!='\n');
    do
    {
        for (i=0; i<numleaves; i++)
        {
            j=k;
            while ((c=tolower((fgetc(infile))))!='\n')
            {
                if ((c!=' ')&&(c!='\t')) {
                    if (c=='a') tree[i+numleaves-1].seq[j] = 0;
                    else if (c=='c') tree[i+numleaves-1].seq[j]= 1;
                    else if (c=='g') tree[i+numleaves-1].seq[j] = 2;
                    else if (c=='t') tree[i+numleaves-1].seq[j] = 3;
                    else if (c=='n'||c=='-'||c=='~') tree[i+numleaves-1].seq[j] = -1;
                    else{
                        printf("\nBAD BASE (%c) in species %i base %i",c,i+1,j+1);
                        scanf("%i",&i);
                        exit(-1);
                    }
                    j++;
                    if (i==(numleaves-1)) k++;
                }
            }
        }
    } while (k<*numbase);
    return numleaves;
}

void freetreememmory(int numleaves)
{
    int i;
    
    for (i=0; i<(numleaves*2-1); i++){
        	free(tree[i].like);
            free(tree[i].condlike);
            free(tree[i].posterior);
        }
    free(tree);
}

void print_likes(int site, int numleaves)
{
    int i, j;

    for (j=0; j<2*numleaves-1; j++){
	printf("Node %i:",j);
        for (i=0; i<4; i++)
	     printf(" %.15f",tree[j].condlike[i]);
        printf (" and ");
        for (i=4; i<8; i++)
             printf(" %.15f",tree[j].condlike[i]);
        printf("\n");
    }
}

void printdata(int numleaves,int numbase)
{
    int i, j;

    for (i=0; i<numleaves;i++){
        printf("seq %i: ",i);
        for (j=0; j<numbase; j++)
            printf("%i",tree[i+numleaves-1].seq[j]);
        printf("\n");
    }
}

void check_likecalc(int site, int numleaves, int root, double pi[4])
{
    int i, j, k;
    double L, loc;

    for (i=0; i<2*numleaves-1;i++)
    {
    	L=0.0;
    	if (i!=root){
    	    maketransitionmatrix(0, tree[i].bl);
    	    for (j=0; j<4; j++){
    	        loc = 0.0;
    	        for (k=0; k<4; k++)
    	            loc += PMAT1[j][k]*exp(tree[i].condlike[k]);
    	        L += loc*exp(tree[i].condlike[j+4])*pi[j];
    	        }
    	    }
    	else {
    	    for (j=0; j<4; j++)
    	        L += exp(tree[i].condlike[j+4])*pi[j];
    	    }
    	printf("Likelihood calculated for node %i, site %i: %.15f\n",i,site,log(L));
    }
}

void checkaccuracy(int numleaves, int numbase, int root, double pi[4])
{
    int i, j, k, site;
    double L, P, loc, e, TL=0.0;
    
    printf("Checking accuracy of calculations:\n");
    for (site=0; site<numbase; site++)
    {
    	P = get_likes(numleaves, root, site, pi);
    	e = 0.0;
    	for (i=0; i<2*numleaves-1;i++)
    	    {
    	    L=0.0;
    	    if (i!=root){
    	        maketransitionmatrix(0, tree[i].bl);
    	        for (j=0; j<4; j++){
    	            loc = 0.0;
    	            for (k=0; k<4; k++)
    	                loc += PMAT1[j][k]*exp(tree[i].condlike[k]);
    	            L += loc*exp(tree[i].condlike[j+4])*pi[j]; 
    	            }
    	        }
    	    else {
    	        for (j=0; j<4; j++)
    	            L += exp(tree[i].condlike[j+4])*pi[j];   
    	        }
    	    e += (P-log(L))*(P-log(L));
    	    }
    	printf("Total squared numerical error for site %i: %.15f (L: %.15f)\n",site,e,P);
    	TL += P;
    }
}

//Yang 1994 median discretization with standardization
//  New gamma discretization
void definegammaquantiles(int k, double alpha, double beta)
{
    int i;
    double sum=0.0;;

    for (i=0; i<k; i++){
        statevector[i] = PointGamma(((double)i+0.5)/(double)k, alpha, beta);
        sum+=statevector[i];
    }
    for (i=0; i<k; i++)
        statevector[i]=(double)k*statevector[i]/sum;
    if (DEBUGGAMMA) {
        printf("Rates for discretized Gamma distribution: ");
        for (i=0; i<k; i++)
            printf("%.15f ",statevector[i]);
        printf("\n");
    }

}

//void definegammaquantiles(int k, double alpha, double beta)
//{
//    int i;
//    double mean=0.0;;
//    
//    for (i=0; i<k; i++){
//        statevector[i] = PointGamma(((double)i+0.5)/(double)k, alpha, beta);
//	    if (DEBUGGAMMA)  mean+=statevector[i];
//    }
//    if (DEBUGGAMMA) printf("mean: %lf\n",mean/NUMCAT);
//}

void store_branch_lengths(int numleaves)
{
    int i;
    for (i=0; i<2*numleaves-1; i++)
        tree[i].blstore=tree[i].bl;
}

void printffraclikelihoods(int numleaves)
{
    int i, j;

    // leaf nodes first
    for (i=numleaves-1; i<2*numleaves-1; i++){
	//fprintf(outfile, "\n%d or %d ", i, i - numleaves + 1);
        for (j=0; j<4; j++)
            fprintf(outfile,"%.15f ",tree[i].condlike[j]);
    }
    // internal nodes
    for (i=0; i<numleaves-1; i++){
	//fprintf(outfile, "\n%d or %d ", i, i + numleaves);
        for (j=0; j<8; j++)
            fprintf(outfile,"%.15f ",tree[i].condlike[j]);
    }
}

double find_ages(int node, int numleaves, int *nodeOrder, double *nodeAge, double *nodebl)
{
    int u1, u2;
    double t, t2, d;
    
    if ((u1=tree[node].up[0])>-1){
        u2=tree[node].up[1];
        t = find_ages(u1,numleaves, nodeOrder, nodeAge, nodebl) + tree[u1].blstore;
        t2 = find_ages(u2,numleaves, nodeOrder, nodeAge, nodebl) + tree[u2].blstore;
        d = t - t2;
        if (d<-0.000001 || d>0.00001)
        {
            printf("Tree not ultrametric");
            exit(-1);    
        }
        t = (t+t2)/2.0;
	nodeOrder[node + numleaves] = node + numleaves;
	nodeAge[node + numleaves] = t;
	nodebl[node + numleaves] = tree[node].bl;
	//fprintf(outfile,"%i: %lf %lf\n",node+numleaves,t,tree[node].bl);
        return t;
    }
    else {
        //fprintf(outfile,"%i: %lf %lf\n",node-numleaves+1,0.0,tree[node].bl);
	nodeOrder[node - numleaves + 1] = node - numleaves + 1;
	nodeAge[node - numleaves + 1] = 0.0;
	nodebl[node - numleaves + 1] = tree[node].bl;
        return 0.0;
    } 
}

//merge functionality, using in-place merge sort to reduce memory load
// Will base merge on nodeAge[i]+nodebl[i] to get age in which reads assigned to i will need to cut in ratePlacer
// for reference: https://www.geeksforgeeks.org/merge-sort/ and https://www.geeksforgeeks.org/in-place-merge-sort/
void merge(int *nodeOrder, double *nodeAge, double *nodebl, int left, int right, int mid)
{
	// If already merged, are there cases where this may be true but not sorted?
	if (nodeAge[mid] + nodebl[mid] <= nodeAge[mid+1] + nodebl[mid+1])
	{
		return;
	}

	int midstart = mid + 1;

	while(left <= mid && midstart <= right)
	{
		// if elements are in right place (ie left < right)
		if(nodeAge[left] + nodebl[left] <= nodeAge[midstart] + nodebl[midstart])
		{
			left++;
		}
		else
		{
			int index = midstart, tempOder = nodeOrder[midstart];
			double tempAge = nodeAge[midstart], tempbl = nodebl[midstart];

			// shift all elements between left and midstart to right by 1
			while (index != left)
			{
				nodeOrder[index] = nodeOrder[index - 1];
				nodeAge[index] = nodeAge[index - 1];
				nodebl[index] = nodebl[index - 1];
				index--;
			}
			nodeOrder[left] = tempOder;
			nodeAge[left] = tempAge;
			nodebl[left] = tempbl;

			// Update
			left++;
			mid++;
			midstart++;	// may be redundant, but removes some extra calc?
		}
	}

}

// merge sort implementation for sort_ages
void mergeSort(int *nodeOrder, double *nodeAge, double *nodebl, int left, int right)
{
	if (left < right)
	{
		// Mid point calculation to hopefully avoid overflow
		int mid = left + (right - left) / 2;

		// recurse down
		mergeSort(nodeOrder, nodeAge, nodebl, left, mid);
		mergeSort(nodeOrder, nodeAge, nodebl, mid + 1, right);

		// merge
		merge(nodeOrder, nodeAge, nodebl, left, right, mid);
	}
}

//Sorting of nodes using ?? sort
void sort_ages(int numleaves, int *nodeOrder, double *nodeAge, double *nodebl)
{
	//This will sort the nodes of the tree from youngest to oldest by nodeAge + nodebl
	//	ie the sorting will be by node parent's age for ratePlacer to know which reads to cut
	
	mergeSort(nodeOrder, nodeAge, nodebl, 0, 2*numleaves-2);

	//print
	// maybe at some point add another quick sort based on node age when nodeage + nodebl are equal
	for(int i = 0; i < 2 * numleaves - 1; i++)
	{
		fprintf(outfile,"%i: %.15f %.15f\n", nodeOrder[i], nodeAge[i], nodebl[i]);
	}
}

// gfl inputDir outputDir numTrees
int main(int argc, char *argv[])
{
    int i, j, numleaves, numbase, root;
    // pi is the prior probability
    double alpha, pi[4], maxAge;
    char referencedatafile[500], GTRAparfile[500], inputDir[500], outputDir[500], strTree[100], outfilePath[500];
	int *nodeOrder;
	double *nodeAge, *nodebl;
   
    /** Check if number of files is correct and any flags **/
    if (argc != 4) {
	//printf("Specify name of two infiles and one outfile: referencedatafile, GTR+Gamma parameterfile, and likelihood outfile\nMaximum path name 500 characters\n");
	printf("Specify name of input and output directory and number of trees: input directory, output directory, number of trees\nMaximum path name 500 characters\n");
	exit(-1);
    }

    //Get number of trees and check its a valid number
    //Reference: https://stackoverflow.com/questions/26080829/detecting-strtol-failure
    const char *nptr = argv[3];
    char *endptr = NULL;  
    long numTrees = 0;

    //reset errno to 0 before call for validation of number
    errno = 0;
    numTrees = strtol (nptr, &endptr, 10);

    /* test return to number and errno values,  */
    if (nptr == endptr || (errno == ERANGE && numTrees == LONG_MIN) || (errno == ERANGE && numTrees == LONG_MAX) || errno == EINVAL || (errno != 0 && numTrees == 0) || (errno == 0 && nptr && *endptr != 0))
    {
        printf ("Error in reading number of trees. Please double check. If you think this is an error, please contact maya_lemmon-kishi@berkeley.edu\n");
	exit(0);
    } 
    
    printf("Number of trees being read %lu\n", numTrees);

    //Check if input directory exists
    sprintf(inputDir, "%s", argv[1]);
	struct stat s;
	int err = stat(inputDir, &s);
	if(-1 == err) {
	    if(ENOENT == errno) {
	        printf("Input directory does not exist. Please double check input\n");
		exit(0);
	    } else {
	        perror("Error in stat");
	        exit(1);
	    }
	} else {
	    if(!S_ISDIR(s.st_mode)) {
	        printf("Input directory is not a directory. Please double check input\n");
		exit(0);
	    }
	}

	//check if output directory exists, else create one
	sprintf(outputDir, "%s", argv[2]);
	err = stat(outputDir, &s);
	if(-1 == err) {
	    if(ENOENT == errno) {
	        printf("Creating output directory\n");
		mkdir(outputDir, 0700);
	    } else {
	        perror("Error in stat");
	        exit(1);
	    }
	}
	if(!S_ISDIR(s.st_mode)) {
		printf("Output directory is not a directory. Please double check\n");
		exit(0);
	}

	//Loop over number of trees next
	for(int treeNum = 0; treeNum < numTrees; treeNum++)
	{
		sprintf(strTree, "%d", treeNum);
		referencedatafile[0] = '\0';
		strcat(referencedatafile, inputDir);
		strcat(referencedatafile, "/");
		strcat(referencedatafile, strTree);
		strcat(referencedatafile, "_reference.txt");
		printf("Opening reference file %s\n", referencedatafile);

		/*READS DATA, TREE, AND PARAMETERS AND INITIALIZES ALL VARIABLES*/
    		if (NULL==(infile=fopen(referencedatafile,"r"))){
    		    puts("Cannot open reference data file.\n");
    		    exit(-1);
		}

		//Globals from Rasmus code, may want to fix...
		tip,comma=0;

    		numleaves = readseq(&numbase); //Reads in the sequence data and allocates memory
    		root = getclade(numleaves) - 1;  //Reads in the Newick tree
    		    
    		//printtree(numleaves, root);

    		fclose(infile);

		GTRAparfile[0] = '\0';
		strcat(GTRAparfile, inputDir);
		strcat(GTRAparfile, "/");
		strcat(GTRAparfile, strTree);
		strcat(GTRAparfile, "_parameter.txt");
		printf("Opening parameter file %s\n", GTRAparfile);

    		if (NULL==(infile2=fopen(GTRAparfile,"r"))){
    		    puts("Cannot open parameter file.\n");
    		    exit(-1);
    		}
    		//else printf("Opened parameter file successfully\n");
    		store_branch_lengths(numleaves);
    		fscanf(infile2,"%lf",&alpha);
    		inittransitionmatrix(pi); //THIS FUNCTION ALSO DOES THE MATRIX DIAGONALIZATION
    		fclose(infile2);
    		
    		/*DISCRETIZES THE GAMMA DISTRIBUTION*/
    		statevector = malloc(NUMCAT*(sizeof(double)));
    		definegammaquantiles(NUMCAT, alpha, alpha);

    		//OPENS THE INFILE FOR PRINITNG THE FRACTIONAL LIKELIHOODS
		outfilePath[0] = '\0';
		strcat(outfilePath, outputDir);
		strcat(outfilePath, "/");
		strcat(outfilePath, strTree);
		strcat(outfilePath, "_likelihood.txt");
		printf("Opening likelihood file %s\n", outfilePath);
    		if (NULL==(outfile=fopen(outfilePath,"w"))){
    		    puts("Hmm, cannot open outfile for likelihoods...");
    		    exit(-1);}
    		fprintf(outfile,"%i %i %i\n",numleaves,numbase,NUMCAT);
    		for (i=0; i<NUMCAT; i++)
    		    fprintf(outfile,"%.15f ",statevector[i]);
    		fprintf(outfile,"\n");

		//int nodeOrder[2*numleaves-1];
    		//double nodebl[2*numleaves-1], nodeAge[2*numleaves-1];

    		nodeOrder = (int*)malloc(sizeof(int) * (2*numleaves-1));
    		nodebl = (double*)malloc(sizeof(double) * (2*numleaves-1));
		nodeAge = (double*)malloc(sizeof(double) * (2*numleaves-1));
    		find_ages(root, numleaves, nodeOrder, nodeAge, nodebl);//PRINT OUT NODE AGES - MODIFY THIS IF CHANGING NODE LABELING

    		maxAge = nodeAge[root + numleaves];

    		//maxAge = getMaxAge(root); //go from root to leaves --- may not actually need this function if we already gather it
    		// Thought: Instead of printing within find_ages, we sort and then print to keep likelihood file clean  
    		sort_ages(numleaves, nodeOrder, nodeAge, nodebl);

    		fprintf(outfile, "%.15f\n", maxAge);
   
    		//printf("Starting calculations\n");
    		//LOOPS OVER ALL CATEGORIES OF THE DISCRETIZED GAMMA DISTRIBUTION AND OVER ALL SITES
    		for (i=0; i<NUMCAT; i++) //LOOPS OVER CATEGORIES
    		{
    		    if (DEBUG) printf("Calculating fractional likelihoods for category %i with rate %.15f\n",i+1,statevector[i]);
    		    fprintf(outfile,"C%i\n",i+1);
    		    scale_branch_lengths(numleaves, statevector[i]); //THIS FUNCTION IS RESPONSIBLE FOR SCALING THE BRANCHLENGTHS ACCORIDNG TO THE GAMMA DISCRETIZATION
    		    for (j=0; j<numbase; j++){//LOOPS OVER SITES
    		        fprintf(outfile,"S%i: ",j+1);
    		        get_likes(numleaves, root, j, pi);//CALCULATES THE FRACTIONAL LIKELIHOODS FOR EACH NODE
    		        if (DEBUG) check_likecalc(j, numleaves, root, pi);
    		        if (DEBUG) printf("Calculated likes for site %i\n",j);
    		        if (DEBUG) print_likes(j, numleaves); //prints fractional likelihoods. 
    		        printffraclikelihoods(numleaves);//PRINTS THE FRACTIONAL LIKELIHOODS
    		        fprintf(outfile,"\n");
    		    }
    		    if (CHECKACCURACY) checkaccuracy(numleaves, numbase, root, pi);//THIS FUNCTION CAN BE USED TO CHECK FOR NUMERICAL PRECISION
    		}
    		//if (CALCFULLLIKE) calculate_full_likelihood(root, numbase, numleaves, pi);

    		//printtree(numleaves, root);

    		fclose(outfile);
    		free(statevector);
    		freetreememmory(numleaves);
		free(nodeOrder);
		free(nodeAge);
		free(nodebl);
	}
}
   ////if you are reading in many trees with the same number of leafs, make sure no to repeatedly allocate and freeing the memory
//}
