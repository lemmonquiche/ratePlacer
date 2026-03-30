#include "jph.h"
#include "minfunc.h"

#define NR_END 1
#define FREE
#define TOLX2 1.0e-7
#define FREE_ARG char*
#define ALF 1.0e-4
#define STPMX 200.0
#define ITMAX 500
#define TOLX (4*EPS)
#define EPS 3.0e-8
#define NMAX 5000
#define INF DBL_MAX

#    define TOLER_PASS_1        0.001
#    define TOLER_PASS_2        praxisTol
#    define MAX_STEP_SIZE_1    1.0
#    define MAX_STEP_SIZE_2    1.0

double praxisTol = 0.00000001;

double *dg,*g,*hdg,*pnew,*xi,**hessin; 
int npar, CENTRALMODE;
FILE *tempfile;

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
		(maxarg1) : (maxarg2))

			     static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

void nrerror(char error_text[])
	/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

double *dvector(long nl, long nh)
	/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

void free_dvector(double *v, long nl)
	/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));

}


double **dmatrix(long nrl, long nrh, long ncl, long nch)
	/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


void free_dmatrix(double **m, long nrl, long ncl)
	/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}




int lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
		double *f, double stpmax, int *check, double (*func)(double [], double []), double lowbound[], double upbound[], double otherstuff[])
{
	int i;
	double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
	       test,tmplam;
	double new1;

	/*printf ("B: "); for (i=1;i<=n;i++) printf("%f ",p[i]);
	  printf("\n");*/
	*check=0;
	for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax)
		for (i=1;i<=n;i++) p[i] *= stpmax/sum;
	for (slope=0.0,i=1;i<=n;i++)
		slope += g[i]*p[i];
	test=0.0;
	for (i=1;i<=n;i++) {
		temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
		if (temp > test) test=temp;
	}
	if (test==0.0)
		return 1;
	alamin=TOLX2/test;
	alam=1.0;
	for (;;) {
		for (i=1;i<=n;i++)
		{
			new1=xold[i]+alam*p[i];
			if (new1 < lowbound[i])	/*my bloody code*/
				new1 = lowbound[i];
			if  (new1 > upbound[i])
				new1 = upbound[i];
			x[i] = new1;
			/*printf("x[i]: %f,newin[i]: %f,xold[i]: %f,alam: %f,p[i]: %f\n",x[i],newin[i],xold[i],alam,p[i]);*/
		}
		*f=(*func)(x,otherstuff);
		if (alam < alamin) {
			for (i=1;i<=n;i++) x[i]=xold[i];
			*check=1;
			return 10;
		} else if (*f <= fold+ALF*alam*slope) return 10;
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(*f-fold-slope));
			else {
				rhs1 = *f-fold-alam*slope;
				rhs2=f2-fold2-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else {
					disc=b*b-3.0*a*slope;
					if (disc<0.0) {/*nrerror("Roundoff problem in lnsrch.") return -1;*/disc=0.0;}
					/*else*/tmplam=(-b+sqrt(disc))/(3.0*a);
				}
				if (tmplam>0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = *f;
		fold2=fold;
		alam=FMAX(tmplam,0.1*alam);
	}
}

void doNRinits(int n)

{
	dg=dvector(1,n);
	g=dvector(1,n);
	hdg=dvector(1,n);
	hessin=dmatrix(1,n,1,n);
	pnew=dvector(1,n);
	xi=dvector(1,n);
}

void freeNRinits(int n)

{
	free_dvector(dg, 1);
	free_dvector(g,1);
	free_dvector(hdg,1);
	free_dmatrix(hessin,1,1);
	free_dvector(pnew,1);
	free_dvector(xi,1);
}

void dfpmin(double p[], int n, double gtol, int *iter, double *fret,
		double(*func)(double [], double[]), void (*dfunc)(double [], double [],double [], double [], double(*fu)(double [], double[]), double[]), double lowbound[], double upbound[], double otherstuff[])
{
	int lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
			double *f, double stpmax, int *check, double (*func)(double [], double[]), double lowbound[], double upbound[], double otherstuff[]);
	int check,i,its,j;
	double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;

	fp=(*func)(p,otherstuff);
	(*dfunc)(p,g, lowbound, upbound, func, otherstuff);
	for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++) hessin[i][j]=0.0;
		hessin[i][i]=1.0;
		xi[i] = -g[i];
		sum += p[i]*p[i];
	}
	stpmax=STPMX*FMAX(sqrt(sum),(double)n);
	for (its=1;its<=ITMAX;its++) {
		//printf("Its: %d\n", its);
		*iter=its;
		/*printf ("A: "); for (i=1;i<=n;i++) printf("%f ",g[i]); printf("\n");*/
		if (lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,&check,func,lowbound,upbound,otherstuff) == -1) /*MY CODE*/
			return;		fp = *fret;
		for (i=1;i<=n;i++) {
			xi[i]=pnew[i]-p[i];
			p[i]=pnew[i];
		}
		test=0.0;
		for (i=1;i<=n;i++) {
			temp=fabs(xi[i])/FMAX(fabs(p[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX) {
			/*FREEALL*/
			return;
		}
		for (i=1;i<=n;i++) dg[i]=g[i];
		(*dfunc)(p,g,lowbound,upbound, func,otherstuff);
		/*	printf ("C: "); for (i=1;i<=n;i++) printf("%f ",g[i]); printf("\n");*/
		test=0.0;
		den=FMAX(*fret,1.0);
		for (i=1;i<=n;i++) {
			temp=fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
			if (temp > test) test=temp;
		}
		if (test < gtol) {
			/*FREEALL*/
			return;
		}
		for (i=1;i<=n;i++) dg[i]=g[i]-dg[i];
		for (i=1;i<=n;i++) {
			hdg[i]=0.0;
			for (j=1;j<=n;j++) hdg[i] += hessin[i][j]*dg[j];
		}
		fac=fae=sumdg=sumxi=0.0;
		for (i=1;i<=n;i++) {
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
			sumdg += SQR(dg[i]);
			sumxi += SQR(xi[i]);
		}
		if (fac*fac > EPS*sumdg*sumxi) {
			fac=1.0/fac;
			fad=1.0/fae;
			for (i=1;i<=n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
			for (i=1;i<=n;i++) {
				for (j=1;j<=n;j++) {
					hessin[i][j] += (double)(fac*xi[i]*xi[j]
							-fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j]);
				}
			}
		}
		for (i=1;i<=n;i++) {
			xi[i]=0.0;
			for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];
		}
	}
	nrerror("too many iterations in dfpmin");
	/*FREEALL*/
}


#define GET_PSUM \
	for (j=1;j<=ndim;j++) {\
		for (sum=0.0,i=1;i<=mpts;i++) sum += p[i][j];\
		psum[j]=sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

double amotry(double **p, double y[], double psum[], int ndim,
		double (*funk)(double [], double []), int ihi, double fac, double *otherstuff)
{
	int j;
	double fac1,fac2,ytry,*ptry;

	ptry=dvector(1,ndim);
	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;
	for (j=1;j<=ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	ytry=(*funk)(ptry,otherstuff);
	if (ytry < y[ihi]) {
		y[ihi]=ytry;
		for (j=1;j<=ndim;j++) {
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}
	free_dvector(ptry,1);
	return ytry;
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */

void amoeba(double **p, double y[], int ndim, double ftol,
		double (*funk)(double [], double []), int *nfunk, double *otherstuff)
{
	double amotry(double **p, double y[], double psum[], int ndim,
			double (*funk)(double [], double []), int ihi, double fac, double otherstuff[]);
	int i,ihi,ilo,inhi,j,mpts=ndim+1;
	double rtol,sum,swap,ysave,ytry,*psum;

	psum=dvector(1,ndim);
	*nfunk=0;
	GET_PSUM
		for (;;) {
			ilo=1;
			ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
			for (i=1;i<=mpts;i++) {
				if (y[i] <= y[ilo]) ilo=i;
				if (y[i] > y[ihi]) {
					inhi=ihi;
					ihi=i;
				} else if (y[i] > y[inhi] && i != ihi) inhi=i;
			}
			rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
			//printf("\t\trtol: %.16f, ftol: %.16f\n", rtol, ftol);
			// instead of resulting in error, return current best
			if (rtol < ftol || *nfunk >= NMAX) {
				SWAP(y[1],y[ilo])
					for (i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i])
						break;
			}
			//if (*nfunk >= NMAX) nrerror("NMAX exceeded");		
			//if (*nfunk >= NMAX)
			//{
			//	SWAP(y[1],y[ilo])
			//		for (i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i])
			//			break;
			//	//break;		//Not sure if this is the correct move, but need to limit the # of iterations. Check with Rasmus later
			//}
			*nfunk += 2;
			ytry=amotry(p,y,psum,ndim,funk,ihi,-1.0, otherstuff);
			if (ytry <= y[ilo])
				ytry=amotry(p,y,psum,ndim,funk,ihi,2.0, otherstuff);
			else if (ytry >= y[inhi]) {
				ysave=y[ihi];
				ytry=amotry(p,y,psum,ndim,funk,ihi,0.5, otherstuff);
				if (ytry >= ysave) {
					for (i=1;i<=mpts;i++) {
						if (i != ilo) {
							for (j=1;j<=ndim;j++)
								p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
							y[i]=(*funk)(psum,otherstuff);
						}
					}
					*nfunk += ndim;
					GET_PSUM
				}
			} else --(*nfunk);
		}
	free_dvector(psum,1);
}

/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */

void Yanggradient (int n, double x[], double f0, double g[],
		double (*fun)(double x[], double z[]), double space[], int central, double lowbound[], double upbound[], double otherstuff[])
{

	/*f0=fun(x) is given for Central=0*/

	int i,j;
	double *x0=space, *x1=space+n, eh0=1e-8, eh01=1e-8, eh;

	if (central) {
		for (i=1;i<=n;i++)  {
			for (j=1;j<=n;j++)  x0[j]=x1[j]=x[j];
			eh=pow(eh01*(fabs(x[i])+1), 0.67);
			x0[i]-=eh; x1[i]+=eh;
			if (x0[i]<lowbound[i])
			{x1[i]+=eh; g[i] = ((*fun)(x1,otherstuff)-f0)/(eh*2.0);}
			else if (x1[i]>upbound[i])
			{x0[i]-=eh; g[i] = (f0-(*fun)(x0,otherstuff))/(eh*2.0);}
			else
				g[i] = ((*fun)(x1,otherstuff) - (*fun)(x0,otherstuff))/(eh*2.0);
			if (x[i] <= lowbound[i] && g[i] > 0.0)
				g[i] = 0.0;
			else if (x[i] >= upbound[i] && g[i] < 0.0)
				g[i] = 0.0; 
		}
	}
	else {/* (C) Copr. 1986-92 Numerical Recipes Software '$&'3$. */

		for (i=1;i<=n;i++)  {
			for (j=1;j<=n;j++)  
				x1[j]=x[j];
			/*eh=eh0*(fabs(x[i])+1);*/
			eh=2.0*pow(eh0*(fabs(x[i])+1), 0.67);
			if (x1[i]+eh>upbound[i])
			{
				x1[i]-=eh;
				g[i] = (f0-(*fun)(x1, otherstuff))/eh;
			}
			else
			{
				x1[i]+=eh;
				g[i] = ((*fun)(x1, otherstuff)-f0)/eh;
			}
			if (x[i] <= lowbound[i] && g[i] > 0.0)
				g[i] = 0.0;
			else if (x[i] >= upbound[i] && g[i] < 0.0)
				g[i] = 0.0;
		}
	}
}


double oldf0;
void getgradient(double invec[], double outvec[], double lowbound[], double upbound[], double(*func)(double [], double[]),double otherstuff[])
{
	int i;
	static double space[200];
	double f0, arbitrarysum = 0.0;

	f0 = func(invec,otherstuff);
	if (CENTRALMODE > 0)
		Yanggradient(npar, invec, f0, outvec, func, space, 1, lowbound, upbound, otherstuff);
	else 
		Yanggradient(npar, invec, f0, outvec, func, space, 0, lowbound, upbound, otherstuff);
	//printf("GRADIENT:\n");
	/*if (NULL==(tempfile=fopen("tempfile","w"))){
	  puts ("Kan ikke aabne tempfile!");
	  exit(-1);}*/
	for(i=1; i<=npar; i++)
	{
		/*	fprintf(tempfile,"%f ",invec[i]);  */
		//	printf(" %5f",outvec[i]);
		if (invec[i] > lowbound[i] && invec[i] < upbound[i])
			arbitrarysum = arbitrarysum + fabs(outvec[i]);
		/*else printf("*");*/
	}
	/* 	fclose(tempfile);*/

	//printf("%lf < %lf\n", arbitrarysum, 0.001 * npar);
	if (CENTRALMODE == 0 && fabs(oldf0-f0) < 0.1)
	{
		CENTRALMODE = 1;
		/*printf("Changing to central method\n");*/
	}
	else if (CENTRALMODE > 0 && fabs(oldf0-f0) > 1.0)
		CENTRALMODE = 0;
	else if (arbitrarysum < 0.001*npar)
		CENTRALMODE = 2;
	oldf0 = f0;
	//printf(" Like (grad): %f\n",-oldf0);
}

/*Call nrinits before calling this the first time*/
double findmax(double newinvecter[], double lowbound[], double upbound[], int n, double (*fun)(double x[],double z[]), double otherstuff[])

{

	double fret;
	int i, iter, whileCounter;

	// Due to potential infinite while loop, to remove in future hopefully
	whileCounter = 0;

	CENTRALMODE = 0;
	oldf0 = 0.0;
	npar = n;
	do    {
		dfpmin(newinvecter, npar, 0.0000000001, &iter, &fret, fun, getgradient, lowbound, upbound, otherstuff);
		whileCounter++;
	}while (CENTRALMODE < 2 && whileCounter < 10);
	return -fret;
}

//put bounds in otherstuff.  Assume first n entries are the bounds
double findmax_amoeba(double newinvecter[], double lowbound[], double upbound[], int n, double (*fun)(double x[],double z[]), double in_otherstuff[], int n_otherstuff)
{
	double **newinmat, *like, *otherstuff, L, eh0=1e-8;
	int nf, i, j, k;

	//printf("findmax_amoeba called with %i parameters\n",n);

	//first make one vector that contains bounds and other parameters
	otherstuff = malloc((2*n+n_otherstuff)*(sizeof(double)));
	for (i=0; i<n_otherstuff; i++)
		otherstuff[i] = in_otherstuff[i];
	for (i=0; i<n; i++){
		otherstuff[2*i+n_otherstuff] = lowbound[i+1]; //the bounds start counting at 1
		otherstuff[2*i+n_otherstuff+1] = upbound[i+1]; //the bounds start counting at 1
							       //printf("Bounds: %lf %lf\n",otherstuff[2*i+n_otherstuff] ,otherstuff[2*i+n_otherstuff+1]);
	}

	//initialize the simplex
	// Very rough multiple restart, work to make better
	like = malloc((n+2)*(sizeof(double)));
	newinmat = malloc((n+2)*(sizeof(double *)));
	for (i=1; i<n+2; i++)
		newinmat[i] = malloc((n+1)*sizeof(double));
	//for (i=1; i<n+2; i++){
	//	//printf("Point %i: ",i);
	//	for (j=1; j<=n; j++){
	//		if (i==j+1 && 2*(i/2)==i) newinmat[i][j] = lowbound[j]+eh0;
	//		else if (i==j+1) newinmat[i][j] = upbound[j]-eh0;
	//		else newinmat[i][j] = (lowbound[j]+upbound[j])/2.0;
	//		//printf("%lf ",newinmat[i][j]);
	//	}
	//	like[i] = fun(newinmat[i],otherstuff);	
	//	//	printf("Like: %lf\n",like[i]);
	//}

	// Hard coding initialization of vertices of the simplex
	// Case 1
	newinmat[1][1] = lowbound[1] + eh0;
	newinmat[1][2] = lowbound[2] + eh0;

	newinmat[2][1] = upbound[1] - eh0;
	newinmat[2][2] = lowbound[2] + eh0;

	newinmat[3][1] = (lowbound[1]+upbound[1])/2.0;
	newinmat[3][2] = upbound[2] - eh0;

	like[1] = fun(newinmat[1],otherstuff);
	like[2] = fun(newinmat[2],otherstuff);
	like[3] = fun(newinmat[3],otherstuff);

	amoeba(newinmat, like, n, 0.00000001, fun, &nf, otherstuff);
	L = like[1];
	for (i=1; i<n+1; i++)
		newinvecter[i] = newinmat[1][i];

	//Case 2
	newinmat[1][1] = upbound[1] - eh0;
	newinmat[1][2] = upbound[2] - eh0;

	newinmat[2][1] = lowbound[1] + eh0;
	newinmat[2][2] = upbound[2] - eh0;

	newinmat[3][1] = (lowbound[1]+upbound[1])/2.0;
	newinmat[3][2] = lowbound[2] + eh0;

	like[1] = fun(newinmat[1],otherstuff);
	like[2] = fun(newinmat[2],otherstuff);
	like[3] = fun(newinmat[3],otherstuff);

	amoeba(newinmat, like, n, 0.00000001, fun, &nf, otherstuff);
	if (like[1] < L)
	{
		L = like[1];
		for (i=1; i<n+1; i++)
			newinvecter[i] = newinmat[1][i];
	}

	//Case 3
	newinmat[1][1] = lowbound[1] + eh0;
	newinmat[1][2] = upbound[2] - eh0;

	newinmat[2][1] = lowbound[1] + eh0;
	newinmat[2][2] = lowbound[2] + eh0;

	newinmat[3][1] = upbound[1] - eh0;
	newinmat[3][2] = (lowbound[2]+upbound[2])/2.0;

	like[1] = fun(newinmat[1],otherstuff);
	like[2] = fun(newinmat[2],otherstuff);
	like[3] = fun(newinmat[3],otherstuff);

	amoeba(newinmat, like, n, 0.00000001, fun, &nf, otherstuff);
	if (like[1] < L)
	{
		L = like[1];
		for (i=1; i<n+1; i++)
			newinvecter[i] = newinmat[1][i];
	}
	//Case 4
	newinmat[1][1] = upbound[1] - eh0;
	newinmat[1][2] = upbound[2] - eh0;

	newinmat[2][1] = upbound[1] - eh0;
	newinmat[2][2] = lowbound[2] + eh0;

	newinmat[3][1] = lowbound[1] + eh0;
	newinmat[3][2] = (lowbound[2]+upbound[2])/2.0;

	like[1] = fun(newinmat[1],otherstuff);
	like[2] = fun(newinmat[2],otherstuff);
	like[3] = fun(newinmat[3],otherstuff);	

	amoeba(newinmat, like, n, 0.00000001, fun, &nf, otherstuff);
	if (like[1] < L)
	{
		L = like[1];
		for (i=1; i<n+1; i++)
			newinvecter[i] = newinmat[1][i];
	}


	for (i=1; i<n+2; i++)
		free(newinmat[i]);
	free(newinmat);
	free(otherstuff);
	free(like);
	return L;
}

//put bounds in otherstuff.  Assume first n entries are the bounds
double findmax_amoeba_rand(double newinvecter[], double lowbound[], double upbound[], int n, double (*fun)(double x[],double z[]), double in_otherstuff[], int n_otherstuff)
{
	double **newinmat, *like, *otherstuff, L, eh0=1e-8;
	int nf, i, j, k;

	//first make one vector that contains bounds and other parameters
	otherstuff = malloc((2*n+n_otherstuff)*(sizeof(double)));
	for (i=0; i<n_otherstuff; i++)
		otherstuff[i] = in_otherstuff[i];
	for (i=0; i<n; i++){
		otherstuff[2*i+n_otherstuff] = lowbound[i+1]; //the bounds start counting at 1
		otherstuff[2*i+n_otherstuff+1] = upbound[i+1]; //the bounds start counting at 1
	}

	like = malloc((n+2)*(sizeof(double)));
	newinmat = malloc((n+2)*(sizeof(double *)));
	for (i=1; i<n+2; i++)
		newinmat[i] = malloc((n+1)*sizeof(double));

	// No tests if they are in a line or not
	L = INF;
	for(i = 0; i < 1000; i++)
	{
		// Case 1
		newinmat[1][1] = lowbound[1] + ((double)rand()/RAND_MAX) * (upbound[1] - lowbound[1]);
		newinmat[1][2] = lowbound[2] + ((double)rand()/RAND_MAX) * (upbound[2] - lowbound[2]);

		newinmat[2][1] = lowbound[1] + ((double)rand()/RAND_MAX) * (upbound[1] - lowbound[1]);
		newinmat[2][2] = lowbound[2] + ((double)rand()/RAND_MAX) * (upbound[2] - lowbound[2]);

		newinmat[3][1] = lowbound[1] + ((double)rand()/RAND_MAX) * (upbound[1] - lowbound[1]);
		newinmat[3][2] = lowbound[2] + ((double)rand()/RAND_MAX) * (upbound[2] - lowbound[2]);

		like[1] = fun(newinmat[1],otherstuff);
		like[2] = fun(newinmat[2],otherstuff);
		like[3] = fun(newinmat[3],otherstuff);

		amoeba(newinmat, like, n, 0.00000001, fun, &nf, otherstuff);

		if (like[1] < L)
		{
			L = like[1];
			for (i=1; i<n+1; i++)
				newinvecter[i] = newinmat[1][i];
		}
	}

	for (i=1; i<n+2; i++)
		free(newinmat[i]);
	free(newinmat);
	free(otherstuff);
	free(like);
	return L;
}

double transformValue(double value, double max)
{
	return(max/(max - value));
}

double transformBack_amoeba(double value, double max)
{
	return(max - (max/value));
}

double findmax_amoeba_rand_trans(double newinvecter[], double lowbound[], double upbound[], int n, double (*fun)(double x[],double z[]), double in_otherstuff[], int n_otherstuff)
{
	double **newinmat, *like, *otherstuff, L, eh0=1e-8;
	int nf, i, j, k;

	//first make one vector that contains bounds and other parameters
	otherstuff = malloc((2*n+n_otherstuff+2)*(sizeof(double)));
	for (i=0; i<n_otherstuff; i++)
		otherstuff[i] = in_otherstuff[i];
	for (i=0; i<n; i++){
		otherstuff[2*i+n_otherstuff] = transformValue(lowbound[i+1], upbound[i+1]); //the bounds start counting at 1
		otherstuff[2*i+n_otherstuff+1] = transformValue(upbound[i+1]-eh0, upbound[i+1]); //the bounds start counting at 1
	}

	// this is hardcoded, be careful
	otherstuff[7] = upbound[1];
	otherstuff[8] = upbound[2];

	like = malloc((n+2)*(sizeof(double)));
	newinmat = malloc((n+2)*(sizeof(double *)));
	for (i=1; i<n+2; i++)
		newinmat[i] = malloc((n+1)*sizeof(double));

	// No tests if they are in a line or not
	L = INF;
	for(i = 0; i < 10000; i++)
	{
		// Case 1
		newinmat[1][1] = otherstuff[3] + ((double)rand()/RAND_MAX) * (otherstuff[4] - otherstuff[3]);
		newinmat[1][2] = otherstuff[5] + ((double)rand()/RAND_MAX) * (otherstuff[6] - otherstuff[5]);

		newinmat[2][1] = otherstuff[3] + ((double)rand()/RAND_MAX) * (otherstuff[4] - otherstuff[3]);
		newinmat[2][2] = otherstuff[5] + ((double)rand()/RAND_MAX) * (otherstuff[6] - otherstuff[5]);

		newinmat[3][1] = otherstuff[3] + ((double)rand()/RAND_MAX) * (otherstuff[4] - otherstuff[3]);
		newinmat[3][2] = otherstuff[5] + ((double)rand()/RAND_MAX) * (otherstuff[6] - otherstuff[5]);

		like[1] = fun(newinmat[1],otherstuff);
		like[2] = fun(newinmat[2],otherstuff);
		like[3] = fun(newinmat[3],otherstuff);

		amoeba(newinmat, like, n, 0.00000001, fun, &nf, otherstuff);

		if (like[1] < L)
		{
			L = like[1];
			for (i=1; i<n+1; i++)
				newinvecter[i] = transformBack_amoeba(newinmat[1][i], upbound[i]);
		}
	}

	for (i=1; i<n+2; i++)
		free(newinmat[i]);
	free(newinmat);
	free(otherstuff);
	free(like);
	return L;
}

// Golden section search from Andrew Vaughn and Wikipedia
// This is the version that reuses function evaluations
// a is left bound, b is right bound, tol is final size of interval
// (double newinvecter[], double lowbound[], double upbound[], int n, double (*fun)(double x[],double z[]), double in_otherstuff[], int n_otherstuff)
// double GoldenSection(double a, double b, double tol) 
double GoldenSection(double newinvecter[], double lowbound[], double upbound[], int n, double (*fun)(double x,double z[]), double in_otherstuff[], int n_otherstuff)
{
	int i;
	double a = lowbound[1], b = upbound[1];
	double invphi = (sqrt(5) - 1) / 2;
	double invphi2 = (3 - sqrt(5)) / 2;
	double h = b - a;
	double c = a + invphi2 * h;
	double d = a + invphi * h;
	double yc, yd, tol = 0.00000001; 	// same as in ameoba
	double *otherstuff;

	//first make one vector that contains bounds and other parameters
	otherstuff = malloc((2*n+n_otherstuff)*(sizeof(double)));
	for (i=0; i<n_otherstuff; i++)
	{
		otherstuff[i] = in_otherstuff[i];
	}
	for (i=0; i<n; i++)
	{
		otherstuff[2*i+n_otherstuff] = lowbound[i+1]; //the bounds start counting at 1
		otherstuff[2*i+n_otherstuff+1] = upbound[i+1]; //the bounds start counting at 1
							       //printf("Bounds: %lf %lf\n",otherstuff[2*i+n_otherstuff] ,otherstuff[2*i+n_otherstuff+1]);
	}

	yc = (*fun)(c, otherstuff);
	yd = (*fun)(d, otherstuff);

	while ( h > tol ){
		if (yc < yd) 
		{
			b = d;
			d = c;
			yd = yc;
			h = invphi * h;
			c = a + invphi2 * h;
			yc = (*fun)(c, otherstuff);
		}
		else 
		{
			a = c;
			c = d;
			yc = yd;
			h = invphi * h;
			d = a + invphi * h;
			yd = (*fun)(d, otherstuff);
		}
	}


	if (yc < yd) 
	{
		newinvecter[1] = (a + d)/2;
		free(otherstuff);
		return(yc);
	}
	else
	{
		newinvecter[1] = (b + c)/2;
		free(otherstuff);
		return(yd);
	}

	free(otherstuff);
}

// Increase the tolerance value to make this a less stringent golden section to make faster for calls of golden section that don't need high accuracy of the likelihood
double GoldenSection_rough(double newinvecter[], double lowbound[], double upbound[], int n, double (*fun)(double x,double z[]), double in_otherstuff[], int n_otherstuff)
{
	int i;
	double a = lowbound[1], b = upbound[1];
	double invphi = (sqrt(5) - 1) / 2;
	double invphi2 = (3 - sqrt(5)) / 2;
	double h = b - a;
	double c = a + invphi2 * h;
	double d = a + invphi * h;
	double yc, yd, tol = 0.000001; 	// same as in ameoba
	double *otherstuff;

	//first make one vector that contains bounds and other parameters
	otherstuff = malloc((2*n+n_otherstuff)*(sizeof(double)));
	for (i=0; i<n_otherstuff; i++)
	{
		otherstuff[i] = in_otherstuff[i];
		//printf("otherstuff %d %.16f\n", i, otherstuff[i]); 
	}
	for (i=0; i<n; i++)
	{
		otherstuff[2*i+n_otherstuff] = lowbound[i+1]; //the bounds start counting at 1
		otherstuff[2*i+n_otherstuff+1] = upbound[i+1]; //the bounds start counting at 1
							       //printf("Bounds: %lf %lf\n",otherstuff[2*i+n_otherstuff] ,otherstuff[2*i+n_otherstuff+1]);
	}

	yc = (*fun)(c, otherstuff);
	yd = (*fun)(d, otherstuff);

	while ( h > tol ){
		if (yc < yd) 
		{
			b = d;
			d = c;
			yd = yc;
			h = invphi * h;
			c = a + invphi2 * h;
			yc = (*fun)(c, otherstuff);
		}
		else 
		{
			a = c;
			c = d;
			yc = yd;
			h = invphi * h;
			d = a + invphi * h;
			yd = (*fun)(d, otherstuff);
		}
	}


	if (yc < yd)
	{
		newinvecter[1] = (a + d)/2;
		return(yc);
	}
	else
	{
		newinvecter[1] = (b + c)/2;
		return(yd);
	}
}

/*
 * Adapter for Brent1D: bridges GoldenSection-style callback
 * double (*fun)(double x, double z[]) to MinimizeFxn double (*)(double *, void *).
 */
typedef struct {
	double (*fun)(double x, double z[]);
	double *otherstuff;
} Brent1D_ctx;

static double Brent1D_bridge(double *x, void *extra_data)
{
	Brent1D_ctx *ctx = (Brent1D_ctx *)extra_data;
	return ctx->fun(*x, ctx->otherstuff);
}

/*
 * 1D Brent's method (LocalMin) with the same calling convention as GoldenSection.
 * Faster on smooth functions due to quadratic interpolation (superlinear convergence).
 * Drop-in replacement: result stored in newinvecter[1], returns minimum value.
 */
double Brent1D(double newinvecter[], double lowbound[], double upbound[], int n,
               double (*fun)(double x, double z[]), double in_otherstuff[], int n_otherstuff)
{
	int i;
	double px, result;

	/* Build the same combined otherstuff array layout as GoldenSection uses */
	double *otherstuff = malloc((2 * n + n_otherstuff) * sizeof(double));
	for (i = 0; i < n_otherstuff; i++)
		otherstuff[i] = in_otherstuff[i];
	for (i = 0; i < n; i++) {
		otherstuff[2 * i + n_otherstuff]     = lowbound[i + 1];
		otherstuff[2 * i + n_otherstuff + 1] = upbound[i + 1];
	}

	Brent1D_ctx ctx = {fun, otherstuff};
	result = LocalMin(lowbound[1], upbound[1], sqrt(DBL_EPSILON), 1e-8,
	                  Brent1D_bridge, &px, &ctx);

	newinvecter[1] = px;
	free(otherstuff);
	return result;
}

//double minimize_brent(double newinvecter[], int n, double (*fun)(double x[]), int maxIterations) {
//        double minusLnL;
//
//        double *directions = (double*)malloc(sizeof(double) * n * n);
//        if (!directions)
//        {
//                printf ("Could not allocate directions (%lu)\n", sizeof(double) * n * n);
//                exit (1);
//        }
//        
//        double *powellWork = (double*)malloc(sizeof(double) * 6 * n);
//        if (!powellWork)
//        {
//                printf ("Could not allocate powellWork (%lu)\n", sizeof(double) * 6 * n);
//                exit (1);
//        }
//
//        minusLnL = PrAxis(TOLER_PASS_2, MAX_STEP_SIZE_2, n, newinvecter, *fun, directions, powellWork, maxIterations);
//
//		free(directions);
//		free(powellWork);
//
//        return(minusLnL);
//}

double minimize_brent(double newinvecter[], int n, double (*fun)(double x[], void*), int maxIterations, void *extra_data) {
    double minusLnL;

    double *directions = (double*)malloc(sizeof(double) * n * n);
    if (!directions)
    {
        printf ("Could not allocate directions (%lu)\n", sizeof(double) * n * n);
        exit (1);
    }
    
    double *powellWork = (double*)malloc(sizeof(double) * 6 * n);
    if (!powellWork)
    {
        printf ("Could not allocate powellWork (%lu)\n", sizeof(double) * 6 * n);
        exit (1);
    }

    minusLnL = PrAxis(TOLER_PASS_2, MAX_STEP_SIZE_2, n, newinvecter, *fun, directions, powellWork, maxIterations, extra_data);

    free(directions);
    free(powellWork);

    return(minusLnL);
}