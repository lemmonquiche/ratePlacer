/*	minfunc.c
|
|	Routines for function minimization based on those of Brent (1973; "Algorithms for
|	Minimization Without Derivatives"; Prentice-Hall).  Brent's routines were translated from
|	Algol to C and cleaned up (goto elimination, provision for dynamic memory allocation, etc.).
|
|	Copyright (c) 1996 by David L. Swofford, Smithsonian Institution.
|	All rights reserved.
*/
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "minfunc.h"
#include "jph.h"


static void   LineMin (int j, int nits, double *pd2, double *px1, double f1, int fk);
static void   Quad (void);
static void   MinFit (double eps, double tol, double ab[], double q[], double e[]);
static void   SortDV (void);
static double FLin (int j, double lambda);
static void   SetIdentityMatrix2 (double a[], int n);

#define MAX_ITER	50
#define CGOLD		0.3819660113		/* = (3 - sqrt(5))/2 */
#define GOLD		1.618034
#define GLIMIT		100.0
#define TINY		1.0e-20
#define SIGN(a, b)	((b) > 0.0 ? fabs(a) : -fabs(a))
#if !defined(MAX)
#	define MAX(x, y)	((x) >= (y)	? (x) : (y))
#endif

/*--------------------------------------------------------------------------------------------------
|
|	LocalMin
|
|	Find a local minimum of a function of one variable using Brent's (1973) method.
|
|	This function is reentrant.
*/

double LocalMin (double a, double b, double eps, double t, MinimizeFxn f, double *px)
	/* double		a, b;	 on input, (a,b) must bracket a local minimum */
	/* double		eps;	 t and eps define tol = eps|x| + t, f is never evaluated at two points */
	/* double		t;		 closer together than tol;  eps should be > sqrt(DBL_EPSILON)        */
	/* MinimizeFxn	f;		 function to minimize */
	/* double		*px;	 value of x when f(x) is minimal */
{

	double		e, m, p, q, r, x, tol, t2, u, v, w, fu, fv, fw, fx,
				d = 0.0;	/* shuts up bogus lint warning */
	int			iter;

	v  = w = x = a + CGOLD*(b - a);
	e  = 0.0;
	fv = fw = fx = (*f)(&x);
	
	/* main loop */

	for (iter = 0; iter < MAX_ITER; iter++)
		{
		m = 0.5*(a + b);
		tol = eps*fabs(x) + t;
		t2 = 2.0*tol;

		/* check stopping criterion */
		if (fabs(x - m) <= t2 - 0.5*(b - a))
			break;

		p = q = r = 0.0;
		if (fabs(e) > tol)
			{
			/* fit parabola (trial) */
			r = (x - w)*(fx - fv);
			q = (x - v)*(fx - fw);
			p = (x - v)*q - (x - w)*r;
			q = 2.0*(q - r);
			if (q > 0.0)
				p = -p;
			else
				q = -q;
			r = e;
			e = d;	/* lint complains about possible use of 'd' before being set, but it's not
			           a problem because e=0 first time through, so "fabs(e) > tol" test fails */
			}
			
		/* Take parabolic-interpolation or golden-section step (note that Brent's Algol procedure
		   had an error, p > q*(a-x) is correct) */
		if ((fabs(p) < fabs(0.5*q*r)) && (p > q*(a - x)) && (p < q*(b - x)))
			{
			/* parabolic interpolation step */
			d = p/q;
			u = x + d;
			/* don't evaluate f too close to a or b */
			if ((u - a < t2) || (b - u < t2))
				d = (x < m) ? tol : -tol;
			}
		else
			{
			/* "golden section" step */
			e = ((x < m) ? b : a) - x;
			d = CGOLD*e;
			}

		/* don't evaluate f too close to x */
		if (fabs(d) >= tol)
			u = x + d;
		else if (d > 0.0)
			u = x + tol;
		else 
			u = x - tol;

		fu = (*f)(&u);
		
		/* update a, b, v, w, and x */
		if (fu <= fx)
			{
			if (fu <= fx)
				{
				if (u < x)
					b = x;
				else
					a = x;
				v  = w;
				fv = fw;
				w  = x;
				fw = fx;
				x  = u;
				fx = fu;
				}
			}
		else
			{
			if (u < x)
				a = u;
			else
				b = u;
			if ((fu <= fw) || (w == x))
				{
				v  = w;
				fv = fw;
				w  = u;
				fw = fu;
				}
			else if ((fu <= fv) || (v == x) || (v == w))
				{
				v  = u;
				fv = fu;
				}
			}
		}
	*px = x;
	return fx;
	
}

/*--------------------------------------------------------------------------------------------------
|
|	PrAxis
|
|	Minimize a function f of n variables using Brent's (1973) principal axis method.  Translated
|	from the Algol program in Brent's book.
|
|	Note: this function is NOT reentrant, due to the need for the globals below.
*/

/* following variables have scope of function 'PrAxis' in Algol, but must be global to file in C */
static long int		praxisSeed;
static int			nl;				/* number of line minimizations performed */
static double		dmin;
static double		ldt;
static double		qf1;
static double		qd0, qd1;
static double		m2, m4;
static double		eps2;
static double		toler, htol;	/* global versions of PrAxis arguments */
static MinimizeFxn	fxn;			/* objective function to minimize */
static int			n;
static double		*vv;
static double		*gx, gfx;
static double		*d, *q0, *q1;
static double		*xnew;

double PrAxis (double tol, double h, int nn, double xx[], MinimizeFxn f, double v[], double work[])
	/* double			tol;	tolerance used for convergence criterion */
	/* double			h;	 	maximum step size */
	/* int				nn;		number of variables */
	/* double			*xx;	on input, must contain initial guess for optimal point */
	/* MinimizeFxn		f;		the function to be minimized */
	/* double			**v;	n x n matrix in row-ptr form */
	/* double			*work;	work vector of size 6*nn */
{
	int				i, j, k, k2, illc, kl, kt, ktm;
	double			vsmall, large, vlarge, scbd, ldfac, sf, df, f1, lds, t2, sl, dn, s, sz, *y, *z;

	/* copy args to global equivalents */
	fxn   = f;
	toler = tol;
	htol  = h;
	n     = nn;
	gx    = xx;
	vv    = v;
	
	/* partition work space provided by caller into arrays used here */
	d    = work;
	q0   = d + nn;
	q1   = q0 + nn;
	xnew = q1 + nn;
	y    = xnew + nn;
	z    = y + nn;

	/* machine-dependent initializations */
	eps2   = DBL_EPSILON*DBL_EPSILON;
	vsmall = eps2*eps2;
	large  = 1.0/eps2;
	vlarge = 1.0/vsmall;
	m2     = sqrt(DBL_EPSILON);
	m4     = sqrt(m2);
	
	praxisSeed = 1;	/* starting seed for random number generator */
	
	/* heuristic numbers:
	   - if axes may be badly scaled (which should be avoided if possible), set scbd=10, otherwise 1
	   - if the problem is known to be ill-conditioned, set illc=TRUE, otherwise FALSE
	   - ktm+1 is the number of iterations without improvement before the algorithm terminates
	     (see section 7.6 of Brent's book).  ktm=4 is very cautious; usually ktm=1 is satisfactory.
	*/
	scbd = 1.0;
	illc = FALSE;
	ktm  = 4;

	ldfac = illc ? 0.1 : 0.01;
	kt    = nl = 0;
	qf1   = gfx = (*f)(gx);
	toler = t2 = eps2 + fabs(toler);
	dmin  = eps2;
	if (htol < 100.0*toler)
		htol = 100.0*toler;
	ldt = htol;
	SetIdentityMatrix2 (v, n);

	d[0] = qd0 = 0.0;
	for (i = 0; i < n; i++)
		{
		q0[i] = 0.0;	/* DLS: this wasn't included in Brent's code, but it's important */
		q1[i] = gx[i];
		}

	/* ------ main loop ------ */
	
	for (;;)
		{
		sf   = d[0];
		d[0] = s = 0.0;
		
		/* minimize along first direction */
	
		LineMin (0, 2, &d[0], &s, gfx, FALSE);
		if (s < 0.0)
			{
			for (i = 0; i < n; i++)
				v[pos1(i,0,nn)] = -v[pos1(i,0,nn)];
			}
		if ((sf <= 0.9*d[0]) || (0.9*sf >= d[0]))
			{
			for (i = 1; i < n; i++)
				d[i] = 0.0;
			}
	
		for (k = 1; k < n; k++)
			{
			for (i = 0; i < n; i++)
				y[i] = gx[i];
			sf   = gfx;
			illc = illc || (kt > 0);
	
			for (;;)
				{
				kl = k;
				df = 0.0;
				if (illc)
					{
					/* random step to get off resolution valley */
					for (i = 0; i < n; i++)
						{
						s = z[i] = (0.1*ldt + t2*pow(10.0, kt))*(RandomNumber(&praxisSeed) - 0.5);
						for (j = 0; j < n; j++)
							gx[j] += s*v[pos1(j,i,nn)];
						}
					gfx = (*f)(gx);
					}
				for (k2 = k; k2 < n; k2++)
					{
					sl = gfx;
					s  = 0.0;
					/* minimize along "non-conjugate" directions */
					LineMin(k2, 2, &d[k2], &s, gfx, FALSE);
					if (illc)
						{
						sz = s + z[k2];
						s  = d[k2]*sz*sz;
						}
					else
						s = sl - gfx;
					if (df < s)
						{
						df = s;
						kl = k2;
						}
					}
				if (!illc && (df < fabs(100.0*DBL_EPSILON*gfx)))
					illc = TRUE;	/* no success with illc=FALSE so try once with illc=TRUE */
				else
					break;
				}
			
			for (k2 = 0; k2 < k; k2++)
				{
				/* minimize along "conjugate" directions */
				s = 0.0;
				LineMin(k2, 2, &d[k2], &s, gfx, FALSE);
				}
	
			f1  = gfx;
			gfx = sf;
			lds = 0.0;
			for (i = 0; i < n; i++)
				{
				sl = gx[i];
				gx[i] = y[i];
				y[i] = (sl -= y[i]);
				lds += sl*sl;
				}
			lds = sqrt(lds);
			if (lds > eps2)
				{
				/* throw away direction kl */
				for (i = kl - 1; i >= k; i--)
					{
					for (j = 0; j < n; j++)
						v[pos1(j,i+1,nn)] = v[pos1(j,i,nn)];
					d[i + 1] = d[i];
					}
					
				/* set new "conjugate" direction ... */
				d[k] = 0.0;
				for (i = 0; i < n; i++)
					v[pos1(i,k,nn)] = y[i]/lds;
				
				/* ... and minimize along it */
				LineMin(k, 4, &d[k], &lds, f1, TRUE);
				if (lds <= 0.0)
					{
					lds = -lds;
					for (i = 0; i < n; i++)
						v[pos1(i,k,nn)] = -v[pos1(i,k,nn)];
					}
				}
			ldt *= ldfac;
			if (ldt < lds)
				ldt = lds;
			t2 = 0.0;
			for (i = 0; i < n; i++)
				t2 += gx[i]*gx[i];
			t2 = m2*sqrt(t2) + toler;
			
			/* see if step length exceeds half the tolerance */
			kt = (ldt > 0.5*t2) ? 0 : kt + 1;
			if (kt > ktm)
				{
				return gfx;
				}
			}
		
		/* try quadratic extrapolation in case we are stuck in a curved valley */		
		Quad();
	
		/* calculate V = U.(D^(-1/2))  (note: 'v' currently contains U) */
	
		dn = 0.0;
		for (i = 0; i < n; i++)
			{
			d[i] = 1.0/sqrt(d[i]);
			if (d[i] > dn)
				dn = d[i];
			}
		for (j = 0; j < n; j++)
			{
			s = d[j]/dn;
			for (i = 0; i < n; i++)
				v[pos1(i,j,nn)] *= s;
			}
	
		if (scbd > 1.0)
			{
			/* scale axes in attempt to reduce condition number */
			s = vlarge;
			for (i = 0; i < n; i++)
				{
				sl = 0.0;
				for (j = 0; j < n; j++)
					sl += v[pos1(i,j,nn)]*v[pos1(i,j,nn)];
				z[i] = sqrt(sl);
				if (z[i] < m4)
					z[i] = m4;
				if (s > z[i])
					s = z[i];
				}
			for (i = 0; i < n; i++)
				{
				sl   = s/z[i];
				z[i] = 1.0/sl;
				if (z[i] > scbd)
					{
#					if 0	/* DLS/POL 1/6/97: Borland C says this is a do-nothing statement,
                                               and it looks like they're right */
					sl = 1.0/scbd;
#					endif
					z[i] = scbd;
					}
				}
			}
			
		/* Find the singular value decomposition of v.  This gives the eigenvalues and principal axes
		   of the approximating quadratic form without squaring the condition number. */
		/* transpose v for MinFit */
		for (i = 1; i < n; i++)
			{
			for (j = 0; j < i; j++)
				{
				s = v[pos1(i,j,nn)];
				v[pos1(i,j,nn)] = v[pos1(j,i,nn)];
				v[pos1(j,i,nn)] = s;
				}
			}
		MinFit (DBL_EPSILON, vsmall, v, d, y);	/* ('y' is just a scratch vector) */
	
		if (scbd > 1.0)
			{
			/* unscaling */
			for (i = 0; i < n; i++)
				{
				s = z[i];
				for (j = 0; j < n; j++)
					v[pos1(i,j,nn)] *= s;
				}
			for (i = 0; i < n; i++)
				{
				s = 0.0;
				for (j = 0; j < n; j++)
					s += v[pos1(j,i,nn)]*v[pos1(j,i,nn)];
				s = sqrt(s);
				d[i] *= s;
				s = 1.0/s;
				for (j = 0; j < n; j++)
					v[pos1(j,i,nn)] *= s;
				}
			}
	
		for (i = 0; i < n; i++)
			{
			s = dn*d[i];
			if (s > large)
				d[i] = vsmall;
			else if (s < eps2)
				d[i] = vlarge;
			else
				d[i] = 1.0/(s*s);
			}
			
		/* sort new eigenvalues and eigenvectors */
		SortDV();
		
		dmin = d[n - 1];
		if (dmin < eps2)
			dmin = eps2;
		
		illc = (m2*d[0] > dmin);
		}
}





/*--------------------------------------------------------------------------------------------------
|
|	MinFit
|
|	Get singular value decomposition using a modified version of Golub and Reinsch's (1969)
|	routine, restricted to m=n.  The singular values of the array 'ab' are returned in 'q', and
|	'ab' is overwritten with the orthogonal matrix V such that U.diag(Q) = AB.V, where U is another
|	orthogonal matrix.
*/

void MinFit (double eps, double tol, double ab[], double q[], double e[])
	/*double		eps;
	double		tol;
	double		**ab;
	double		*q;
	double		*e;	*/	/* work vector of size n */
{
	int			i, j, k, l, l2, kt;
	double		c, f, g, h, s, x, y, z;

	/* Householder's reduction to bidiagonal form */
	g = x = 0.0;
	for (i = 0; i < n; i++)
		{
		e[i] = g;
		s = 0.0;
		l = i + 1;
		for (j = i; j < n; j++)
			s += ab[pos1(j,i,n)]*ab[pos1(j,i,n)];
		if (s < tol)
			g = 0.0;
		else
			{
			f = ab[pos1(i,i,n)];
			g = (f < 0.0) ? sqrt(s) : -sqrt(s);
			h = f*g - s;
			ab[pos1(i,i,n)] = f - g;
			for (j = l; j < n; j++)
				{
				f = 0.0;
				for (k = i; k < n; k++)
					f += ab[pos1(k,i,n)]*ab[pos1(k,j,n)];
				f /= h;
				for (k = i; k < n; k++)
					ab[pos1(k,j,n)] += f*ab[pos1(k,i,n)];
				}
			}
		q[i] = g;
		s = 0.0;
		for (j = l; j < n; j++)
			s += ab[pos1(i,j,n)]*ab[pos1(i,j,n)];
		
		if (s < tol)
			g = 0.0;
		else
			{
			f = ab[pos1(i,i+1,n)];
			g = (f < 0.0) ? sqrt(s) : -sqrt(s);
			h = f*g - s;
			ab[pos1(i,i+1,n)] = f - g;
			for (j = l; j < n; j++)
				e[j] = ab[pos1(i,j,n)]/h;
			for (j = l; j < n; j++)
				{
				s = 0.0;
				for (k = l; k < n; k++)
					s += ab[pos1(j,k,n)]*ab[pos1(i,k,n)];
				for (k = l; k < n; k++)
					ab[pos1(j,k,n)] += s*e[k];
				}
			}
		y = fabs(q[i]) + fabs(e[i]);
		if (y > x)
			x = y;
		}

	/* accumulation of right-hand transformations */
	for (i = n-1; i >= 0; i--)
		{
		if (g != 0.0)
			{
			h = ab[pos1(i,i+1,n)]*g;
			for (j = l; j < n; j++)
				ab[pos1(j,i,n)] = ab[pos1(i,j,n)]/h;
			for (j = l; j < n; j++)
				{
				s = 0.0;
				for (k = l; k < n; k++)
					s += ab[pos1(i,k,n)]*ab[pos1(k,j,n)];
				for (k = l; k < n; k++)
					ab[pos1(k,j,n)] += s*ab[pos1(k,i,n)];
				}
			}
		for (j = l; j < n; j++)
			ab[pos1(i,j,n)] = ab[pos1(j,i,n)] = 0.0;
		ab[pos1(i,i,n)] = 1.0;
		g = e[i];
		l = i;
		}
		
	/* diagonalization of the bidiagonal form */
	eps *= x;
	for (k = n-1; k >= 0; k--)
		{
		kt = 0;
		
		test_splitting:

		if (++kt > 30)
			{
			e[k] = 0.0;
			}
		for (l2 = k; l2 >= 0; l2--)
			{
			l = l2;
			if (fabs(e[l]) <= eps)
				goto test_convergence;
			if (fabs(q[l-1]) <= eps)
				break;
			}

		/* cancellation of e[l] if l > 1 */
		c = 0.0;
		s = 1.0;
		for (i = l; i <= k; i++)
			{
			f = s*e[i];
			e[i] *= c;
			if (fabs(f) <= eps)
				break;
			g = q[i];
			if (fabs(f) < fabs(g))
				h = fabs(g)*sqrt(1.0 + (f/g)*(f/g));
			else if (f != 0.0)
				h = fabs(f)*sqrt(1.0 + (g/f)*(g/f));
			else
				h = 0.0;
			q[i] = h;
			if (h == 0.0)
				g = h = 1.0;	/* note: this replaces q[i]=h=sqrt(g*g+f*f) which may give
				                   incorrect results if the squares underflow or if f=g=0 */
			c = g/h;
			s = -f/h;
			}
			
		test_convergence:
		
		z = q[k];
		if (l != k)
			{			
			/* shift from bottom 2*2 minor */
			
			x = q[l];
			y = q[k-1];
			g = e[k-1];
			h = e[k];
			f = ((y-z)*(y+z) + (g-h)*(g+h))/(2.0*h*y);
			g = sqrt(f*f + 1.0);
			s = (f < 0.0) ? f - g : f + g;
			f = ((x-z)*(x+z) + h*(y/s - h))/x;
				
			/* next QR transformation */
			
			c = s = 1.0;
			for (i = l + 1; i <= k; i++)
				{
				g = e[i];
				y = q[i];
				h = s*g;
				g *= c;
				if (fabs(f) < fabs(h))
					z = fabs(h)*sqrt(1.0 + (f/h)*(f/h));
				else if (f != 0.0)
					z = fabs(f)*sqrt(1.0 + (h/f)*(h/f));
				else
					z = 0.0;
				e[i-1] = z;
				if (z == 0.0)
					z = f = 1.0;
				c = f/z;
				s = h/z;
				f = x*c + g*s;
				g = -x*s + g*c;
				h = y*s;
				y *= c;
				for (j = 0; j < n; j++)
					{
					x = ab[pos1(j,i-1,n)];
					z = ab[pos1(j,i,n)];
					ab[pos1(j,i-1,n)] = x*c + z*s;
					ab[pos1(j,i,n)] = -x*s + z*c;
					}
				if (fabs(f) < fabs(h))
					z = fabs(h)*sqrt(1.0 + (f/h)*(f/h));
				else if (f != 0.0)
					z = fabs(f)*sqrt(1.0 + (h/f)*(h/f));
				else
					z = 0.0;
				q[i-1] = z;
				if (z == 0.0)
					z = f = 1.0;
				c = f/z;
				s = h/z;
				f = c*g + s*y;
				x = -s*g + c*y;
				}
			e[l] = 0.0;
			e[k] = f;
			q[k] = x;
			goto test_splitting;
			}
		
		if (z < 0.0)
			{
			/* q[k] is made non-negative */
			q[k] = -z;
			for (j = 0; j < n; j++)
				ab[pos1(j,k,n)] = -ab[pos1(j,k,n)];
			}
		}
}





/*--------------------------------------------------------------------------------------------------
|
|	SortDV
|
|	Sorts the elements of 'd' and corresponding elements of 'v' into descending order.
*/

void SortDV (void)
{
	int		i, j, k;
	double	s;

	for (i = 0; i < n - 1; i++)
		{
		k = i;
		s = d[i];
		for (j = i + 1; j < n; j++)
			{
			if (d[j] > s)
				{
				k = j;
				s = d[j];
				}
			}
		if (k > i)
			{
			d[k] = d[i];
			d[i] = s;
			for (j = 0; j < n; j++)
				{
				s = vv[pos1(j,i,n)];
				vv[pos1(j,i,n)] = vv[pos1(j,k,n)];
				vv[pos1(j,k,n)] = s;
				}
			}
		}
}

/*--------------------------------------------------------------------------------------------------
|
|	LineMin
|
|	Minimize the objective function from point x in the direction v(*,j) (if j >= 0), or
|	perform a quadratic search in the plane defined by q0, q1, and x (if j < 0).
*/

void LineMin (int j, int nits, double *pd2, double *px1, double f1, int fk)
	/* int		j;		column of direction matrix (or <0 flag for quadratic search) */
	/* int		nits;	number of times an attempt is made to halve the interval */
	/* double	*pd2;	approximation to half f'' (or zero) */
	/* double	*px1;	x1 = estimate of distance to minimum, returned as actual distance found */
	/* double	f1;		if fk=TRUE, FLin(x1), otherwise ignored */
	/* int		fk;		flag (see above) */
{
	int		i, k, need_d2z, success;
	double	x1, x2, xm, f0, f2, fm, d1, d2, t2, s, sf1, sx1;

	/* copy args passed by reference to locals (will pass back at end) */
	d2 = *pd2;
	x1 = *px1;
	
	sf1 = f1;
	sx1 = x1;
	k = 0;
	xm = 0.0;
	f0 = fm = gfx;
	need_d2z = (d2 < DBL_EPSILON);	/* if TRUE, we need f''(0) */
	
	/* find step size */
	s = 0.0;
	for (i = 0; i < n; i++)
		s += gx[i]*gx[i];
	s = sqrt(s);
	t2 = m4*sqrt(fabs(gfx)/(need_d2z ? dmin : d2) + s*ldt) + m2*ldt;
	s = m4*s + toler;
	if (need_d2z && (t2 > s))
		t2 = s;
	if (t2 < eps2)
		t2 = eps2;
	if (t2 > 0.01*htol)
		t2 = 0.01*htol;
	if (fk && (f1 <= fm))
		{
		xm = x1;
		fm = f1;
		}
	if (!fk || (fabs(x1) < t2))
		{
		x1 = (x1 >= 0.0) ? t2 : -t2;
		f1 = FLin(j, x1);
		}
	if (f1 <= fm)
		{
		xm = x1;
		fm = f1;
		}
	
	/* find a distance x2 ("lambda*") that approximately minimizes f in the chosen direction */
	
	do	{
		if (need_d2z)
			{
			/* evaluate FLin at another point and estimate the second derivative */
			x2 = (f0 < f1) ? -x1 : 2.0*x1;
			f2 = FLin(j, x2);
			if (f2 <= fm)
				{
				xm = x2;
				fm = f2;
				}
			d2 = (x2*(f1 - f0) - x1*(f2 - f0))/(x1*x2*(x1 - x2));
			}
	
		/* estimate first derivative at 0 */	
		d1 = (f1 - f0)/x1 - x1*d2;
		need_d2z = TRUE;			/* reset flag in case we don't exit loop */
		
		/* predict minimum */
		if (d2 <= eps2)
			x2  = (d1 < 0.0) ? htol : -htol;
		else
			x2 = -0.5*d1/d2;
		if (fabs(x2) > htol)
			x2 = (x2 > 0.0) ? htol : -htol;
		
		/* evaluate f at predicted minimum */
		do	{
			f2 = FLin(j, x2);
			success = TRUE;
			if ((k < nits) && (f2 > f0))
				{
				/* no success so halve interval and try again */
				success = FALSE;
				k++;
				if ((f0 < f1) && (x1*x2 > 0.0))
					break;
				x2 *= 0.5;
				}
			}
			while (!success);
		}		
		while (!success);

	nl++;	/* increment one-dimensional search counter */
	if (f2 > fm)
		x2 = xm;
	else
		fm = f2;
		
	/* get new estimate of second derivative */
	if (fabs(x2*(x2 - x1)) > eps2)
		d2 = (x2*(f1 - f0) - x1*(fm - f0))/(x1*x2*(x1 - x2));
	else if (k > 0)
		d2 = 0.0;
	if (d2 < eps2)
		d2 = eps2;
	x1 = x2;
	gfx = fm;
	if (sf1 < gfx)
		{
		gfx = sf1;
		x1 = sx1;
		}
		
	/* update x for linear search but not for parabolic search */
	if (j >= 0)
		{
		for (i = 0; i < n; i++)
			gx[i] += x1*vv[pos1(i,j,n)];
		}

	*px1 = x1;
	*pd2 = d2;
}





/*--------------------------------------------------------------------------------------------------
|
|	FLin
|
|	Evaluate a function of one variable after moving trial point a distance lambda from the
|	initial point in the direction v[j] (if j >= 0) or perform a curvilinear extrapolation
|	(if j < 0).
*/

double FLin (int j, double lambda)

{
	int		i;
	double	qa, qb, qc;

	if (j >= 0)
		{
		/* linear search */
		for (i = 0; i < n; i++)
			xnew[i] = gx[i] + lambda*vv[pos1(i,j,n)];
		}
	else
		{
		/* search along a parabolic space curve */
		qa = lambda*(lambda - qd1)/(qd0*(qd0 + qd1));
		qb = (lambda + qd0)*(qd1 - lambda)/(qd0*qd1);
		qc = lambda*(lambda + qd0)/(qd1*(qd0 + qd1));
		
		/* previous three points were stored as follows: x' in q0, x'' in gx, and x''' in q1;
		   see comments in 'Quad' */
		for (i = 0; i < n; i++)
			xnew[i] = qa*q0[i] + qb*gx[i] + qc*q1[i];
		}

	return (*fxn)(xnew);
}





/*--------------------------------------------------------------------------------------------------
|
|	Quad
|
|	Look for the function minimum along a curve defined by q0, q1, and x.
*/

void Quad (void)

{
	int		i;
	double	lambda, s;
	double	qa, qb, qc;
	
	/* q0 and q1 contain previous two points */

	s   = gfx;
	gfx = qf1;
	qf1 = s;
	qd1 = 0.0;
	for (i = 0; i < n; i++)
		{
		/* copy x to q1 for use in next cycle (but save current q1 in x so we can calculate the
		   norm below) */
		s     = gx[i];
		gx[i] = q1[i];
		q1[i] = s;			/* copy original x to q1 */
		qd1  += (q1[i] - gx[i])*(q1[i] - gx[i]);
		}
	qd1 = sqrt(qd1);
	if ((qd0 > 0.0) && (qd1 > 0.0) && (nl >= 3*n*n))
		{
		s = 0.0;
		lambda = qd1;
		LineMin(-1, 2, &s, &lambda, qf1, TRUE);
		qa = lambda*(lambda - qd1)/(qd0*(qd0 + qd1));
		qb = (lambda + qd0)*(qd1 - lambda)/(qd0*qd1);
		qc = lambda*(lambda + qd0)/(qd1*(qd0 + qd1));
		}
	else
		{
		gfx = qf1;
		qa  = qb = 0.0;
		qc  = 1.0;
		}
	qd0 = qd1;
	for (i = 0; i < n; i++)
		{
		s = q0[i];							/* save current q0 for calculation below */
		q0[i] = gx[i];						/* copy current q1 (now in gx) to q0 for next cycle */
		gx[i] = qa*s + qb*gx[i] + qc*q1[i];	/* gx now contains q1, and q1 contains gx */
		}
}





/*--------------------------------------------------------------------------------------------------
|
|	BracketMinimum
|
|	Bracket a function minimum using method of Press et al.
|
|	This function is reentrant.
*/

void BracketMinimum (double *pA, double *pB, MinimizeFxn func)

{
	double		a, b, c, x, fa, fb, fc, fx, xlim, r, q, temp;
	
	fa = (*func)(pA);
	a  = *pA;
	fb = (*func)(pB);
	b  = *pB;
	if (fb > fa)
		{
		/* swap a and b so we can go downhill in direction of a to b */
		temp = a;
		a    = b;
		b    = temp;
		temp = fa;
		fa   = fb;
		fb   = temp;
		}
	c  = b + GOLD*(b - a);		/* first guess for c */
	fc = (*func)(&c);

	while (fb > fc)
		{
		/* use inverse parabolic interpolation to compute a new trial point x (= the abcissa
		   which is the minimum of a parabola through f(a), f(b), and f(c))  */
		r = (b - a)*(fb - fc);
		q = (b - c)*(fb - fa);
		x = b - ((b - c)*q - (b - a)*r) / (2.0*SIGN(MAX(fabs(q-r), TINY), q-r));
		xlim = b + GLIMIT*(c - b);
		if ((b - x)*(x - c) > 0.0)
			{
			/* x is between b and c */
			fx = (*func)(&x);
			if (fx < fc)
				{
				/* there's a minimum between b and c */
				a = b;
#				if 0	/* lines below are apparently do-nothing statements */
				fa = fb;
				b  = x;
				fb = fx;
#				endif
				break;
				}
			else if (fx > fb)
				{
				/* there's a minimum between a and x */
				c = x;
#				if 0	/* line below is apparently a do-nothing statement */
				fc = fx;
#				endif
				break;
				}
			/* parabolic fit failed; "magnify" using default magnification */
			x  = c + GOLD*(c - b);
			fx = (*func)(&x);
			}
		else if ((c - x)*(x - xlim) > 0.0)
			{
			/* x from parabolic fit is between c and its allowed limit */
			fx = (*func)(&x);
			if (fx < fc)
				{
				b = c; fb = fc;
				c = x; fc = fx;
				x = c + GOLD*(c - b);
				fx = (*func)(&x);
				}
			}
		else if ((x - xlim)*(xlim - c) > 0.0)
			{
			/* limit parabolic fit to its maximum allowed value */
			x  = xlim;
			fx = (*func)(&x);
			}
		else
			{
			/* reject parabolic x; use default magnification */
			x  = c + GOLD*(c - b);
			fx = (*func)(&x);
			}
		
		/* new bracket is (b,c,x) */
		a = b; fa = fb;
		b = c; fb = fc;
		c = x; fc = fx;		
		}
	if (a < c)
		{
		*pA = a;
		*pB = c;
		}
	else
		{
		*pA = c;
		*pB = a;
		}

}





/*-------| rnum |------------------------------------------------
|   This pseudorandom number generator is described in:
|   Park, S. K. and K. W. Miller.  1988.  Random number generators: good
|      ones are hard to find.  Communications of the ACM, 31(10):1192-1201.
*/
#define A_A				16807L
#define M_M				2147483647L
#define Q_Q				127773L
#define R_R				2836L

double RandomNumber (long int *seed)

{

	long int	lo, hi, test;
	
	hi = (*seed) / Q_Q;
	lo = (*seed) % Q_Q;
	
	test = A_A * lo - R_R * hi;
	if (test > 0.0) 
		*seed = test;
	else	      
		*seed = test + M_M;
	
	return (double)(*seed) / (double)M_M;

}


/*--------------------------------------------------------------------------------------------------
|
|	SetIdentityMatrix2
|
|	Initialize a matrix to the identity matrix.
*/

void SetIdentityMatrix2 (double a[], int n)

{
	int			i, j;

	for (i = 0; i < n; i++)
		{
		for (j = 0; j < n; j++)
			a[pos1(i,j,n)] = 0.0;
		a[pos1(i,i,n)] = 1.0;
		}
}

