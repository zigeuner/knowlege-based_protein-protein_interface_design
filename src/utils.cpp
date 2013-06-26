// These are general routines for doing numerical manipulation.
// Also included here are general utility routines.

#ifdef _USE_OLD_LIBRARIES
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#else
#include <cstdlib>
#include <cmath>
#include <cstdio>
#endif
#include <iostream>
#include <fstream>
#include <string>
//#include "include/mrqmin.h"
//#include "include/gaussj.h"
//#include "include/powell.h"
#include "include/utils.h"

using std::ifstream;
using std::ofstream;
using std::ostream;
using std::string;
using std::cout;
using std::ios;
using std::endl;

int     ncom;
double  *pcom, *xicom, (*nrfunc)(double []);

/* Moved these here from utils.h because they were causing linker */
/* problems. */
double *vectorp(long nl, long nh)
{
  double *v;

  v = (double *)malloc((size_t)((nh-nl+2)*sizeof(double)));
  if (!v) printf("allocation error in vector\n");
  return v-nl+1;
}

void free_vector(double *v, long nl, long nh)
{
  free((char *)(v+nl-1));
}

int Min(const int& a, const int& b){
  if (a<b) {return(a);}
  return(b);
}

double Min(const double& a, const double& b){
  if (a<b) {return(a);}
  return(b);
}

/*
 *  Levenberg-Marquardt method, attempting to reduce the value chi^2 of
 *  a fit between a set of data points x[0..ndata-1], y[0..ndata-1] with
 *  individual standard deviations sig[0..ndata-1], and a nonlinear
 *  function dependent on ma coefficients a[0..ma-1].  The input array
 *  ia[0..ma-1] indicates by nonzero entries those components of a that
 *  should be fitted for, and by zero entries those components that
 *  should be held fixed at their input values.  The program returns
 *  current best-fit values for the parameters a[0..ma-1], and chi^2.
 *  The arrays covar[0..ma-1][0..ma-1], alpha[0..ma-1][0..ma-1] are used
 *  as working space during most iterations.  Supply a routine 
 *  funcs(x,a,yfit,dyda,ma) that evaluates the fitting function yfit, and
 *  its derivatives dyda[0..ma-1] with respect to the fitting parameters
 *  a at x.  On the first call provide an initial guess for the
 *  parameters a, and set alamda < 0 for initialization (which then sets
 *  alamda = 0.001).  If a step succeeds chi^2 becomes smaller and
 *  alamda decreases by a factor of 10.  If a step fails, alamda grows by
 *  a factor of 10.  You must call this routine repeatedly until convergence
 *  is achieved.  Then, make one final call with alamda=0, so that
 *  covar[0..ma-1][0..ma-1] returns the covariance matrix, and alpha the
 *  curvature matrix.  (Parameters held fixed will return zero
 *  covariances). (NRC)
 */

#define SWAP(a,b) {swap=(a); (a)=(b); (b)=swap;}

void mrqmin(double x[], double y[], double sig[], int ndata, double a[],
	    int ia[], int ma, double **covar, double **alpha, double *chisq,
	    void (*funcs)(double,double[],double*,double[],int), 
            double *alamda)
{
  static int	mfit;
  static double	ochisq, *atry, *beta, *da, **oneda;
  int		j, k, l, m;

  if( *alamda < 0.0 ) {
    atry = (double*) malloc(ma * sizeof(double));
    beta = (double*) malloc(ma * sizeof(double));
    da = (double*) malloc(ma * sizeof(double));
    for(mfit=0, j=0; j<ma; j++)
      if( ia[j] ) mfit++;
    oneda = (double**) malloc(mfit * sizeof(double *));
    for(m=0; m<mfit; m++) oneda[m] = (double*) malloc(sizeof(double));
    *alamda = 0.001;
    mrqcof(x, y, sig, ndata, a, ia, ma, alpha, beta, chisq, funcs);
    ochisq = (*chisq);
    for(j=0; j<ma; j++) atry[j] = a[j];
  }
  for(j=0; j<mfit; j++) {
    for(k=0; k<mfit; k++) covar[j][k] = alpha[j][k];
    covar[j][j] = alpha[j][j] * (1.0 + (*alamda));
    oneda[j][0] = beta[j];
  }

  gaussj(covar, mfit, oneda, 1);
  for(j=0; j<mfit; j++) da[j] = oneda[j][0];

  if( *alamda == 0.0 ) {
    covsrt(covar, ma, ia, mfit);
    covsrt(alpha, ma, ia, mfit);
    for(m=0; m<mfit; m++)
      free(oneda[m]);
    free(oneda);
    free(da);
    free(beta);
    free(atry);
    return;
  }

  for(j=-1, l=0; l<ma; l++)
    if (ia[l]) atry[l] = a[l] + da[++j];
  mrqcof(x, y, sig, ndata, atry, ia, ma, covar, da, chisq, funcs);
  if( *chisq < ochisq ) {
    *alamda *= 0.1;
    ochisq = (*chisq);
    for(j=0; j<mfit; j++) {
      for(k=0; k<mfit; k++) alpha[j][k] = covar[j][k];
      beta[j] = da[j];
    }
    for(l=0; l<ma; l++) a[l] = atry[l];
  } else {
    *alamda *= 10.0;
    *chisq = ochisq;
  }
}


/*
 *  Used by mrqmin to evaluate the linearized fitting matrix alpha, and
 *  vector beta as in (15.5.8) and calculate chi^2. (NRC)
 */

void mrqcof(double x[], double y[], double sig[], int ndata, double a[],
	    int ia[], int ma, double **alpha, double beta[], double *chisq,
	    void (*funcs)(double, double[], double *, double [], int))
{
  int		i, j, k, l, m, mfit=0;
  double	ymod, wt, sig2i, dy, *dyda;

  dyda = (double*) malloc(ma * sizeof(double));
  for(j=0; j<ma; j++)
    if (ia[j]) mfit++;
  for(j=0; j<mfit; j++) {
    for(k=0; k<=j; k++) alpha[j][k] = 0.0;
    beta[j] = 0.0;
  }
  *chisq = 0.0;
  for(i=0; i<ndata; i++) {
    (*funcs)(x[i], a, &ymod, dyda, ma);
    sig2i = 1.0 / (sig[i] * sig[i]);
    dy = y[i] - ymod;
    for(j=-1, l=0; l<ma; l++) {
      if (ia[l]) {
	wt = dyda[l] * sig2i;
	for(j++, k=-1, m=0; m<=l; m++)
	  if (ia[m]) alpha[j][++k] += wt*dyda[m];
	beta[j] += dy * wt;
      }
    }
    *chisq += dy * dy * sig2i;
  }
  for(j=1; j<mfit; j++)
    for(k=0; k<j; k++) alpha[k][j] = alpha[j][k];
  free(dyda);
}


/*
 *  Expand in storage the covariance matrix covar, so as to take into
 *  account parameters that are being held fixed. (For the latter, return
 *  zero covariances.) (NRC)
 */

void covsrt(double **covar, int ma, int ia[], int mfit)
{
  int		i, j, k;
  double	swap;

  for(i=mfit; i<ma; i++)
    for(j=0; j<=i; j++) covar[i][j] = covar[j][i] = 0.0;
  k = mfit-1;
  for(j=ma-1; j>=0; j--) {
    if (ia[j]) {
      for(i=0; i<ma; i++) SWAP(covar[i][k], covar[i][j]);
      for(i=0; i<ma; i++) SWAP(covar[k][i], covar[j][i]);
      k--;
    }
  }
}


/* Powell's method for finding minimum of multi-dimensional function */

/*	
 *  Minimization of a function func of n variables.  Input consists
 *  of an initial starting point p[1..n]; an initial matrix 
 *  xi[1..n][1..n], whose columns contain the initial set of directions
 *  (usually the n unit vectors); and ftol, the fractional tolerance
 *  in the function value such that failure to decrease by more than
 *  this amount on one iteration signals doneness. On output,
 *  p is set to the best point found, xi is the then-current direction
 *  set, fret is the returned function value at p, and iter is the
 *  number of iterations taken.  The routine linmin is used.  (NRC)
 */

void powell(double p[], double **xi, int n, double ftol, int *iter,
            double *fret, double (*func)(double []))
{
  int   i, ibig, j;
  double del, fp, fptt, t, *pt, *ptt, *xit;

  //pt = vectorp(1,n);
  //ptt= vectorp(1,n);
  //xit= vectorp(1,n);
  pt = (double*)malloc(n*sizeof(double));
  ptt= (double*)malloc(n*sizeof(double));
  xit= (double*)malloc(n*sizeof(double));
  *fret = (*func)(p);
  //for(j=1;j<=n;j++) pt[j]=p[j];
  for(j=0;j<n;j++) pt[j]=p[j];
  for(*iter=1;;++(*iter)) {
    fp=(*fret);
    ibig=0;
    del=0.0;
    //for(i=1;i<=n;i++) {
    for(i=0;i<n;i++) {
      //for(j=1;j<=n;j++) xit[j]=xi[j][i];
      for(j=0;j<n;j++) xit[j]=xi[j][i];
      fptt=(*fret);
      linmin(p,xit,n,fret,func);
      if (fabs(fptt-(*fret)) > del) {
        del = fabs(fptt-(*fret));
        ibig = i;
      }
    }
    if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
      //free_vector(xit,1,n);
      //free_vector(ptt,1,n);
      //free_vector(pt,1,n);
      free(xit);
      free(ptt);
      free(pt);
      return;
    }
    if (*iter==ITMAXP) printf("Maximum iterations exceeded\n");
    //for(j=1;j<=n;j++) {
    for(j=0;j<n;j++) {
      ptt[j]=2.0*p[j]-pt[j];
      xit[j]=p[j]-pt[j];
      pt[j]=p[j];
    }
    fptt=(*func)(ptt);
    if (fptt < fp) {
      t = 2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
      if (t<0.0) {
        linmin(p,xit,n,fret,func);
        //for(j=1;j<=n;j++) {
        for(j=0;j<n;j++) {
          xi[j][ibig]=xi[j][n];
          xi[j][n]=xit[j];
        }
      }
    }
  }
}

void linmin(double p[], double xi[], int n, double *fret,
            double (*func)(double []))
{
  int   j;
  double xx, xmin, fx, fb, fa, bx, ax;

  ncom=n;
  //pcom=vectorp(1,n);
  //xicom=vectorp(1,n);
  pcom=(double*)malloc(n*sizeof(double));
  xicom=(double*)malloc(n*sizeof(double));
  nrfunc=func;
  //for(j=1;j<=n;j++) {
  for(j=0;j<n;j++) {
    pcom[j]=p[j];
    xicom[j]=xi[j];
  }
  ax = 0.0;
  xx = 1.0;
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
  *fret = brent(ax,xx,bx,f1dim,TOL,&xmin);
  //for(j=1;j<=n;j++) {
  for(j=0;j<n;j++) {
    xi[j] *= xmin;
    p[j] += xi[j];
  }
  //free_vector(xicom,1,n);
  //free_vector(pcom,1,n);
  free(xicom);
  free(pcom);
}

double f1dim(double x)
{
  int   j;
  double f, *xt;

  //xt = vectorp(1,ncom);
  xt = (double*)malloc(ncom*sizeof(double));;
  //for(j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
  for(j=0;j<ncom;j++) xt[j]=pcom[j]+x*xicom[j];
  f=(*nrfunc)(xt);
  //free_vector(xt,1,ncom);
  free(xt);
  return f;
}

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
            double *fc, double (*func)(double))
{
  double        ulim, u, r, q, fu, dum;

  *fa=(*func)(*ax);
  *fb=(*func)(*bx);
  if (*fb > *fa) {
    SHFT(dum,*ax,*bx,dum)
    SHFT(dum,*fb,*fa,dum)
  }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=(*func)(*cx);
  while (*fb>*fc) {
    r = (*bx-*ax)*(*fb-*fc);
    q = (*bx-*cx)*(*fb-*fa);
    u = (*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
        (2.0*SIGN(DMAX(fabs(q-r),TINY),q-r));
    ulim = (*bx)+GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu=(*func)(u);
      if (fu < *fc) {
        *ax = (*bx);
        *bx = u;
        *fa = (*fb);
        *fb = fu;
        return;
      } else if (fu > *fb) {
        *cx = u;
        *fc = fu;
        return;
      }
      u = (*cx)+GOLD*(*cx-*bx);
      fu = (*func)(u);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu = (*func)(u);
      if (fu<*fc) {
        SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
        SHFT(*fb,*fc,fu,(*func)(u))
      }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u = ulim;
      fu = (*func)(u);
    } else {
      u = (*cx)+GOLD*(*cx-*bx);
      fu = (*func)(u);
    }
    SHFT(*ax,*bx,*cx,u)
    SHFT(*fa,*fb,*fc,fu)
  }
}

double brent(double ax, double bx, double cx, double (*f)(double),
             double tol, double *xmin)
{
  int   iter;
  double a, b, d, etemp, fu, fv, fw, fx, p, q;
  double r, tol1, tol2, u, v, w, x, xm;
  double e=0.0;

  a = (ax < cx ? ax : cx);
  b = (ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x);
  for(iter=1;iter<=ITMAXP;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin = x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r = (x-w)*(fx-fv);
      q = (x-v)*(fx-fw);
      p = (x-v)*q-(x-w)*r;
      q = 2.0*(q-r);
      if (q>0.0) p=-p;
      q = fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) ||
          p >= q*(b-x))
        d = CGOLD*(e=(x>=xm?a-x:b-x));
      else {
        d = p/q;
        u = x+d;
        if (u-a < tol2 || b-u < tol2)
          d = SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e= (x>=xm ? a-x : b-x));
    }
    u = (fabs(d)>=tol1 ? x+d : x+SIGN(tol1,d));
    fu = (*f)(u);
    if (fu<=fx) {
      if (u>=x) a=x; else b=x;
      SHFT(v,w,x,u)
      SHFT(fv,fw,fx,fu)
    } else {
      if (u < x) a=u; else b=u;
      if (fu<=fw || w==x) {
        v = w;
        w = u;
        fv = fw;
        fw = fu;
      } else if (fu<=fv || v==x || v==w) {
        v = u;
        fv = fu;
      }
    }
  }
  printf("Too many iterations in brent\n");
  *xmin = x;
  return fx;
}


/*
 *  Linear equation solution by Gauss-Jordan elimination, equation (2.1.1)
 *  above.  a[0..n-1][0..n-1] is the input matrix.  b[0..n-1][0..m-1] is
 *  input containing the m right-hand side vectors.  On output, a is
 *  replaced by its matrix inverse, and b is replaced by the corresponding
 *  set of solution vectors. (NRC) 
 */ 

#define SWAPG(a,b)	{temp=(a);(a)=(b);(b)=temp;}

void gaussj(double **a, int n, double **b, int m)
{
  int		*indxc, *indxr, *ipiv;
  int		i, icol = 0, irow = 0, j, k, l, ll;
  double	big, dum, pivinv, temp;

  indxc = (int*) malloc(n * sizeof(int));
  indxr = (int*) malloc(n * sizeof(int));
  ipiv = (int*) malloc(n * sizeof(int));
 
  for(j=0; j<n; j++) ipiv[j] = 0;

  for(i=0; i<n; i++) {
    big = 0.0; 
    for(j=0; j<n; j++)
      if (ipiv[j] != 1)
	for(k=0; k<n; k++) {
	  if (ipiv[k] == 0) {
	    if (fabs(a[j][k]) >= big) {
	      big = fabs(a[j][k]);
	      irow = j;
	      icol = k;
	    }
	  } else
	    if (ipiv[k] > 1) 
	      fprintf(stderr,"gaussj: Singular Matrix (1)\n");
	}
    ++(ipiv[icol]);
    if (irow != icol) {
      for(l=0; l<n; l++) SWAPG(a[irow][l],a[icol][l])
      for(l=0; l<m; l++) SWAPG(b[irow][l],b[icol][l])
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if (a[icol][icol] == 0.0) {
      fprintf(stderr,"gaussj: Singular Matrix (2)\n");
      exit(1);
    }
    pivinv = 1.0 / a[icol][icol];
    a[icol][icol] = 1.0;
    for(l=0; l<n; l++) a[icol][l] *= pivinv;
    for(l=0; l<m; l++) b[icol][l] *= pivinv;
    for(ll=0; ll<n; ll++)
      if (ll != icol) {
        dum = a[ll][icol];
	a[ll][icol] = 0.0;
	for(l=0; l<n; l++) a[ll][l] -= a[icol][l]*dum;
	for(l=0; l<m; l++) b[ll][l] -= b[icol][l]*dum;
      }
  }

  for(l=n-1; l>=0; l--) {
    if (indxr[l] != indxc[l])
      for(k=0; k<n; k++) 
        SWAPG(a[k][indxr[l]],a[k][indxc[l]]);
  }

  free(ipiv);
  free(indxr);
  free(indxc);
}

void lubksb (double** a, int n, int* indx, double* b)
{
  /*  Solves the set of n linear equations A * X = B.  Here a[1..n][1..n] 
      is input, not as the matrix A but rather as its LU decomposition, 
      determined by the routine ludcmp.  indx[1..n] is input as the 
      permutation vector returned by ludcmp.  b[1..n] is input as the 
      right-hand side vector B, and returns with the solution vector X. 
      a, n, and indx are not modifed by this routine and can be left in 
      place for successive calls with different right-hand sides b. This 
      routine takes into account the possibility that b will begin with 
      many zero elements,  so it is efficient for use in matrix inversion.  
      */ 

  int   i, ii = 0, ip, j;
  double sum;

  for (i=0; i<n; i++) {         // When ii is set to a positive value, 
    ip = indx[i];               // it will become the index of the first 
    sum = b[ip];                // nonvanishing element of b. We now do 
    b[ip] = b[i];               // the forward substitution, equation 
    if( ii ) {                  // (2.3.6). The only new wrinkle is to 
      for(j=ii; j<i; j++)       // unscramble the permutation as we go.
        sum -= a[i][j]*b[j];
    }
    else if( sum )          // A nonzero element was encountered, so 
      ii = i;               // from now on we will have to do the sums 
    b[i] = sum;             // in the loop above. 
  }
  for(i=n-1; i>=0; i--) {   // Now we do the backsubstitution, equation 
    sum = b[i];             // (2.3.7).
    for(j=i+1; j<n; j++) 
      sum -= a[i][j]*b[j];
    b[i] = sum/a[i][i];     // Store a component of the solution vector X.
  }
}


void ludcmp (double** a, int n, int* indx, double* d)
{
  /*  Given a matrix a[1..n][1..n], this routine replaces it by the LU 
      decomposition of a rowwise permutation of itself.  a and n are input.  
      a is output, arranged as in equation (2.3.14) above; indx[1..n] is 
      an output vector that records the row permutation effected by the 
      partial pivoting; d is output as +/-1 depending on whether the number 
      of row interchanges was even or odd, respectively. This routine is 
      used in combination with lubksb to solve linear equations or invert a 
      matrix.  
      */
  
  int     i, imax = 0, j, k;
  double  big, dum, sum, temp;
  double* vv; // vv stores the implicit scaling of each row.

#ifdef DEBUG
  cout << "DEBUG at " << __FILE__ << ":" << __LINE__ << endl;
  cout << " matrix size: " << n << endl;
  for (i=0; i<n; i++) { 
    for (j=0; j<n; j++) { 
      cout << a[i][j] << " ";
    }
    cout << endl;
  }
#endif
  
  vv = new double[n];
  
  *d = 1.0;   // No row interchanges yet.
  for (i=0; i<n; i++) {     // Loop over rows to get the implicit 
    big = 0.0;            // scaling information. 
    for (j=0; j<n; j++)
      if ( (temp=fabs(a[i][j])) > big ) big = temp;
    if (big == 0.0) {     // No nonzero largest element.
      cout << "Singular matrix in routine ludcmp." << endl;
      exit(1);
    } 
    vv[i] = 1.0/big;      // Save the scaling.
  }
  for (j=0; j<n; j++) {     // This is the loop over columns of 
    // Crout's method.
    for (i=0; i<j; i++) { // This is equation (2.3.12) except 
      sum = a[i][j];    // for i = j.
      for (k=0; k<i; k++) 
        sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;            // Initialize for the search for largest 
    // pivot element.
    for (i=j; i<n; i++) { // This is i = j of equation (2.3.12) and 
      // i = j + 1...N of equation (2.3.13). 
      sum = a[i][j];
      for (k=0; k<j; k++)
        sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      if( (dum=vv[i]*fabs(sum)) >= big ) { // Is the figure of merit 
        big = dum;                       // for the pivot better 
        imax = i;                        // than the best so far? 
      }
    }
    if (j != imax) { // Do we need to interchange rows?
      for (k=0; k<n; k++) { // Yes, do so...
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      *d = -(*d);           // ...and change the parity of d.
      vv[imax] = vv[j];     // Also interchange the scale factor.
    }
    indx[j] = imax;
    if( a[j][j] == 0.0 )  // If the pivot element is zero the matrix 
      a[j][j] = 1E-20;  // is singular (at least to the precision of 
    // the algorithm). For some applications on 
    // singular matrices, it is desirable to 
    // substitute 1E-20 for zero.
    if( j != n ) {        // Now, finally, divide by the pivot element.
      dum = 1.0/(a[j][j]);
      for(i=j+1; i<n; i++) 
        a[i][j] *= dum;
    }
  }  // Go back for the next column in the reduction.
  delete [] vv;
}

/* 
 * Given an n-dimensional point xold[1..n], the value of the function and 
 * gradient there, fold and g[1..n], and a direction p[1..n], finds a new 
 * point x[1..n] along the direction p from xold where the function func has 
 * decreased sufficiently.  The new function value is returned in f. 
 * stpmax is an input quantity that limits the length of the steps so that 
 * you do not try to evaluate the function in regions where it is undefined 
 * or subject to overflow.  p is usually the Newton direction.  The output 
 * quantity check is false (0) on a normal exit.  It is true (1) when x is 
 * too close to xold.  In a minimization algorithm, this usually signals 
 * convergence and can be ignored.  However, in a zero-finding algorithm the 
 * calling program should check whether the convergence is spurious.  Some 
 * "difficult" problems may require double precision in this routine. 
 */

#define ALF 1.0e-4  // Ensures sufficient decrease in function value. 
#define TOLX 1.0e-7 // Convergence criterion on delx . 

void lnsrch(int n, float xold[], float fold, float g[], float p[], float x[], 
            float *f, float stpmax, int *check, float (*func)(float [])) 
{ 
  double a, alam, alam2, alamin, b, disc, f2, rhs1, rhs2, slope, sum, temp;
  double test,tmplam; 

  *check=0; 

  sum = 0.0;
  for(int i=0; i<n; i++) 
    sum += p[i]*p[i];
  sum = sqrt(sum);

  if( sum > stpmax ) 
    for(int i=0; i<n; i++) 
      p[i] *= stpmax/sum; // Scale if attempted step is too big

  slope = 0.0;
  for(int i=0; i<n; i++) 
    slope += g[i]*p[i]; 

  if( slope >= 0.0 ) {
    cout << " Error: roundoff problem in lnsrch" << endl << endl;
    exit(1);
  }

  test = 0.0; // Compute lambda_min
  for(int i=0; i<n; i++) {
    temp = fabs(p[i])/DMAX(fabs(xold[i]),1.0); 
    if( temp > test ) test = temp; 
  } 
  alamin = TOLX/test; 

  alam = 1.0;  // Always try full Newton step first
  for( ; ; ) { // Start of iteration loop
    for(int i=0; i<n; i++) x[i] = xold[i] + alam*p[i]; 

    *f = (*func)(x); 

    if( alam < alamin ) { // Convergence on delx. For zero finding, the 
                          // calling program should verify the convergence
      for(int i=0; i<n; i++) x[i] = xold[i]; 
      *check = 1; 
      return; 
    } 
    else if( *f <= fold+ALF*alam*slope ) return; // Sufficient decrease 
    else { // Backtrack
      if( alam == 1.0 ) 
	tmplam = -slope / ( 2.0*(*f - fold - slope) ); // First time
      else { // Subsequent backtracks
	rhs1 = *f - fold - alam*slope; 
	rhs2 = f2 - fold - alam2*slope; 

	a = (rhs1/(alam*alam) - rhs2/(alam2*alam2))/(alam-alam2);
	b = (-alam2*rhs1/(alam*alam) + alam*rhs2/(alam2*alam2))/(alam-alam2);

	if( a == 0.0 ) tmplam = -slope/(2.0*b); 
	else { 
	  disc = b*b - 3.0*a*slope; 
	  if( disc < 0.0 ) 
	    tmplam = 0.5*alam; 
	  else if( b <= 0.0 ) 
	    tmplam = (-b + sqrt(disc))/(3.0*a); 
	  else 
            tmplam = -slope/(b + sqrt(disc)); 
	} 
	if( tmplam > 0.5*alam ) 
	  tmplam = 0.5*alam; // lamda <= 0.5*lambda_1 
      } 
    } 
    alam2 = alam; 
    f2 = *f; 
    alam = DMAX(tmplam, 0.1*alam); // lambda >= 0.1*lambda_1 
  } // Try again. 
}

