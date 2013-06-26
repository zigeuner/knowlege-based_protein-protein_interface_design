// These are the prototype statements for the numerical utilities

#ifndef __UTILS_H
#define __UTILS_H

#define SCALE   1.0
#define FRACT   1.0e-6
#define SQR(x)  ((x)*(x))
#define TINY    1.0e-20
#define GOLD    1.618034
#define GLIMIT  100.0

#define SHFT(a,b,c,d)   (a)=(b);(b)=(c);(c)=(d);

#define ITMAXP  500
#define TOL     2.0e-4

static double dmax1, dmax2;
#define DMAX(a,b) ((dmax1=(a)) > (dmax2=(b)) ? dmax1 : dmax2)

// this is a HACK to avoid compiler warnings
//static double dmin1, dmin2; 
#define DMIN(a,b) ((dmin1=(a)) < (dmin2=(b)) ? dmin1 : dmin2)

#define SIGN(x,y) ((y)>=0.0 ? fabs(x) : -fabs(x))
#define CGOLD   0.3819660
#define ZEPS    1.0e-10
#define SQRTEPS 1.1e-6
#define ALF     1.0e-4
#define TOLX    1.0e-7


//int     ncom;
//double  *pcom, *xicom, (*nrfunc)(double []);
extern int     ncom;
extern double  *pcom, *xicom, (*nrfunc)(double []);

void powell(double p[], double **xi, int n, double ftol, int *iter,
                double *fret, double (*func)(double []));
void linmin(double p[], double xi[], int n, double *fret,
                double (*func)(double []));
void mnbrak(double *ax, double *bx, double *cx, double *fa,
               double *fb, double *fc, double (*func)(double));
double brent(double ax, double bx, double cx, double (*f)(double),
              double tol, double *xmin);
double f1dim(double x);

void gaussj(double **, int, double **, int);

void lubksb(double** a, int n, int* indx, double* b);

void ludcmp(double** a, int n, int* indx, double* d); 

void mrqmin(double [], double [], double [], int, double [], int [],
	    int, double **, double **, double *, void (*f)(double,
	    double[], double *, double [], int), double *);
void mrqcof(double [], double [], double sig[], int, double [], int[],
	    int, double **, double [], double *, void (*f)(double,
	    double[], double *, double [], int));
void covsrt(double **, int, int [], int);

double *vectorp(long nl, long nh);

void free_vector(double *v, long nl, long nh);

int Min(const int& a, const int& b);

double Min(const double& a, const double& b);

#endif
