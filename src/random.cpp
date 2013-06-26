// Routines for picking random numbers and numbers from distributions

#ifdef _USE_OLD_LIBRARIES
#include <math.h>
#else
#include <cmath>
#endif
#include <iostream>
#include "include/config.h"
#include "include/random.h"

using std::cout;
using std::endl;

/* Reference?!?  */
double ran2(long& idum) {
  int         j;
  static long idum2 = 123456789, iy = 0, iv[NTAB];
  long        k;
  double      temp;

  if (idum <= 0) {
    if (-(idum) < 1) idum = 1;
    else idum = -(idum);
    idum2 = (idum);
    for (j=NTAB+7; j>=0; j--) {
      k = (idum)/IQ1;
      idum = IA1*(idum - k*IQ1) - k*IR1;
      if (idum < 0) idum += IM1;
      if (j < NTAB) iv[j] = idum;
    }
    iy = iv[0];
  }
  k = (idum)/IQ1;
  idum = IA1*(idum - k*IQ1) - k*IR1;
  if (idum < 0) idum += IM1;
  k = idum2/IQ2;
  idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j = iy/NDIV;
  iy = iv[j] - idum2;
  iv[j] = idum;
  if (iy < 1) iy += IMM1;
  if ((temp = AM*iy) > RNMX) return RNMX;
  else return temp;
}

/* Reference?!?  */
/* Returns a normally distributed deviate with zero mean and unit variance */
/* using ran2(idum) as the source of uniform deviates. (NRC)               */ 
double gasdev(long& idum) {
  static int iset=0;
  static double gset;
  double fac, rsq, v1, v2;

  if(idum < 0) iset=0;                   /* Reinitialize. */
  if(iset == 0) {                         // We don't have an extra deviate
    // handy, so pick two uniform numbers
    do {                                // in the square extending from -1
      v1=2.0*ran2(idum) - 1.0;        // to +1 in each direction, see if
      // they are in the unit circle, and 
      v2=2.0*ran2(idum) - 1.0;        // if they are not, try again.
                                            
      rsq=v1*v1+v2*v2;

    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);        // Now make the Box-Muller
    // transformation to get two normal
    gset=v1*fac;                        // deviates.  Return one and save 
    // the other for next time.

    iset=1;                             /* Set flag. */ 
    return v2*fac;
  }
  else {                                  // We have an extra deviate handy,
    iset=0;                             // so unset the flag, and return it.
    return gset;
  }
}
