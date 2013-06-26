/* A 3x3 matrix class */

#include <iostream>
#include "include/matrix3.h"

using std::ostream;
using std::istream;
using std::cout;
using std::endl;

/*-------------------------------------------------------------------------*/
/* Construct an Euler rotation matrix using the 3 Eulerian angles */
/* Requires:  phi,theta,psi -- Euler rotation angles */
/*------------------------------------------------------------------------*/
CMatrix3::CMatrix3(const double phi, const double theta, const double psi) {
  ConstructEulerMtx(phi,theta,psi);
}

/*-------------------------------------------------------------------------*/
/* Construct an Euler rotation matrix using the 3 Eulerian angles */
/* Requires:  phi,theta,psi -- Euler rotation angles */
/*------------------------------------------------------------------------*/
void CMatrix3::ConstructEulerMtx(const double phi, const double theta, 
				 const double psi) {
  double    sinphi, cosphi, sintheta, costheta, sinpsi, cospsi;

  ndim = 3;
  sinphi = sin(phi);
  cosphi = cos(phi);
  sintheta = sin(theta);
  costheta = cos(theta);
  sinpsi = sin(psi);
  cospsi = cos(psi);
  m[0][0] = cospsi*cosphi - costheta*sinphi*sinpsi;
  m[0][1] = cospsi*sinphi + costheta*cosphi*sinpsi;
  m[0][2] = sinpsi*sintheta;
  m[1][0] = -sinpsi*cosphi - costheta*sinphi*cospsi;
  m[1][1] = -sinpsi*sinphi + costheta*cosphi*cospsi;
  m[1][2] = cospsi*sintheta;
  m[2][0] = sintheta*sinphi;
  m[2][1] = -sintheta*cosphi;
  m[2][2] = costheta;
}


ostream& operator<<(ostream& os, const CMatrix3& matrix) {
  for(int j=0; j<matrix.ndim; j++) {
    os << matrix.m[0][j] << " " << matrix.m[1][j] 
       << " " << matrix.m[2][j] << endl;
  }
  return os;
}

istream& operator>>(istream& is, CMatrix3& matrix) {
  for(int j=0; j<matrix.ndim; j++) {
    is >> matrix.m[0][j] >> matrix.m[1][j] >> matrix.m[2][j];
  }
  return is;
}
