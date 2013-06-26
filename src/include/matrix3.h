#ifndef __MATRIX3_H
#define __MATRIX3_H

#ifdef _USE_OLD_LIBRARIES
#include <math.h>
#else
#include <cmath>
#endif
#include <string>
#include "config.h"
#include "stdio.h"
#include "vector3.h"

class CMatrix3 {

protected:
  int    ndim;
  double m[3][3];

public:
  /* Constructors */
  CMatrix3(): ndim(3), m() {};
  CMatrix3(const double value) : ndim(3) {
    for(int i=0; i<ndim; i++) {
      for(int j=0; j<ndim; j++) {
	m[i][j] = value;
      }
    }
  };
  CMatrix3(const double, const double, const double);

  /* Functions declared in the .cpp file */
  void ConstructEulerMtx(const double, const double, const double);

  /* Operator overloads */
  CMatrix3& operator=(const CMatrix3& rh) {
    ndim = rh.ndim;
    for(int i=0; i<ndim; i++) {
      for(int j=0; j<ndim; j++) {
	m[i][j] = rh.m[i][j];
      }
    }
    return *this;
  }
  CMatrix3& operator=(const double value) {
    for(int i=0; i<ndim; i++) {
      for(int j=0; j<ndim; j++) {
	m[i][j] = value;
      }
    }
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& os, const CMatrix3& matrix);
  friend std::istream& operator>>(std::istream& is, CMatrix3& matrix);

  friend CMatrix3 operator-(const CMatrix3& a, const CMatrix3& b) { 
    CMatrix3  obj;
    for(int i=0; i<a.ndim; i++) {
      for(int j=0; j<a.ndim; j++) {
	obj.m[i][j] = a.m[i][j] - b.m[i][j];
      }
    }    
    return (obj);
  }

  friend CMatrix3 operator+(const CMatrix3& a, const CMatrix3& b) { 
    CMatrix3  obj;
    for(int i=0; i<a.ndim; i++) {
      for(int j=0; j<a.ndim; j++) {
	obj.m[i][j] = a.m[i][j] + b.m[i][j];
      }
    }    
    return (obj);
  }

  friend CMatrix3 operator*(const double a, const CMatrix3& b) { 
    CMatrix3  obj;
    for(int i=0; i<b.ndim; i++) {
      for(int j=0; j<b.ndim; j++) {
	obj.m[i][j] = a * b.m[i][j];
      }
    }    
    return (obj);
  }

  friend CMatrix3 operator*(const CMatrix3& a, const double b) { 
    CMatrix3  obj;
    for(int i=0; i<a.ndim; i++) {
      for(int j=0; j<a.ndim; j++) {
	obj.m[i][j] = b * a.m[i][j];
      }
    }    
    return (obj);
  }

  friend CVector3 operator*(const CMatrix3& a, const CVector3& b) { 
    double    x,y,z;
    x = b.X() * a.m[0][0] + b.Y() * a.m[0][1] + b.Z() * a.m[0][2];
    y = b.X() * a.m[1][0] + b.Y() * a.m[1][1] + b.Z() * a.m[1][2];
    z = b.X() * a.m[2][0] + b.Y() * a.m[2][1] + b.Z() * a.m[2][2];
    return CVector3(x,y,z);
  }

};

#endif

