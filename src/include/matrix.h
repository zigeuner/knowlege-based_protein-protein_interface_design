#ifndef __MATRIX_H
#define __MATRIX_H

#ifdef _USE_OLD_LIBRARIES
#include <math.h>
#else
#include <cmath>
#endif
#include <string>
#include "config.h"
#include "stdio.h"

class CMatrix {

protected:
  int       ndim;
  double**  m;

public:
  /* Constructors */
  CMatrix(): ndim(0), m() {};
  CMatrix(const int value) {
    Init(value);
  };

  /* Destructor */
  ~CMatrix() {
    for(int i=0; i<ndim; i++) {
      delete[] m[i];
    }
    delete[] m;
    ndim = 0;
  }

  /* Functions for returning the state variables of the object*/
  int Size() {return (ndim);}

  /* size the matrix */
  void Init(const int value) {
    ndim = value;
    m = new double*[value];
    for(int i=0; i<ndim; i++) {
      m[i] = new double[value];
      for(int j=0; j<ndim; j++) {
	m[i][j] = 0.0;
      }
    }
    //    std::cout << "init dimensions: " << ndim << "x" << ndim << std::endl;
  };

  /* Operator overloads */
  CMatrix& operator=(const CMatrix& rh) {
    ndim = rh.ndim;
    for(int i=0; i<ndim; i++) {
      for(int j=0; j<ndim; j++) {
	m[i][j] = rh.m[i][j];
      }
    }
    return *this;
  }
  CMatrix& operator=(const double value) {
    for(int i=0; i<ndim; i++) {
      for(int j=0; j<ndim; j++) {
	m[i][j] = value;
      }
    }
    return *this;
  }

  friend CMatrix operator-(const CMatrix& a, const CMatrix& b) { 
    CMatrix  obj;
    for(int i=0; i<a.ndim; i++) {
      for(int j=0; j<a.ndim; j++) {
	obj.m[i][j] = a.m[i][j] - b.m[i][j];
      }
    }    
    return (obj);
  }

  friend CMatrix operator+(const CMatrix& a, const CMatrix& b) { 
    CMatrix  obj;
    for(int i=0; i<a.ndim; i++) {
      for(int j=0; j<a.ndim; j++) {
	obj.m[i][j] = a.m[i][j] + b.m[i][j];
      }
    }    
    return (obj);
  }

  friend CMatrix operator*(const double a, const CMatrix& b) { 
    CMatrix  obj;
    for(int i=0; i<b.ndim; i++) {
      for(int j=0; j<b.ndim; j++) {
	obj.m[i][j] = a * b.m[i][j];
      }
    }    
    return (obj);
  }

  friend CMatrix operator*(const CMatrix& a, const double b) { 
    CMatrix  obj;
    for(int i=0; i<a.ndim; i++) {
      for(int j=0; j<a.ndim; j++) {
	obj.m[i][j] = b * a.m[i][j];
      }
    }    
    return (obj);
  }

  /* set an entry of the matrix */
  void Set(const int i, const int j, const double value) {
    m[i][j] = value;
  };

  /* add to a matrix entry */
  void Add(const int i, const int j, const double value) {
    m[i][j] += value;
  };

  /* Divide a matrix entry by a constant */
  void Divide(const int i, const int j, const double value) {
    m[i][j] /= value;
  };

  /* Multiply a matrix entry by a constant */
  void Multiply(const int i, const int j, const double value) {
    m[i][j] *= value;
  };

  /* get an entry of the matrix */
  double Entry(const int i, const int j) {
    return(m[i][j]);
  };

  friend std::ostream& operator<<(std::ostream& os, const CMatrix& matrix);
  friend std::istream& operator>>(std::istream& is, CMatrix& matrix);

  /* Functions defined in .cpp file */
  void Symmeterize();
  void TriangleNormalize();
  void Write(const std::string, const std::string);
  void Display(int);

};

#endif

