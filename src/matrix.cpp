/* A generic mxm matrix class */

#include <iostream>
#include <fstream>
#include <string>
#include "include/matrix.h"

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::istream;
using std::string;
using std::cout;
using std::ios;
using std::endl;

ostream& operator<<(ostream& os, const CMatrix& matrix) {
  os << "matrix dimensions: " << matrix.ndim << "x" << matrix.ndim << endl;
  for(int i=0; i<matrix.ndim; i++) {
    for(int j=0; j<matrix.ndim; j++) {
      os << matrix.m[i][j] << " "; 
    }
    os << endl;
  }
  return os;
}

istream& operator>>(istream& is, CMatrix& matrix) {
  for(int i=0; i<matrix.ndim; i++) {
    for(int j=0; j<matrix.ndim; j++) {
      is >> matrix.m[i][j]; 
    }
  }
  return is;
}

/*--------------------------------------------------------------------*/
/* Symmeterize the matrix so that Mij = Mji */
/*--------------------------------------------------------------------*/
void CMatrix::Symmeterize() {
  double save;
  for(int i=0; i<ndim; i++) {
    for(int j=i+1; j<ndim; j++) {
      save = m[i][j];
      m[i][j] = (save + m[j][i])/2.0;
    }
  }
  for(int i=0; i<ndim; i++) {
    for(int j=0; j<i; j++) {
      m[i][j] = m[j][i];
    }
  }
}

/*--------------------------------------------------------------------*/
/* Normalize the upper triangle (including the diagonal) and then */
/* copy these values into the lower triangle  */
/*--------------------------------------------------------------------*/
void CMatrix::TriangleNormalize() {
  double sum=0;

  /* get the normalization factor */
  for(int i=0; i<ndim; i++) {
    for(int j=i; j<ndim; j++) {
      sum += m[i][j];
    }
  }
  /* normalize */
  for(int i=0; i<ndim; i++) {
    for(int j=i; j<ndim; j++) {
      m[i][j] /= sum;
    }
  }
  /* reflect */
  for(int i=0; i<ndim; i++) {
    for(int j=0; j<i; j++) {
      m[i][j] = m[j][i];
    }
  }
}

/*--------------------------------------------------------------------*/
/* A routine to dump a matrix to file */
/* Dumps the supplied comment line as the last comment before the */
/* matrix if the comment line isn't empty */
/* Requires:  filename -- filename to write into */
/*            commentline -- input comment line */
/*--------------------------------------------------------------------*/
void CMatrix::Write(const string filename, const string commentline) {
  ofstream outfile;

  /* Open the outfile for output*/
  outfile.open(filename.c_str(), ios::out);
  if(!outfile) {
    cout << endl << " ERROR: cannot open file '" << filename << "'" << endl;
    exit(-1);
  }
  printf("Writing matrix to filename '%s'\n",filename.c_str());

  outfile << "# matrix dimensions: " << ndim << "x" << ndim << endl;
  if (commentline.length() > 0) {
    outfile << commentline << endl;
  }
  for(int i=0; i<ndim; i++) {
    for(int j=0; j<ndim; j++) {
      outfile << m[i][j] << " "; 
    }
    outfile << endl;
  }
}

/*--------------------------------------------------------------------*/
/* A display routine to dump forcefield parameters to screen */
/* Requires:  indent -- number of spaces to indent */
/*--------------------------------------------------------------------*/
void CMatrix::Display(int indent) {
  int   i;
  char  spacing[100];

  for(i=0;i<indent;i++) {spacing[i] = ' ';}
  spacing[indent] = '\0';

  cout << spacing << "# matrix dimensions: " << ndim << "x" << ndim << endl;
  for(int i=0; i<ndim; i++) {
    cout << spacing; 
    for(int j=0; j<ndim; j++) {
      cout << m[i][j] << " "; 
    }
    cout << endl;
  }
}

