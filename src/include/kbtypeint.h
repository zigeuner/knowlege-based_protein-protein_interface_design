#ifndef __KBTYPEINT_H
#define __KBTYPEINT_H
// A residue-residue type interaction class based on knowledge-based 
// potential evaluations.  This class is intended to identify the pair of 
// residue types and calculate any appropriate interaction potentials

#include <string>
#include <iostream>
#include <cstdlib>
#include "config.h"
#include "pairint.h"
#include "mtxintbase.h"
#include "matrix.h"
 
class CKBTypeInt : public CMtxIntBase {
 private:
  std::string   atype1,atype2;
  /* Chi factors for residue type pairs */
  CMatrix       chimtx;
  /* array of residue type x residue type interaction pointers */
  int           arraydim;
  CFFBase***    ff;
  /* residue type index in AAIntType array to use as reference state */
  int           refresindex;

 public:
  /* constructors */
  CKBTypeInt() 
    : atype1(""), atype2(""), chimtx(0) , arraydim(0), ff(0) {};
  CKBTypeInt(const std::string potentialsdir, const std::string
	     atomtype1, const std::string atomtype2, int refindex,
	     const std::string excludelist, bool verbose) {
    Init(potentialsdir,atomtype1,atomtype2,refindex,excludelist,verbose);
  }

  /* Destructor */
  ~CKBTypeInt() {
    for(int i=0; i<arraydim; i++) {
      for(int j=i; j<arraydim; j++) {
        delete ff[i][j];
      }
      delete[] ff[i];
    }
    delete[] ff;
  }

  /* Functions */
  void Init(const std::string, const std::string, 
	    const std::string, int, const std::string, bool);
  void Display(int);
  bool Potential(CResidue*, const int, CResidue*, const int, 
		 const double&, double&);
  bool ResResPot(CResidue*, const int, CResidue*, const int, 
		 const double&, double&);
  int PairProb(CResidue*, const int, CResidue*, const int, double*&);
  double ResResPairProb(CResidue*, const int, CResidue*, const int);
  bool Param(CResidue*, const int, CResidue*, const int, double&);

};

#endif
