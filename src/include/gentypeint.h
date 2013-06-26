#ifndef __GENTYPEINT_H
#define __GENTYPEINT_H
// A generic residue-residue type interaction class.  This class is 
// intended to take a pair of residue types and residue coordinate and
// calculate an appropriate interaction potential.

#include <string>
#include <iostream>
#include "config.h"
#include "ljint.h"
#include "mtxintbase.h"
 
class CGenTypeInt : public CMtxIntBase {
 private:
  std::string   atype1,atype2;
  /* array of residue type x residue type interaction pointers */
  int           arraydim;
  CFFBase***    ff;
  /* residue type index in AAIntType array to use as reference state */
  int           refresindex;

 public:
  /* constructors */
  CGenTypeInt() 
    : atype1(""), atype2(""), arraydim(0), ff(0) {};
  CGenTypeInt(const std::string atomtype1, const std::string atomtype2, 
	      const std::string excludelist) {
    Init(atomtype1,atomtype2,excludelist);
  }

  /* Destructor */
  ~CGenTypeInt() {
    for(int i=0; i<arraydim; i++) {
      for(int j=i; j<arraydim; j++) {
        delete ff[i][j];
      }
      delete[] ff[i];
    }
    delete[] ff;
  }

  /* Functions */
  void Init(const std::string, const std::string, const std::string);
  void Display(int);
  bool Potential(CResidue*, const int, CResidue*, const int, 
		 const double&, double&);
  bool ResResPot(CResidue*, const int, CResidue*, const int, 
		 const double&, double&);
  int PairProb(CResidue*, const int, CResidue*, const int, double*&);
  double ResResPairProb(CResidue*, const int, CResidue*, const int);
  bool Param(CResidue*, const int, CResidue*, const int, double&);
			

  std::string ID() {return("GN");}

};

#endif
