#ifndef __MTXINTBASE_H
#define __MTXINTBASE_H
// This file provides an abstract parent class intended to
// serve as the base of more specific classes that use matrices
// of "forcefield base" CFFBase classes to evaluate a set of
// interactions for a specific residue-residue pair.

#include <cstdio>
//#include <cmath>
#include "config.h"
#include "vector3.h"
#include "ffbase.h"
#include "residue.h"  /* has to be here rather than in pairint or ljint */

class CMtxIntBase {
 private:

 public:
  /* constructors */
  CMtxIntBase() {};
  //  CMtxIntBase(double) {};
    
  /* destructor */
  //virtual ~CMtxIntBase() = 0;
  virtual ~CMtxIntBase() { };

  /* display */
  virtual void Display(int) = 0;

  /* Generic Energy and Force evaluations */
  virtual bool Potential(CResidue*, const int, 
			 CResidue*, const int, const double&, double&) = 0;
  virtual bool ResResPot(CResidue*, const int, 
			 CResidue*, const int, const double&, double&) = 0;

  /* Routine for evaluating pair probabilities */
  virtual int PairProb(CResidue*, const int, 
		       CResidue*, const int, double*&) = 0;
  virtual double ResResPairProb(CResidue*, const int, CResidue*, const int) = 0;
  virtual bool Param(CResidue*, const int, CResidue*, const int, double&) = 0;

};

#endif
