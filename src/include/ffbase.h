#ifndef __FFBASE_H
#define __FFBASE_H
// This file provides an abstract parent class intended to
// serve as the base of more specific energy and force evaluation 
// classes.

#include <cstdio>
//#include <cmath>
#include "config.h"
#include "vector3.h"
//#include "system.h"  /* has to be here rather than in pairint or ljint */

class CFFBase {
 private:

 public:
  /* constructors */
  CFFBase() {};
  CFFBase(double) {};
    
  /* destructor */
  //virtual ~CFFBase() = 0;
  virtual ~CFFBase() { };

  /* display */
  virtual void Display(int) = 0;
  virtual std::string ID() = 0;
  virtual void Check(const double, int&, int&, bool) = 0;

  /* Generic Energy and Force evaluations */
  virtual int Potential(const CVector3&, double&) = 0;
  virtual bool Force(const CVector3&, double, CVector3&) = 0;

};

#endif
