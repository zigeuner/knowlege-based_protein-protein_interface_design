#ifndef __MOVEBASE_H
#define __MOVEBASE_H
// This file provides an abstract parent class intended to
// serve as the base of more specific move-type classes

#include <cstdio>
//#include <cmath>
#include "config.h"
#include "vector3.h"
#include "matrix3.h"
#include "system.h"      /* has to be here at the base */
#include "forcefield.h"  /* has to be here at the base */

class CMoveBase {
 private:

 public:
  /* constructors */
  CMoveBase() {};
  CMoveBase(double) {};
    
  /* destructor */
  //virtual ~CMoveBase() = 0;
  virtual ~CMoveBase() { };

  /* display */
  virtual void Display() = 0;
  virtual std::string ID() = 0;

  /* Generic system move routines */
  virtual bool Move(CSystem&, const CForcefield&) = 0;
  virtual void ResizeStep(const double) = 0;
  virtual void EnforceMinStep(const double) = 0;
  virtual void EnforceMaxStep(const double) = 0;
  virtual void Reset() = 0;
  virtual std::string Disp() = 0;
  virtual void Disp(std::string&, std::string&) = 0;
};

#endif

