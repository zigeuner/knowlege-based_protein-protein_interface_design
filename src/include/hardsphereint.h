// This file provides the core routines for performing hard sphere
// energy and force evaluations.  Routines here are inline
// to allow the compiler to put them into the loop in forcefield.

#include <cmath>
#include <string>
#include "config.h"
#include "ffbase.h"
#include <stdlib.h>

class CHardSphere : public CFFBase {
 private:
  double       epsilon, radius, radius2;
  std::string  id;

 public:
  /* Constructors */
  CHardSphere() 
    : epsilon(0.0), radius(0.0), radius2(0.0), id("") {};
  CHardSphere(double epsilon, double radius, std::string label) {
    Init(epsilon, radius, label);
  }
  
  /* Destructor */
  ~CHardSphere() { /* cout << " CHardSphere destructor.." << endl; */ }

  /* Functions declared in .cpp file */
  void Init(double, double, std::string);
  void Init(double, std::string);
  void Check(const double, int&, int&, bool);
  void Display(int);

  std::string ID() {return(id);}

  /*------------------------------------------------------------------------*/
  /* Energy routine */
  /* Requires:  sepvec -- Rij separation vector */
  /*            potential -- returned potential */
  /* Returns: -1 -- failed (can't currently happen) */
  /*           0 -- under range (can't currently happen) */
  /*           1 -- ok, within range, returns epsilon */
  /*           2 -- above range, returns zero */
  /*------------------------------------------------------------------------*/
  int Potential(const CVector3& sepvec, double& potential) {
    double r2 = sepvec.SqLength();

    if (r2 > radius2) {
      potential = 0.0;      
      return(2);
    } else {
      potential = epsilon;
      return(1);
    }
  }

  /*------------------------------------------------------------------------*/
  /* Force routine */
  /* Requires:  sepvec -- Rij separation vector */
  /*            potential -- returned potential */
  /*            force -- returned force vector */
  /*------------------------------------------------------------------------*/
  bool Force(const CVector3& sepvec, double potential, CVector3& force) {
    std::cout << "Force routine not yet ready" << std::endl;
    exit(-1);    
    return(false);
  }
};

