// This file provides the core routines for performing Lennard-Jones
// energy and force evaluations.  Routines here are inline
// to allow the compiler to put them into the loop in forcefield.

#include <cmath>
#include <string>
#include <cstdlib>
#include "config.h"
#include "ffbase.h"

class CLJInt : public CFFBase {
 private:
  double       sigma, epsilon, A, B;
  double       locut, hicut, locut2, hicut2;
  std::string  id;

 public:
  /* Constructors */
  CLJInt() 
    : sigma(0.0), epsilon(0.0), A(0.0), B(0.0), locut(0.0), hicut(0.0),
      locut2(0.0), hicut2(0.0), id("") {};
  CLJInt(double sigma, double epsilon, std::string label) {
    Init(sigma, epsilon, label);
  }
  
  /* Destructor */
  ~CLJInt() { /* cout << " CLJInt destructor.." << endl; */ }

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
  /*           0 -- under range, returns big number */
  /*           1 -- ok, within range */
  /*           2 -- above range, returns zero */
  /*------------------------------------------------------------------------*/
  int Potential(const CVector3& sepvec, double& potential) {
    double r2,r2i,r6i,r12i;

    potential = 0.0;
    r2 = sepvec.SqLength();
    if (r2 > hicut2) {return(2);}
    if (r2 >= locut2) {
      r2i = 1.0/r2;
      r6i = r2i*r2i*r2i;
      r12i= r6i*r6i;
      potential = A*r12i - B*r6i;
      return(1);
    }

    potential = 1.0e30;    
    return(0);
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

    double r2,r2i,r6i,r12i;

    force = CVector3(0.0,0.0,0.0);
    r2 = sepvec.SqLength();
    if (r2 > hicut2) {return(true);}
    if (r2 >= locut2) {
      r2i = 1.0/r2;
      r6i = r2i*r2i*r2i;
      r12i= r6i*r6i;
      potential = A*r12i - B*r6i;
      force = sepvec * (12.0*A*r12i - 6.0*B*r6i)*r2i;
      return(true);
    }
    
    return(false);
  }
};

