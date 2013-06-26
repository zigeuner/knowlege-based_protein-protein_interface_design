// This file provides the core routines for performing sub-system rotations

#include <string>
#include "config.h"
#include "movebase.h"
#include "random.h"

class CRotate : public CMoveBase {
 private:
  int           nmobile;  /* number of mobile chains */
  std::string*  mobile;   /* a list of mobile chains ID's */
  double        stepsize;
  double        phi, theta, psi;
  CMatrix3      rotmtx;   /* current rotation matrix */

 public:
  /* Constructors */
  CRotate() 
    : nmobile(0), mobile(0), stepsize(0.05), 
      phi(0.0), theta(0.0), psi(0.0), rotmtx(0.0) {};
  CRotate(const std::string chains) {
    Init(chains);
  }
  
  /* Destructor */
  ~CRotate() { /* cout << " CRotate destructor.." << endl; */ }

  /* Functions for returning the state variables of the object*/
  double StepSize() {return (stepsize);}  
  std::string ID() {return("ROTATE");}

  /* Functions declared in .cpp file */
  void Init(const std::string);
  bool Move(CSystem&, const CForcefield&);
  std::string Disp();
  void Disp(std::string&, std::string&);
  void Display();

  /* Simple Inline functions */
  void ResizeStep(const double factor) {stepsize *= factor;}
  void EnforceMinStep(const double minstepsize) {
    double   newstepsize = minstepsize*DEG2RAD;  
    if (stepsize < newstepsize) {stepsize = newstepsize;}
  }
  void EnforceMaxStep(const double maxstepsize) {
    double   newstepsize = maxstepsize*9.0*DEG2RAD;  /* max=10 => 90degrees */
    if (stepsize > newstepsize) {stepsize = newstepsize;}
  }
  void Reset() {stepsize = 0.05;}
};

