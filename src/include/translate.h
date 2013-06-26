// This file provides the core routines for performing sub-system translation

#include <cmath>
#include <string>
#include "config.h"
#include "movebase.h"
#include "random.h"

class CTranslate : public CMoveBase {
 private:
  int           nmobile;  /* number of mobile chains */
  std::string*  mobile;   /* a list of mobile chains ID's */
  double        stepsize;
  double        dx, dy, dz;
  CVector3      dispvec;

 public:
  /* Constructors */
  CTranslate() 
    : nmobile(0), mobile(0), stepsize(0.1),
      dx(0.0), dy(0.0), dz(0.0), dispvec(0.0) {};
  CTranslate(const std::string chains) {
    Init(chains);
  }
  
  /* Destructor */
  ~CTranslate() { /* cout << " CTranslate destructor.." << endl; */ }

  /* Functions for returning the state variables of the object*/
  double StepSize() {return (stepsize);}  
  std::string ID() {return("TRANSLATE");}

  /* Functions declared in .cpp file */
  void Init(const std::string);
  bool Move(CSystem&, const CForcefield&);
  std::string Disp();
  void Disp(std::string&, std::string&);
  void Display();

  /* Simple Inline functions */
  void ResizeStep(const double factor) {stepsize *= factor;}
  void EnforceMinStep(const double minstepsize) {
    if (stepsize < minstepsize) {stepsize = minstepsize;}
  }
  void EnforceMaxStep(const double maxstepsize) {
    if (stepsize > maxstepsize) {stepsize = maxstepsize;}
  }
  void Reset() {stepsize = 0.1;}
};

