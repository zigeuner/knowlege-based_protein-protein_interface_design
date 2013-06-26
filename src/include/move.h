#ifndef __MOVE_H
#define __MOVE_H
// This move class is responsible for containing (directly or 
// indirectly) any necessary parameters and all move performing
// routines.

#include <string>
#include "config.h"
#include "rotate.h"
#include "translate.h"
 
class CMove {
 private:
  int           nmoves;
  CMoveBase**   movetype;
  int           nmobile;  /* number of mobile chains */
  std::string*  mobile;   /* a list of mobile chains ID's */

 public:
  /* constructors */
  CMove() 
    : nmoves(0), movetype(0), nmobile(0), mobile(0) {};
  CMove(std::string*, int, const std::string, bool);

  /* Destructor */
  ~CMove() {
    for(int i=0; i<nmoves; i++) {
      delete movetype[i];
    }
    delete[] movetype;
  }

  /* Functions */
  void Init(std::string*, int, const std::string, bool);
  void Move(CSystem&, CForcefield&, const double, const bool);
  void Move(CSystem&, CForcefield&, const bool);
  std::string Disp();
  void Disp(std::string&, std::string&);
  void Display(int);

  /* Simple inline functions */
  void ResizeSteps(const double factor) {
    for(int i=0; i<nmoves; i++) {    
      movetype[i]->ResizeStep(factor);
    }
  }
  void Reset() {
    for(int i=0; i<nmoves; i++) {    
      movetype[i]->Reset();
    }
  }

};

#endif
