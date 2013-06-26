#ifndef __FORCEFIELD_H
#define __FORCEFIELD_H
// This forcefield class is responsible for containing (directly or 
// indirectly) any necessary parameters and all forcefield evaluation
// routines.

#include <string>
#include "config.h"
#include "system.h"
#include "subset.h"
#include "gentypeint.h"
#include "kbtypeint.h"
#include "hardsphereint.h"
 
class CForcefield {
 private:
  double        sysnrg;  /* system energy from last evaluation */
  std::string   fftype;
  /* an array of pointers to specific pair type interaction matrices */
  int           npairtypes;
  CMtxIntBase** intset;
  bool          subsetseq;  /* enable use of alternative seq from subset */
  bool          evalchipot; /* enable chi-derived res-res potential calc */
  double        chidistcutoff;  /* CA-CA distance cutoff for chi pots */
  /* optional specifications for specific residue interactions */
  bool          intsubset, useunk;
  CSubset       subset;
  /* optional specifications for hardsphere interactions */
  CHardSphere   hspotential;
  /* residue type index in AAIntType array to use as reference state */
  int           refresindex;

 public:
  /* constructors */
  CForcefield() 
    : sysnrg(0.0), fftype(""), intset(), subsetseq(false), evalchipot(false),
      chidistcutoff(12.0), intsubset(false), useunk(false), subset(),
      hspotential() {};
  CForcefield(const std::string, const std::string, const std::string, bool);

  /* Destructor */
  ~CForcefield() {
    for(int i=0; i<npairtypes; i++) {
      delete intset[i];
    }
    delete[] intset;
  }

  double CurrentPotential() {return(sysnrg);}
  void SetCurrentPotential(const double nrg) {sysnrg = nrg;}
  void SetUnknownInt(const bool flag) {useunk = flag;}
  void UseSubsetSeq(const bool flag) {subsetseq = flag;}
  void EvalChiPot(const bool flag) {evalchipot = flag;}
  void SetChiDistCutoff(const double cutoff) {chidistcutoff = cutoff;}

  /* Functions declared in .cpp file */
  void Init(const std::string, const std::string, bool);
  void SetSubset(const CSubset&);
  void SetSubset(const std::string, bool);
  void SetChainChainInt(const std::string, const std::string, bool);
  int EvalAlphabet(CSystem&, std::string*&, const int);
  void OptimizeSeq(CSystem& sys);
  void CompareSeqs(CSystem& sys, const bool);
  void SetChainSeq(CSystem& sys, const std::string, const std::string);
  std::string AltSeq();
  void Display(int);
  double Potential(CSystem&);
  double Potential(CSystem&, CMatrix&, std::string&);
  double RepulsivePot(CSystem&);

};

#endif
