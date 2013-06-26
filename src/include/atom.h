#ifndef __ATOM_H
#define __ATOM_H
// This atom class contains all information about atom coordinates
// and characterization information.  It also contains utility functions for 
// reading, writing and manipulating these variables.

#include <fstream>
#include <string>
#include <cmath>
#include "vector3.h"
#include "stringutils.h"

/* This character string function is here because I was having trouble */
/* with multiply defined functions when I put it in utils */
//int stripcomment(char*, char*);

class CAtom {
  /* NOTE: resid, chainid, resno and segid are NOT maintained in this */
  /* class, they are present here only for file input and output */
 private:
  int           atomno;
  bool          repul;   /* if true, can also be evaluated as a repulsive ctr */
  double        occ, bfac, mass;
  std::string   recordname, atomname, resid, chainid, resno, segid, element;
  CVector3      r;

 public:  
  /* constructors */
  CAtom() 
    : atomno(0), repul(0), occ(0.0), bfac(0.0), mass(0.0), recordname(), 
      atomname(), resid(), chainid(), resno(), segid(), element(), r() {};
  CAtom(const CAtom& old) {
    InitParameters();
    Copy(old);
  };
    
  /* destructor */
  ~CAtom() {
    //    r = 0;
  }

  /* initialize default parameters all at once by calling this */
  void InitParameters() {
  }

  /* routines for returning state information*/  
  int AtomNo() { return (atomno); }
  bool Repulsive() { return (repul); }
  double Occ() { return (occ); }
  double BFac() { return (bfac); }
  double Mass() { return (mass); }
  std::string ResNo() { return (resno); }
  std::string RecordName() { return (recordname); }
  std::string AtomName() { return (atomname); }
  std::string Element() { return (element); }
  std::string ChainID() { return (chainid); }
  std::string ResID() { return (resid); }
  std::string SegID() { return (segid); }
  CVector3 R() { return (r); }
  CVector3* GetR() { return (&r); }

  /* routines for setting state information*/  
  void SetMass(const double newmass) {mass = newmass;}
  void SetCoords(const CVector3 newcoords) {r = newcoords;}
  void SetAtomName(const std::string newname) {atomname = newname;}
  void SetRepulsive(const bool flag) {repul = flag;}

  /* Functions declared in .cpp file */
  void Init(std::string, bool);
  void InitfromXYZ(std::string, bool);
  std::string PDBLine();
  std::string XYZLine(const bool);
  void Display(int);

  friend std::ostream& operator<<(std::ostream& os, const CAtom& atom);

  /* Operator overloads */
  CAtom& operator=(const CAtom& old) {
    if (this != &old) {  /* avoid old = old problems */
      InitParameters();    
      /* Copy in the old values */
      Copy(old);
    }
    return *this;
  }

  /* Inline functions, here for speed */

  /*------------------------------------------------------------------------*/
  /* Copies one structure into another, without initializing the first */
  /* IMPORTANT: allocation sizes are assumed to match! */
  /* Requires:  old -- old structure to copy from */
  /*------------------------------------------------------------------------*/
  void Copy(const CAtom& old) {
    atomno = old.atomno;
    repul = old.repul;
    occ = old.occ;
    bfac = old.bfac;
    recordname = old.recordname;
    atomname = old.atomname;
    resid = old.resid;
    chainid = old.chainid;
    resno = old.resno;
    segid = old.segid;
    element = old.element;  
    mass = old.mass;
    r = old.r;
  }

};

#endif
