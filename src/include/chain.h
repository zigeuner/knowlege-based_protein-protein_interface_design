#ifndef __CHAIN_H
#define __CHAIN_H
// This chain class contains all information about chain coordinates
// and associated derivatives.  It also contains utility functions for 
// reading, writing and manipulating these variables.

#include <fstream>
#include <string>
#include <cmath>
//#include "vector3.h"
#include "residue.h"
#include "subset.h"

class CChain {
 private:
  int         natoms;   /* number of atoms in chain */
  double      mass;
  std::string chainid;  /* ID of chain from pdb file*/
  CVector3    com;      /* Center of Mass */
  /* These variables store pointers to the 
     r - positions of chain atoms
     v - velocities of chain atoms
     f - forces on chain atoms */
  CVector3**  r;
  CVector3**  v;
  CVector3**  f;

  /* Information about repulsive-only atoms */
  int      nrepulatoms;
  CAtom**  repulatom;
  //  CAtom*  repulatom[10];

  /* An array of residues, each containing atomic specifications */
  int          nres;  /* number of residues in chain */
  int          maxres;  /* maximum residue number */
  CResidue*    res;

 public:  
  /* constructors */
  CChain() 
    : natoms(0), mass(0.0), chainid(" "), com(0.0), r(0), v(0), f(0), 
      nrepulatoms(0), repulatom(0), nres(0), maxres(0), res(0) {
    //    std::cout << "calling chain constructor" << std::endl;
  };
  CChain(const CChain& old) {
    /* set number of atoms and allocate arrays */
    natoms = old.natoms;
    nres = old.nres;
    maxres = old.maxres;
    InitParameters();    

    /* Copy in the old values */
    Copy(old);
  };
    
  /* destructor */
  ~CChain() {
    //    std::cout << "calling chain destructor" << std::endl;
    delete[] r;
    delete[] v;
    delete[] f;
    delete[] res;
    delete[] repulatom;
  }

  /* initialize default parameters all at once by calling this */
  /* NOTE: repulatoms is NOT allocated here! (see PopulateArrays) */
  void InitParameters() {
    /* Initialize space for the coordinates and their derivatives */
    r = new CVector3* [natoms];
    v = new CVector3* [natoms];
    f = new CVector3* [natoms];
    res = new CResidue[nres];
  }

  /* routines for returning state information*/  
  int NAtoms() {return natoms;}
  int NRes() {return nres;}
  int MaxRes() {return maxres;}
  double Mass() {return mass;}
  std::string ResType(const int& num) { return (res[num].ResType()); }
  std::string ResNo(const int& num) { return (res[num].ResNo()); }
  std::string ResSymbol(const int& num) { 
    return (AAType2Symbol(res[num].ResType())); 
  }
  std::string ChainID() { return (chainid); }
  void SetResType(const int& num, const std::string newtype) {
    res[num].SetResType(newtype); 
  }
  CVector3 COM() {return (com);}
  int NRepulAtoms() {return nrepulatoms;}
  void GetRepulAtoms(CAtom**& atoms) {atoms = repulatom;}
  CResidue* GetResidues() { return res; }
  CResidue GetResidue(const int& num) { return res[num]; }
  CResidue* GetResiduePtr(const int& num) { return &(res[num]); }
  CVector3 GetR(const int& num) { return (*r[num]); }
  CVector3 GetV(const int& num) { return (*v[num]); }
  CVector3 GetF(const int& num) { return (*f[num]); }
  CVector3** GetR() { return (r); }
  CVector3** GetV() { return (v); }
  CVector3** GetF() { return (f); }

  /* Functions declared in .cpp file */
  void Init(CAtom*,int,int,bool);
  bool AtomRWType(const int&, const std::string&, CVector3&);
  bool GetAtomRWType(const int&, const std::string&, CVector3*&);
  void setN(int);
  bool FindResWNum(std::string, int&);
  bool FindResWNum(std::string, CResidue*&);
  bool ResSymbolWNum(const std::string, std::string&);
  std::string Sequence();
  void CopyReduce(const CChain&, std::string*, int);
  void CopySubset(CChain&, CSubset&, std::string*, int);
  void CopySubset(CChain&, CSubset&);
  void CopyCentroids(CChain&, CSubset&);
  CVector3 CalcCOM();
  void WritePDB(std::ofstream&);
  void WriteXYZ(std::ofstream&);
  void WritePDB(char*);
  void Display(int);

  friend std::ostream& operator<<(std::ostream& os, const CChain& chain);

  /* Operator overloads */
  CChain& operator=(const CChain& old) {
    //    std::cout << "calling chain =" << std::endl;
    if (this != &old) {  /* avoid old = old problems */
      /* set N and allocate arrays */
      natoms = old.natoms;
      nrepulatoms = old.nrepulatoms;
      nres = old.nres;
      maxres = old.maxres;
      InitParameters();    

      /* Copy in the old values */
      Copy(old);
    }
    return *this;
  }

  friend CChain operator+(CChain& chain1, CChain& chain2) {
    int     nresidues = chain1.nres, resnum = -1;
    CChain  chain;
    //    std::cout << "calling chain +" << std::endl;
    /* look for additional residues in the second chain */
    for(int r=0; r<chain2.nres; r++) {
      if (!chain1.FindResWNum(chain2.res[r].ResNo(),resnum)) {
	nresidues += 1;
      }
    }

    /* size the internals of the new chain */
    chain.chainid = chain1.chainid;
    chain.nres = nresidues;
    chain.natoms = chain1.natoms + chain2.natoms;
    chain.InitParameters();

    /* copy chains using residue number list from first chain */
    /* if both have the same residue number, just add atoms together */
    for(int r=0; r<chain1.nres; r++) {
      if (chain2.FindResWNum(chain1.res[r].ResNo(),resnum)) {
	chain.res[r] = chain1.res[r] + chain2.res[resnum];
      } else {
	chain.res[r] = chain1.res[r];
      }
    }

    /* copy chains with unique residue numbers from second chain */
    nresidues = chain1.nres;
    for(int r=0; r<chain2.nres; r++) {
      if (!chain1.FindResWNum(chain2.res[r].ResNo(),resnum)) {
	chain.res[nresidues] = chain2.res[r];
	nresidues += 1;
      }
    }

    /* set the maximum residue */
    chain.SetMaxRes();

    /* Recalculate the masses */
    chain.CalcCOM();    

    /* Reset the pointers to the physical quantities */
    chain.PopulateArrays();    

    return chain;
  }

  /* Inline functions, here for speed */

  /*------------------------------------------------------------------------*/
  /* Copies one structure into another, without initializing the first */
  /* IMPORTANT: allocation sizes are assumed to match! */
  /* Requires:  old -- old structure to copy from */
  /*------------------------------------------------------------------------*/
  void Copy(const CChain& old) {
    int i;
    //    std::cout << "calling chain.Copy" << std::endl;
    chainid = old.chainid;
    mass = old.mass;
    com = old.com;

    /* Copy residues */
    for(i=0; i<nres; i++) {    
      res[i] = old.res[i];
    }

    /* Reset the pointers to the physical quantities */
    PopulateArrays();    
  }

  /*------------------------------------------------------------------------*/
  /* Set all atoms to be evaluated as repulsive centers */
  /*------------------------------------------------------------------------*/
  void SetAllRepulsive() {
    for(int i=0; i<nres; i++) {    
      res[i].SetAllRepulsive();
    }
  }

  /*------------------------------------------------------------------------*/
  /* Set atoms of a give type to be evaluated as repulsive centers */
  /* Requires:  atomtypes -- array of atom types */
  /*            ntypes -- number of atom types specified */
  /*------------------------------------------------------------------------*/
  void SetAtomTypeRepulsive(std::string* atomtypes, int ntypes) {
    for(int i=0; i<nres; i++) {    
      res[i].SetAtomTypeRepulsive(atomtypes,ntypes);
    }
  }

  /*------------------------------------------------------------------------*/
  /* Get maximum residue number in chain */
  /*------------------------------------------------------------------------*/
  void SetMaxRes() {
    int   num;
    maxres = -10000;
    for(int r=0; r<nres; r++) {
      num = FilterToInt(res[r].ResNo());
      //      std::cout << num << "   " << maxres << std::endl;
      if (num > maxres) {maxres = num;}
    }
  }

  /*------------------------------------------------------------------------*/
  /* Populate the coordinates, velocities and forces pointer arrays */
  /* NOTE1: repulatoms is allocated here! */
  /* NOTE2: velocities and forces aren't currently used */
  /*------------------------------------------------------------------------*/
  void PopulateArrays() {
    int   atom = 0;

    /* setup pointers to each atom in chain */
    for(int i=0; i<nres; i++) {    
      for(int j=0; j<res[i].NAtoms(); j++) {    
	r[atom] = res[i].GetAtomR(j);
	atom += 1;
	if (atom > natoms) {
	  std::cout << 
	    "ERROR: exceeded number of atoms in pointer array, unexpected"
		    << std::endl << "  array size= " << natoms << std::endl;
	  exit(-1);
	}
      }
    }

    /* find repulsive-only atoms */
    nrepulatoms = 0;
    for(int i=0; i<nres; i++) {    
      for(int j=0; j<res[i].NAtoms(); j++) {    
	if (res[i].Repulsive(j)) {nrepulatoms += 1;}
      }
    }

    /* allocate space */
    repulatom = new CAtom*[nrepulatoms];

    /* populate pointers with repulsive-only atoms */
    atom = 0;
    for(int i=0; i<nres; i++) {    
      for(int j=0; j<res[i].NAtoms(); j++) {    
	if (res[i].Repulsive(j)) {
	  repulatom[atom] = res[i].GetAtom(j);
	  atom += 1;
	  if (atom > nrepulatoms) {
	    std::cout << 
	      "ERROR: exceeded number of repulsive atoms in pointer array, "
		      << "unexpected" << std::endl 
		      << "  array size= " << nrepulatoms << std::endl;
	    exit(-1);
	  }
	}
      }
    }

  }

};

#endif
