#ifndef __RESIDUE_H
#define __RESIDUE_H
// This residue class contains all information about residue coordinates
// and associated derivatives.  It also contains utility functions for 
// reading, writing and manipulating these variables.

#include <fstream>
#include <string>
#include <cmath>
//#include "vector3.h"
#include "atom.h"

class CResidue {
 private:
  /* An array of atoms, each containing atomic specifications */
  int           natoms;    /* number of atoms in residue */
  int           resindex;  /* numeric type of residue */
  double        mass;    
  std::string   restype;   /* 3-letter residue type */
  std::string   resno;    /* residue number */
  CVector3      com;     /* Center of Mass */
  CAtom*        atom;

 public:  
  /* constructors */
  CResidue() 
    : natoms(0), resindex(0), mass(0.0), restype(), resno(), 
      com(0.0), atom(0) {};
  CResidue(const CResidue& old) {
    /* set number of atoms and allocate arrays */
    natoms = old.natoms;
    InitParameters();    

    /* Copy in the old values */
    Copy(old);
  };
    
  /* destructor */
  ~CResidue() {
    delete[] atom;
    //    atom = 0;
  }

  /* initialize default parameters all at once by calling this */
  void InitParameters() {
    /* Initialize space for the atoms */
    atom = new CAtom[natoms];
  }

  /* routines for returning state information*/  
  int NAtoms() {return natoms;}
  int ResIndex() {return resindex;}
  double Mass() {return mass;}
  std::string ResType() {return restype;}
  std::string ResNo() {return resno;}
  std::string ResSymbol() {return (AAType2Symbol(restype));}
  CVector3 COM() {return (com);}
  CVector3* GetAtomR(int num) { return (atom[num].GetR()); }
  CAtom* GetAtom(int num) { return &(atom[num]); }
  void SetResType(const std::string newtype) {restype = newtype;}
  void SetResNo(const std::string newresno) {resno = newresno;}
  bool Repulsive(const int num) {return(atom[num].Repulsive());}

  /* Functions declared in .cpp file */
  void Init(char*, int*, bool);
  void Init(CAtom*, int, int, bool);
  bool AtomRWType(const std::string&, CVector3&);
  bool GetAtomRWType(const std::string&, CVector3*&);
  int EndCarbonIndex();
  int CopyAtoms(const CResidue&, std::string*, int);
  int CountAtoms(std::string*, int);
  int NSCAtoms();
  bool CopyCentroid(CResidue&);
  CVector3 CalcCOM();
  void Display(int);
  void WritePDB(std::ofstream&);
  void WritePDB(char*);
  void WriteXYZ(std::ofstream&, std::string);

  friend std::ostream& operator<<(std::ostream& os, const CResidue& residue);

  /* Operator overloads */
  CResidue& operator=(const CResidue& old) {
    if (this != &old) {  /* avoid old = old problems */
      /* set N and allocate arrays */
      natoms = old.natoms;
      InitParameters();    

      /* Copy in the old values */
      Copy(old);
    }
    return *this;
  }

  /* Simply adds atoms together */
  /* NOTE1: takes definitions from res1! */
  /* NOTE2: COM not claculated here */
  friend CResidue operator+(CResidue& res1, CResidue& res2) {
    int       natoms = 0;
    CResidue  res;

    //    std::cout << "calling residue +" << std::endl;

    /* check residue types to make sure they're the same */
    if (res1.restype != res2.restype) {
      std::cout << "WARNING: (CResidue+) residue types are not the same: " 
		<< "res1= " << res1.restype 
		<< "  res2= " << res2.restype << std::endl;
    }

    /* size the internals of the new residue */
    res.natoms = res1.natoms + res2.natoms;
    res.resindex = res1.resindex;
    res.mass = res1.mass + res2.mass;
    res.restype = res1.restype;
    res.resno = res1.resno;
    res.com = 0.0;
    res.InitParameters();

    /* copy the atoms */
    for(int i=0; i<res1.natoms; i++) {    
      res.atom[natoms] = res1.atom[i];
      natoms += 1;
    }
    for(int i=0; i<res2.natoms; i++) {    
      res.atom[natoms] = res2.atom[i];
      natoms += 1;
    }

    return res;
  }

  /* Inline functions, here for speed */

  /*------------------------------------------------------------------------*/
  /* Copies one structure into another, without initializing the first */
  /* IMPORTANT: allocation sizes are assumed to match! natoms set elsewhere! */
  /* Requires:  old -- old structure to copy from */
  /*------------------------------------------------------------------------*/
  void Copy(const CResidue& old) {
    resindex = old.resindex;
    mass = old.mass;
    restype = old.restype;
    resno = old.resno;
    for(int i=0; i<natoms; i++) {    
      atom[i] = old.atom[i];
    }
    com = old.com;
  }

  /*------------------------------------------------------------------------*/
  /* Set all atoms to be evaluated as repulsive centers */
  /*------------------------------------------------------------------------*/
  void SetAllRepulsive() {
    for(int i=0; i<natoms; i++) {    
      atom[i].SetRepulsive(true);
    }
  }

  /*------------------------------------------------------------------------*/
  /* Set atoms of a give type to be evaluated as repulsive centers */
  /* Requires:  atomtypes -- array of atom types */
  /*            ntypes -- number of atom types specified */
  /*------------------------------------------------------------------------*/
  void SetAtomTypeRepulsive(std::string* atomtypes, int ntypes) {
    for(int i=0; i<natoms; i++) {    
      for(int j=0; j<ntypes; j++) {    
	if (atom[i].AtomName() == atomtypes[j]) {
	  atom[i].SetRepulsive(true);
	  continue;
	}
      }
    }
  }

};

#endif
