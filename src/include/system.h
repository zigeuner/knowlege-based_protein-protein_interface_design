#ifndef __SYSTEM_H
#define __SYSTEM_H

// This system class contains all information about system coordinates
// and associated derivatives.  It also contains utility functions for 
// read, writing and manipulating these variables.  Derived information
// for the system as a whole is also here

#include <iostream>
#include "vector3.h"
#include "chain.h"
#include "matrix.h"

class CSystem {
 private:
  /* system size specifications */
  int        nchains;
  double     mass;

  double     temperature;  /* system temperature */
  long       iseed;

  CChain*    chain;
  CVector3   com;          /* Center of Mass */

 public:
  /* Constructors */
  CSystem() : nchains(0), mass(0.0), temperature(300.0), 
	      iseed(371973), chain(0), com(0.0) {};
  CSystem(const CSystem& old) {
    /* set nchains and allocate arrays */
    SizeNChains(old.nchains);

    /* Copy in the old values */
    Copy(old);
  };
  CSystem(const char*, bool);

  /* Destructor */
  ~CSystem() {
    delete[] chain;
    //    chain = 0;
  }

  /* initialize default parameters all at once by calling this */
  void InitParameters() {
    nchains = 0;
  }

  /* Functions for returning the state variables of the object*/
  int NChains() {return (nchains);}
  double Mass() {return (mass);}
  double T() {return (temperature);}
  std::string ChainID(const int& num) { return (chain[num].ChainID()); }
  CVector3 COM() {return (com);}
  CVector3** GetChainR(const int& num) { return (chain[num].GetR()); }
  CVector3** GetChainV(const int& num) { return (chain[num].GetV()); }
  CVector3** GetChainF(const int& num) { return (chain[num].GetF()); }
  int NAtoms(const int& num) { return (chain[num].NAtoms()); }
  int NAtoms() { 
    int natoms = 0;
    for(int i=0; i<nchains; i++) {    
      natoms += chain[i].NAtoms();
    }
    return(natoms);
  }
  int MaxNAtoms() { 
    int max = 0, natoms;
    for(int i=0; i<nchains; i++) {    
      natoms = chain[i].NAtoms();
      if (natoms > max) max = natoms;
    }
    return(max);
  }
  int ChainNo(std::string id) { 
    for(int i=0; i<nchains; i++) {    
      if (id == chain[i].ChainID()) {return (i);}
    }
    return(-1);
  }
  long& GetSeed() { return iseed; } 
  void SetSeed(const long newseed) { iseed = newseed; } 
  double GetRan() { return ran2(iseed); }
  void SetT(const double value) {temperature = value;}

  /* Functions declared in the .cpp file */
  void ReadInit(const char*, bool);
  void ReadPDB(std::ifstream&, bool);
  void ReadXYZ(std::ifstream&, bool);
  void Pairs(std::string, std::string, std::string, std::string, 
	     double, std::string);
  void Pairs(std::string, std::string, std::string, std::string, 
	     double, std::string, std::ofstream&);
  int Pairs(std::string, std::string, std::string, 
	    std::string, double, std::string*, int&);
  void SomePairs(std::string, std::string, std::string, std::string, double, 
		 std::string, std::string*, int, std::ofstream&);
  int SomePairs(std::string, std::string, std::string, 
		std::string, double, std::string*, int, std::string*, int&);
  int SomePairAngles(std::string, std::string, double, std::string*,
		     int, std::string*, int&);
  void SomePairAngles(std::string, std::string, double cutoff, std::string label,
		      std::string* iresline, int nires, std::ofstream&);
  int GetCoords(std::string, std::string, std::string*, CVector3**, const int);
  int GetResidues(std::string, std::string*, CResidue**, const int);
  int GetResidueIndices(std::string, std::string*, int*, const int);
  bool GetResidue(std::string&, std::string&, CResidue*&);
  int GetRepulsiveAtoms(const std::string chainid, CAtom**& ptr) {
    chain[ChainNo(chainid)].GetRepulAtoms(ptr);
    return(chain[ChainNo(chainid)].NRepulAtoms());
  }
  void CopyReduce(const CSystem&, std::string, std::string*, int);
  void CopySubset(const CSystem&, CSubset&, std::string, std::string*, int);
  void CopySubset(const CSystem&, CSubset&);
  void CopyCentroids(const CSystem&, CSubset&);
  std::string ResSymbol(const std::string, const std::string);
  CVector3 CalcCOM();
  CVector3 COM(const std::string*, const int);
  void DistMtx(std::string, std::string, CSubset&, CMatrix&, std::string&);
  void DistSq(CSystem&, CSubset&, std::string, double*);
  void Dist(CSystem&, CSubset&, std::string, double*);
  double RMSD(CSystem&, CSubset&, std::string);
  bool CheckSubset(CSubset&);
  void Mutate(const std::string);
  void WritePDB(char*);
  void WritePDB(char*,const std::string);
  void WritePDB(std::ofstream&);
  void WritePDBSnapshot(const std::string);
  void WritePDBSnapshot();
  void WriteXYZ(char*, std::string);
  void WriteXYZ(std::ofstream&, std::string);
  void WriteProfitScripts(const std::string, const std::string);
  CSubset GenSubset();
  CSubset GenSubset(const std::string);
  void Display(int);

  /* Operator overloads */
  CSystem& operator=(const CSystem& old) {
    if (this != &old) {  /* avoid old = old problems */
      if (chain != NULL) {
	delete[] chain;
      }

      /* set nchains and allocate arrays */
      SizeNChains(old.nchains);

      /* Copy in the old values */
      Copy(old);
    }

    return *this;
  }

  friend CSystem operator+(CSystem& sys1, CSystem& sys2) {
    int      nchains = sys1.nchains, chainno;
    CSystem  sys;

    /* set basics from first system */
    sys.temperature = sys1.temperature;
    sys.iseed = sys1.iseed;

    /* look for additional chains in the second system */
    for(int c=0; c<sys2.nchains; c++) {
      if (sys1.ChainNo(sys2.chain[c].ChainID()) < 0) {nchains += 1;}
    }
    sys.SizeNChains(nchains);

    /* copy chains using ID list from first system */
    for(int c=0; c<sys1.nchains; c++) {
      if ((chainno = sys2.ChainNo(sys1.chain[c].ChainID())) >= 0) {
	sys.chain[c] = sys1.chain[c] + sys2.chain[chainno];
      } else {
	sys.chain[c] = sys1.chain[c];
      }
    }

    /* copy chains with unique ID from second system */
    nchains = sys1.nchains;
    for(int c=0; c<sys2.nchains; c++) {
      if (sys1.ChainNo(sys2.chain[c].ChainID()) < 0) {
	sys.chain[nchains] = sys2.chain[c];
	nchains += 1;
      }
    }

    return sys;
  }

  /* Inline Functions */

  /*------------------------------------------------------------------------*/
  /* Copies one structure into another, without initializing the first */
  /* IMPORTANT: allocation sizes are assumed to match! */
  /* Requires:  old -- old structure to copy from */
  /*------------------------------------------------------------------------*/
  void Copy(const CSystem& old) {
    nchains = old.nchains;
    iseed = old.iseed;
    mass = old.mass;
    temperature = old.temperature;
    com = old.com;

    for(int i=0; i<nchains; i++) {    
      chain[i] = old.chain[i];
    }
  }

  /*------------------------------------------------------------------------*/
  /* Size the chain array */
  /*------------------------------------------------------------------------*/
  void SizeNChains(int n) {
    nchains = n;
    chain = new CChain[n];
  } 

  /*------------------------------------------------------------------------*/
  /* Populate the coordinates, velocities and forces pointer arrays */
  /* NOTE: velocities and forces aren't currently used */
  /*------------------------------------------------------------------------*/
  void PopulateArrays() {
    for(int i=0; i<nchains; i++) {    
      chain[i].PopulateArrays();
    }
  }

  /*------------------------------------------------------------------------*/
  /* Set all atoms to be evaluated as repulsive centers */
  /*------------------------------------------------------------------------*/
  void SetAllRepulsive() {
    for(int i=0; i<nchains; i++) {    
      chain[i].SetAllRepulsive();
    }
  }

  /*------------------------------------------------------------------------*/
  /* Set atoms of a give type to be evaluated as repulsive centers */
  /* Requires:  atomtypes -- array of atom types */
  /*            ntypes -- number of atom types specified */
  /*------------------------------------------------------------------------*/
  void SetAtomTypeRepulsive(std::string* atomtypes, int ntypes) {
    for(int i=0; i<nchains; i++) {    
      chain[i].SetAtomTypeRepulsive(atomtypes,ntypes);
    }
  }

  /*------------------------------------------------------------------------*/
  /* Set atoms of a give type to be evaluated as repulsive centers */
  /* Requires:  chains -- comma-separated chain identifiers */
  /*            atomtypes -- array of atom types */
  /*            ntypes -- number of atom types specified */
  /*------------------------------------------------------------------------*/
  void SetSomeRepulsive(const std::string chains, 
			std::string* atomtypes, int ntypes) {
    std::string    chainlist[100];

    int nset = Split(chains,',',chainlist,100);
    for(int i=0; i<nchains; i++) {    
      for(int j=0; j<nset; j++) {    
	if (chain[i].ChainID() == chainlist[j]) {
	  chain[i].SetAtomTypeRepulsive(atomtypes,ntypes);
	  continue;
	}
      }
    }
  }

};

#endif
