// Functions for the chain class.  These routines allow manual manipulation
// of the chain structure.  

#ifdef _USE_OLD_LIBRARIES
#include <math.h>
#else
#include <cmath>
#endif
#include <iostream>
#include <fstream>
#include "include/config.h"
#include "include/chain.h"

using std::ifstream;
using std::ofstream;
using std::ostream;
using std::string;
using std::cout;
using std::ios;
using std::endl;

/*------------------------------------------------------------------------*/
/* Initializes chain using an array of atoms and a specification of which */
/* ones to copy into the class */
/* Requires:  atom -- an array of atoms */
/*            start -- first index to copy */
/*            end -- last index to copy */
/*            verbose -- flag indicating output to be written to screen */
/*------------------------------------------------------------------------*/
void CChain::Init(CAtom* atom, int start, int end, bool verbose) {
  int     i,startres,endres,resno;
  string  current;

  natoms = end - start + 1;

  /* count the number of residues, looks for changes to identify residues */
  nres = 0;
  current = atom[start].ResNo();
  for(i=start; i<=end; i++) {
    if (atom[i].ResNo() != current) {
      //      cout << current << " -> " << atom[i].ResNo() << endl;
      nres += 1;
      current = atom[i].ResNo();
      //      cout << nres << ":" << atom[i-1].PDBLine() << endl;
    }
    if (i == end) {
      //      cout << current << " -> " << atom[i].ResNo() << endl;
      nres += 1;
      //      cout << nres << ":" << atom[i-1].PDBLine() << endl;
    }
  }

  /* size arrays and set defaults */
  InitParameters();
  chainid = atom[start].ChainID();

  /* skip setting the chain coordinates for now */
  /*  for(i=0; i<natoms; i++) {
    r[i] = atom[i].R();
  } 
  */ 

  if (verbose) {
    cout << "found " << nres << " residues in chain '" << chainid 
	 << "'" << endl;
  }

  /* initialize the residues */
  resno = 0;
  current = atom[start].ResNo();
  startres = start;
  for(i=start; i<=end; i++) {
    if (atom[i].ResNo() != current) {
      endres = i - 1;
      if (verbose) {
	cout << "initializing residue " << resno + 1 << " with atoms " 
	     << startres + 1 << " -> " << endres + 1 << endl;
      }
      res[resno].Init(atom,startres,endres,verbose);
      current = atom[i].ResNo();
      startres = i;
      resno += 1;
    }
    if (i == end) {
      endres = i;
      if (verbose) {
	cout << "initializing residue " << resno + 1 << " with atoms " 
	     << startres + 1 << " -> " << endres + 1 << endl;
      }
      res[resno].Init(atom,startres,endres,verbose);
      current = atom[i].ResNo();
      startres = i;
      resno += 1;
    }
  }

  /* Get maximum residue number in chain */
  SetMaxRes();

  /* Check for duplicate residue numbers, if found eliminate */
  string  resnumber,newresnumber;
  for(i=nres-1; i>=0; i--) {
    resnumber = res[i].ResNo();
    for(int j=i-1; j>=0; j--) {
      if (res[j].ResNo() == resnumber) {
	newresnumber = Int2String(MaxRes() + 1);
	/* if duplicated, relabel as last residue number plus one */
	cout << "WARNING: found duplicate residue number '" << resnumber
	     << "' in residue " << i << " of " << nres << " in chain "
	     << chainid << endl;
	cout << "  relabeling as residue number '" << newresnumber 
	     << "'" << endl;
	res[i].SetResNo(newresnumber);
	SetMaxRes();
	break;
      }
    }
  }

  /* Populate the coordinates pointer arrays */
  PopulateArrays();
}

/*--------------------------------------------------------------------*/
/* Return the coordinates of an atom in a given residue number */
/* of a specific atom type */
/* Returns false if it couldn't find an atom of that type */
/* Requires:  resno -- residue number */
/*            atomtype -- atom type (eg CA) */
/*            coord -- output coordinates */
/*--------------------------------------------------------------------*/
bool CChain::AtomRWType(const int& resno, const string& atomtype, 
			CVector3& coord) {
  if (! res[resno].AtomRWType(atomtype,coord) ) {return(false);}
  return (true);
}

/*--------------------------------------------------------------------*/
/* Return the pointer to coordinates of an atom in a given residue */
/* number of a specific atom type */
/* Returns false if it couldn't find an atom of that type */
/* Requires:  resno -- residue number */
/*            atomtype -- atom type (eg CA) */
/*            coord -- output coordinates */
/*--------------------------------------------------------------------*/
bool CChain::GetAtomRWType(const int& resno, const string& atomtype, 
			   CVector3*& coord) {
  if (! res[resno].GetAtomRWType(atomtype,coord) ) {return(false);}
  return (true);
}

/*------------------------------------------------------------------------*/
/* Copies one system structure into another, reducing the complexity of */
/* specified chains to only Calpha and Cbeta atoms */
/* Requires:  old -- old structure to copy from */
/*            atomtypes -- array of atom types to copy */
/*            ntypes -- number of atom types to copy */
/*------------------------------------------------------------------------*/
void CChain::CopyReduce(const CChain& old, string* atomtypes, int ntypes) {
  int     i,nnew;

  chainid = old.chainid;
  mass = old.mass;
  com = old.com;

  /* set number of residues and allocate array */
  nres = old.nres;
  res  = new CResidue[nres];

  /* Copy only some atoms of each residue */
  natoms = 0;
  for(i=0; i<nres; i++) {    
    nnew += res[i].CopyAtoms(old.res[i],atomtypes,ntypes);
    if (nnew == 0) {
      cout << "WARNING: did not copy any atoms for residue '" << 
	old.res[i].NAtoms() << "' in chain '" << chainid << "'" << endl;
    }
    natoms += nnew;
  }

  /* Now size the number of atoms */
  r  = new CVector3* [natoms];
  v  = new CVector3* [natoms];
  f  = new CVector3* [natoms];
}

/*------------------------------------------------------------------------*/
/* Copies one chain structure into another, retaining only the residues */
/* specified in the subset and optionally reducing the complexity to */
/* only given atom types */
/* Requires:  old -- old structure to copy from */
/*            subset -- system subset definition */
/*            atomtypes -- array of atom types to copy */
/*            ntypes -- number of atom types to copy */
/*------------------------------------------------------------------------*/
void CChain::CopySubset(CChain& old, CSubset& subset, 
			string* atomtypes, int ntypes) {
  int   i, expectnres;

  chainid = old.chainid;
  mass = old.mass;
  com = old.com;
  nres = 0;

  /* set number of residues (accounting for changes) and allocate array */
  expectnres = subset.NResidues(chainid);
  if (ntypes > 0) {
    for(i=0; i<old.NRes(); i++) {    
      if (!subset.InSubset(chainid,old.res[i].ResNo())) {continue;}
      if (old.res[i].CountAtoms(atomtypes,ntypes) == 0) {expectnres -= 1;}
    }
  }
  res = new CResidue[expectnres];

  //  cout << expectnres << " residues to copy in " << chainid << endl;
  //  cout << " natomtypes= " << ntypes << endl;

  /* Copy only some residues and atoms */
  natoms = 0;
  for(i=0; i<old.NRes(); i++) {    
    if (!subset.InSubset(chainid,old.res[i].ResNo())) {continue;}
    if (ntypes > 0) {
      if (old.res[i].CountAtoms(atomtypes,ntypes) == 0) {continue;}
      natoms += res[nres].CopyAtoms(old.res[i],atomtypes,ntypes);
    } else {
      res[nres] = old.res[i];
      natoms += res[nres].NAtoms();
    }
    nres += 1;
  }

  /* Now size the number of atoms in the coordinate pointer arrays */
  r  = new CVector3* [natoms];
  v  = new CVector3* [natoms];
  f  = new CVector3* [natoms];

  if (nres != expectnres) {
    cout << "ERROR: (chain.CopySubset) unexpected number of copied residues: " 
	 << "copied= " << nres << "  expected= " << expectnres << endl;
    cout << " Attempted to copy all the residues defined in the subset" << endl;
    cout << " Some residues may not have any of the atomtypes?  ";
    for(i=0; i<ntypes; i++) {    
      cout << atomtypes[i];
      if (i < ntypes - 1) {cout << ",";}
    }
    cout << endl;
    Display(0);
    exit(-1);
  }
}

/*------------------------------------------------------------------------*/
/* Copies one chain structure into another, retaining only the residues */
/* specified in the subset */
/* Requires:  old -- old structure to copy from */
/*            subset -- system subset definition */
/*------------------------------------------------------------------------*/
void CChain::CopySubset(CChain& old, CSubset& subset) { 
  int i, expectnres;

  chainid = old.chainid;
  mass = old.mass;
  com = old.com;
  nres = 0;

  /* set number of residues (accounting for changes) and allocate array */
  expectnres = subset.NResidues(chainid);
  res  = new CResidue[expectnres];

  //  cout << expectnres << " residues to copy in " << chainid << endl;

  /* Copy only some residues and atoms */
  natoms = 0;
  for(i=0; i<old.NRes(); i++) {    
    if (!subset.InSubset(chainid,old.res[i].ResNo())) {continue;}
    res[nres] = old.res[i];
    natoms += res[nres].NAtoms();
    nres += 1;
  }

  /* Now size the number of atoms in the coordinate pointer arrays */
  r  = new CVector3* [natoms];
  v  = new CVector3* [natoms];
  f  = new CVector3* [natoms];

  if (nres != expectnres) {
    cout << "ERROR: (chain.CopySubset) unexpected number of copied residues: " 
	 << "copied= " << nres << "  expected= " << expectnres << endl;
    cout << " Attempted to copy all the residues defined in the subset" << endl;
    Display(0);
    exit(-1);
  }
}

/*------------------------------------------------------------------------*/
/* Creates a new system with sidechain centroids from the old system for */
/* those residues that are specified in the subset */
/* Requires:  old -- old structure to copy from */
/*            subset -- system subset definition */
/*------------------------------------------------------------------------*/
void CChain::CopyCentroids(CChain& old, CSubset& subset) { 
  int i, resnum = 0;

  chainid = old.chainid;
  mass = old.mass;
  com = old.com;

  /* get number of residues for which centroids could be made */
  nres = 0;
  for(i=0; i<old.NRes(); i++) {    
    if (!subset.InSubset(chainid,old.res[i].ResNo())) {continue;}
    if (old.res[i].NSCAtoms() > 0) {nres += 1;}
  }

  /* set number of residues (accounting for changes) and allocate array */
  res  = new CResidue[nres];

  //  cout << nres << " residues to copy in " << chainid << endl;

  /* Copy only some residues and atoms */
  natoms = 0;
  for(i=0; i<old.NRes(); i++) {    
    if (!subset.InSubset(chainid,old.res[i].ResNo())) {continue;}
    if (old.res[i].NSCAtoms() == 0) {continue;}
    if (!res[resnum].CopyCentroid(old.res[i])) {
      cout << "WARNING: unexpected, could not create centroid for residue " 
	   << old.res[i].ResNo() << " in chain '" << chainid << "'" << endl;
      old.res[i].Display(0);
      continue;
    }
    natoms += 1;
    resnum += 1;
  }
  if (resnum != nres) {
    cout << "ERROR: (chain.CopyCentroids) unexpected number of copied centers: " 
	 << "copied= " << resnum << "  expected= " << nres << endl;
    exit(-1);
  }

  /* Now size the number of atoms in the coordinate pointer arrays */
  r  = new CVector3* [natoms];
  v  = new CVector3* [natoms];
  f  = new CVector3* [natoms];
}

/*-------------------------------------------------------------------------*/
/* Recalculate the Center of Mass (COM) of the individual residues and the */
/* chain as a whole.  Calculation of the masses is included here. */
/*------------------------------------------------------------------------*/
CVector3 CChain::CalcCOM() {
  double     resmass[nres];
  CVector3   vec(0.0), rescom[nres];

  mass = 0.0;
  com = 0.0;
  if (natoms == 0) {return(com);}
  for(int i=0; i<nres; i++) {    
    //    cout << "calculating COM for residue " << res[i].ResNo() << endl;
    rescom[i] = res[i].CalcCOM();
    resmass[i] = res[i].Mass();
    vec += resmass[i]*rescom[i];
    mass += resmass[i];
  }
  if (mass <= 0.0) {
    cout << "WARNING: (chain.CalcCOM) chain mass is zero, " 
	 << "cannot calculate COM" << endl;
    Display(1);
    com = 0.0;
  } else {
    com = vec/mass;
  }

  return(com);
}

/*------------------------------------------------------------------------*/
/* Write chain to a pdb file */
/* Requires:  pdbfile -- open file for configuration */
/*------------------------------------------------------------------------*/
void CChain::WritePDB(ofstream& pdbfile) {
  /* Write each residue to pdb formatted file */
  for(int i=0; i<nres; i++) {  
    res[i].WritePDB(pdbfile);
  }
}

/*------------------------------------------------------------------------*/
/* Write chain to an xyz file */
/* Requires:  xyzfile -- open file for configuration */
/*------------------------------------------------------------------------*/
void CChain::WriteXYZ(ofstream& xyzfile) {
  /* Write each residue to xyz formatted file */
  for(int i=0; i<nres; i++) {  
    res[i].WriteXYZ(xyzfile,chainid);
  }
}

/*------------------------------------------------------------------------*/
/* Write chain to a pdb file */
/* Requires:  file -- filename for configuration */
/*------------------------------------------------------------------------*/
void CChain::WritePDB(char* file) {
  ofstream pdbfile;

  /* Open the pdb file for output*/
  pdbfile.open(file, ios::out);
  if(!pdbfile) {
    cout << endl << " Error: cannot open file '" << file << 
      "'" << endl;
    exit(-1);
  }
  //  printf("Writing chain configuration to pdb file '%s'\n",file);

  /* Write chain to pdb formatted file */
  WritePDB(pdbfile);
  pdbfile << "TER" << endl;
  pdbfile << "END" << endl;

  pdbfile.close();
}

/*--------------------------------------------------------------------*/
/* Returns the 1-letter amino acid sequence of the chain */
/*--------------------------------------------------------------------*/
string CChain::Sequence() {
  int     i;
  string  seq = "";

  for(i=0; i<nres; i++) {
    seq += AAType2Symbol(res[i].ResType());
  }
  return(seq);
}

/*--------------------------------------------------------------------*/
/* Finds the residue index in array given the residue number string */
/* Returns false if it couldn't find it */
/*--------------------------------------------------------------------*/
bool CChain::FindResWNum(string resno, int& resnum) {
  resnum = -1;
  for(int i=0; i<nres; i++) {
    //    cout << i << "  '" << res[i].ResNo() << "'" << endl;
    if (res[i].ResNo() == resno) {
      resnum = i;
      return(true);
    }
  }
  return(false);
}

/*--------------------------------------------------------------------*/
/* Returns a pointer to a residue given the residue number string */
/* Returns false if it couldn't find it */
/* Requires:  resno -- residue number (string) */
/*            ptr -- pointer to residue class */
/*--------------------------------------------------------------------*/
bool CChain::FindResWNum(string resno, CResidue*& ptr) {
  for(int i=0; i<nres; i++) {
    //    cout << i << "  '" << res[i].ResNo() << "'" << endl;
    if (res[i].ResNo() == resno) {
      ptr = &(res[i]);
      return(true);
    }
  }
  return(false);
}

/*--------------------------------------------------------------------*/
/* Returns the symbol of a residue given the residue "number" */
/* Returns false if it couldn't find it */
/*--------------------------------------------------------------------*/
bool CChain::ResSymbolWNum(const string resno, string& symbol) {
  int resnum;
  symbol = "";
  if (!FindResWNum(resno,resnum)) {
    return (false);
  } else {
    symbol = ResSymbol(resnum);
  }
  return(true);
}

/*------------------------------------------------------------------------*/
/* Overloads the << operator to dump chain coordinates to screen */
/*------------------------------------------------------------------------*/
ostream& operator<<(ostream& os, const CChain& chain) {
  register int i;

  for(i=0; i<chain.natoms; i++) {
    os << i+1 << "  " << chain.r[i];
  }

  return os;
}

/*--------------------------------------------------------------------*/
/* A display routine to dump chain information to screen */
/* Requires:  indent -- number of spaces to indent */
/*--------------------------------------------------------------------*/
void CChain::Display(int indent) {
  int   i;
  int   count = 0;
  char  spacing[100];

  for(int i=0;i<indent;i++) {spacing[i] = ' ';}
  spacing[indent] = '\0';

  cout << spacing << "Chain '" << chainid << "' information: natoms = " 
       << natoms << ", nresidues = " << nres << endl;
  cout << spacing << "maximum residue number = " << maxres << endl;
  for(i=0; i<nres; i++) {
    res[i].Display(indent + 1);
    count += res[i].NAtoms();
    //    cout << spacing << "thru atom " << count << endl;
  }

  cout << spacing << "Chain '" << chainid << "' coordinates:" << endl;
  for(i=0; i<natoms; i++) {
    cout << spacing << " " << i << "/" << natoms << "  " << *(r[i]) << endl;
  }

  cout << spacing << "chain '" << chainid 
       << "' number of repulsive atoms: " << nrepulatoms << endl;
  for(i=0; i<nrepulatoms; i++) {
    cout << spacing << i << "/" << nrepulatoms << ": ";
    (*(repulatom[i])).Display(0);
  }

}

