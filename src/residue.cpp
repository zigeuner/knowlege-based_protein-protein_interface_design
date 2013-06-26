// Functions for the residue class.  These routines allow manual manipulation
// of the residue structure.  

#ifdef _USE_OLD_LIBRARIES
#include <math.h>
#else
#include <cmath>
#endif
#include <iostream>
#include <fstream>
#include "include/config.h"
#include "include/residue.h"

using std::ifstream;
using std::ofstream;
using std::ostream;
using std::string;
using std::cout;
using std::ios;
using std::endl;

/*------------------------------------------------------------------------*/
/* Initializes, takes a set of atoms and a pair of numbers */
/* specifying which atoms belong to the residue */
/* Requires:  atomset -- an array of atoms */
/*            start -- first index to copy */
/*            end -- last index to copy */
/*            verbose -- flag indicating output to be written to screen */
/*------------------------------------------------------------------------*/
void CResidue::Init(CAtom* atomset, int start, int end, bool verbose) {
  int    i, atomno = 0;
  string line;

  natoms = end - start + 1;
  InitParameters();
  restype = atomset[start].ResID();
  resno = atomset[start].ResNo();
  resindex = AAType2IntIndex(restype);
  if (resindex < 0) {
    cout << "WARNING: residue with type '" << restype << 
      "' cannot be assigned a unambiguous residue index" << endl;
    cout << " it will be assigned index for residue type 'HET'" << endl;
    resindex = AAType2IntIndex("HET");
  }

  for(i=start; i<=end; i++) {
    atom[atomno] = atomset[i];
    atomno += 1;
  }

  if (verbose) {
    cout << "residue (" << restype << ") has " << natoms << " atoms" << endl;
    for(i=0; i<natoms; i++) {
      line = atom[i].PDBLine();
      cout << line << endl;
    }
  }
}

/*--------------------------------------------------------------------*/
/* Return the coordinates of the atom with a specific atom type */
/* Return false if it couldn't find an atom of that type */
/* Requires:  atomtype -- atom type (eg CA) */
/*--------------------------------------------------------------------*/
bool CResidue::AtomRWType(const string& atomtype, CVector3& coord) {

  if (atomtype == "MAXC") {
    int index = EndCarbonIndex();
    if (index < 0) {return(false);}
    coord = atom[index].R();
    return(true);
  } else {
    for(int i=0; i<natoms; i++) {
      //    cout << "'" << atom[i].AtomName() << "' =? " << atomtype << endl;
      //    atom[i].PDBLine();
      if (atom[i].AtomName() == atomtype) {
	coord = atom[i].R();
	return(true);
      }
    }
  } 
  return (false);
}

/*--------------------------------------------------------------------*/
/* Return the pointer to coordinates of the atom with a specific atom type */
/* Return false if it couldn't find an atom of that type */
/* Requires:  atomtype -- atom type (eg CA) */
/*--------------------------------------------------------------------*/
bool CResidue::GetAtomRWType(const string& atomtype, CVector3*& coord) {

  if (atomtype == "MAXC") {
    int index = EndCarbonIndex();
    if (index < 0) {return(false);}
    coord = atom[index].GetR();
    return(true);
  } else {
    for(int i=0; i<natoms; i++) {
      //    cout << "'" << atom[i].AtomName() << "' =? " << atomtype << endl;
      //    atom[i].PDBLine();
      if (atom[i].AtomName() == atomtype) {
	coord = atom[i].GetR();
	return(true);
      }
    }
  } 
  return (false);
}

/*--------------------------------------------------------------------*/
/* Return the index of the carbon atom at the "end" of the sidechain */
/* Return -1 if it couldn't find an atom of that type */
/*--------------------------------------------------------------------*/
int CResidue::EndCarbonIndex() {
  const int ntypes = 14;
  const std::string atype[ntypes] = 
  {"CB",
   "CG",
   "CG1",
   "CG2",
   "CD",
   "CD1",
   "CD2",
   "CE",
   "CE1",
   "CE2",
   "CE3",
   "CZ",
   "CZ2",
   "CZ3"};

  for(int a=ntypes-1; a>=0; a--) {
    for(int i=natoms-1; i>=0; i--) {
      if (atom[i].AtomName() == atype[a]) {
	//	cout << "found type '" << atype[a] << "' in " << restype << endl;
        return(i);
      }
    } 
  } 
  return(-1);

}

/*-------------------------------------------------------------------------*/
/* Recalculate the Center of Mass (COM) of the individual chains and the */
/* system as a whole.  Calculation of the masses is included here. */
/*------------------------------------------------------------------------*/
CVector3 CResidue::CalcCOM() {
  CVector3   vec(0.0);

  mass = 0.0;
  for(int i=0; i<natoms; i++) {    
    mass += atom[i].Mass();
    vec += atom[i].Mass()*atom[i].R();
  }
  if (mass <= 0.0) {
    cout << "ERROR: (residue.CalcCOM) residue mass is zero, something failed" 
	 << endl;
    Display(1);
    exit(-1);
  }
  com = vec/mass;

  return(com);
}

/*------------------------------------------------------------------------*/
/* Copies selected atom types of one residue structure into another */
/* Requires:  old -- old structure to copy from */
/*            atomtypes -- array of atom types to copy */
/*            ntypes -- number of atom types to copy */
/*------------------------------------------------------------------------*/
int CResidue::CopyAtoms(const CResidue& old, string* atomtypes, int ntypes) {
  int  ncopied = 0;

  restype = old.restype;
  resindex = old.resindex;
  resno = old.resno;

  /* count number of atoms of the proper type(s) */
  natoms = 0;
  for(int i=0; i<old.natoms; i++) {    
    for(int j=0; j<ntypes; j++) {
      if (old.atom[i].AtomName() == atomtypes[j]) {
	natoms += 1;
	break;
      }
    }
  }

  /* allocate arrays */
  InitParameters();    

  for(int i=0; i<old.natoms; i++) {    
    for(int j=0; j<ntypes; j++) {
      if (old.atom[i].AtomName() == atomtypes[j]) {
	atom[ncopied] = old.atom[i]; 
	ncopied += 1;
	break;
      }
    }
  }
  natoms = ncopied;

  return(ncopied);
}

/*------------------------------------------------------------------------*/
/* Copies only the sidechain of a residues as a single centroid atom */
/* Returns false if it could not find any sidechain atoms */
/* Requires:  old -- old structure to copy from */
/*------------------------------------------------------------------------*/
bool CResidue::CopyCentroid(CResidue& old) {
  double    mass;
  string    atomtype;
  CVector3  centroid,vec;

  if (old.natoms == 0) {return(false);}
  restype = old.restype;
  resindex = old.resindex;
  resno = old.resno;

  /* calculate the mass and centroid for the sidechain atoms */
  mass = 0.0;
  for(int i=0; i<old.natoms; i++) {    
    atomtype = old.atom[i].AtomName();
    if ((atomtype == "N")||(atomtype == "CA")||
	(atomtype == "C")||(atomtype == "O")) {continue;}
    mass += old.atom[i].Mass();
    vec += old.atom[i].Mass()*old.atom[i].R();
  }
  if (mass <= 0.0) {return(false);}
  centroid = vec/mass;

  /* allocate arrays */
  natoms = 1;
  InitParameters();    

  /* create the centroid by copying first atom, then changing it */
  atom[0] = old.atom[0]; 
  atom[0].SetMass(mass);
  atom[0].SetCoords(centroid);
  atom[0].SetAtomName("CTR");

  return(true);
}

/*------------------------------------------------------------------------*/
/* Count the number of atoms of given atom types in residue */
/* Requires:  atomtypes -- array of atom types to copy */
/*            ntypes -- number of atom types to copy */
/*------------------------------------------------------------------------*/
int CResidue::CountAtoms(string* atomtypes, int ntypes) {
  int  count = 0;

  /* count number of atoms of the proper type(s) */
  for(int i=0; i<natoms; i++) {    
    for(int j=0; j<ntypes; j++) {
      if (atom[i].AtomName() == atomtypes[j]) {
	count += 1;
	break;
      }
    }
  }

  return(count);
}

/*------------------------------------------------------------------------*/
/* Count the number of atoms sidechain atom types in residue */
/*------------------------------------------------------------------------*/
int CResidue::NSCAtoms() {
  int    count = 0;
  string atomtype;

  /* count number of atoms of the proper type(s) */
  for(int i=0; i<natoms; i++) {    
    atomtype = atom[i].AtomName();
    if ((atomtype == "N")||(atomtype == "CA")||
	(atomtype == "C")||(atomtype == "O")) {continue;}
    count += 1;
  }

  return(count);
}

/*------------------------------------------------------------------------*/
/* Write residue to a pdb file */
/* Requires:  pdbfile -- open file for configuration */
/*------------------------------------------------------------------------*/
void CResidue::WritePDB(ofstream& pdbfile) {
  string   line;

  /* Write each atom to pdb formatted file */
  for(int i=0; i<natoms; i++) {  
    line = atom[i].PDBLine();
    pdbfile << line << endl;
  }
}

/*------------------------------------------------------------------------*/
/* Write residue to a pdb file */
/* Requires:  file -- filename for configuration */
/*------------------------------------------------------------------------*/
void CResidue::WritePDB(char* file) {
  ofstream pdbfile;

  /* Open the pdb file for output*/
  pdbfile.open(file, ios::out);
  if(!pdbfile) {
    cout << endl << " ERROR: cannot open file '" << file << 
      "'" << endl;
    exit(-1);
  }
  //  printf("Writing residue configuration to pdb file '%s'\n",file);

  /* Write residue to pdb formatted file */
  WritePDB(pdbfile);
  pdbfile.close();
}

/*------------------------------------------------------------------------*/
/* Write residue to a xyz file */
/* Requires:  xyzfile -- open file for configuration */
/*            xtra -- extra string for end of line */
/*------------------------------------------------------------------------*/
void CResidue::WriteXYZ(ofstream& xyzfile, const string xtra) {
  string   line, info;

  /* Write each atom to xyz formatted file */
  for(int i=0; i<natoms; i++) {  
    if (natoms < 3) {
      line = atom[i].XYZLine(true);
    } else {
      line = atom[i].XYZLine(false);
    }
    info = "   " + restype + "  " + resno;
    if (xtra.length() > 0) {
      xyzfile << line << info << "  " << xtra << endl;
    } else {
      xyzfile << line << info << endl;
    }
  }
}

/*------------------------------------------------------------------------*/
/* Overloads the << operator to dump residue coordinates to screen */
/*------------------------------------------------------------------------*/
ostream& operator<<(ostream& os, const CResidue& residue) {
  register int i;

  for(i=0; i<residue.natoms; i++) {
    os << i+1 << "  " << residue.atom[i];
  }

  return os;
}

/*--------------------------------------------------------------------*/
/* A display routine to dump residue information to screen */
/* Requires:  indent -- number of spaces to indent */
/*--------------------------------------------------------------------*/
void CResidue::Display(int indent) {
  int   i;
  char  spacing[100];

  for(int i=0;i<indent;i++) {spacing[i] = ' ';}
  spacing[indent] = '\0';

  cout << spacing << "residue '" << resno << "' (" 
       << restype << "," << resindex << ") with " 
       << natoms << " atom(s):" << endl;
  for(i=0; i<natoms; i++) {
    atom[i].Display(indent + 1);
  }
}
