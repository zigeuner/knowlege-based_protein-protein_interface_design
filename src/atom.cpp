// Functions for the atom class.  These routines allow manual manipulation
// of the atom structure.  

#ifdef _USE_OLD_LIBRARIES
#include <math.h>
#else
#include <cmath>
#endif
#include <iostream>
#include <fstream>
#include "include/config.h"
#include "include/atom.h"

using std::ifstream;
using std::ofstream;
using std::ostream;
using std::string;
using std::cout;
using std::ios;
using std::endl;

/*------------------------------------------------------------------------*/
/* Initializes from a PDB ATOM record line */
/* Requires:  line -- PDB line from which to initialize */
/*            verbose -- flag indicating if output should be written to scrn */
/*------------------------------------------------------------------------*/
void CAtom::Init(string line, bool verbose) {
  double    x,y,z;
  int       len = line.length();

  /* pad string with spaces to 80 characters long */
  if (len < 80) {
    for(int i=0;i<(80-len);i++) {line += " ";}    
  }

  repul = false;
  recordname = line.substr(0,6);
  Trim(recordname);
  atomno = ToInt(line.substr(6,5));
  atomname = line.substr(12,4);
  Trim(atomname);
  resid = line.substr(17,3);
  Trim(resid);
  if (resid.length() == 0) {resid = 'z';}
  chainid = line.substr(21,1);
  Trim(chainid);
  /* HACK, need to check for non-numerical residue numbers here? */
  resno = line.substr(22,5);
  Trim(resno);
  occ = ToDbl(line.substr(54,6));
  bfac = ToDbl(line.substr(60,6));
  segid = line.substr(72,4);
  Trim(segid);

  /* use segment ID as chain ID if chain ID is absent */
  if (chainid == "") {chainid = segid;}
  if (chainid == "") {
    //    cout << "ERROR: could not determine chain ID for atom" << endl;
    cout << "WARNING: could not determine chain ID for line, using 'X' " << endl;
    cout << line << endl;
    chainid = "X";
  }

  /* get the element from 77th, 13th or 12th columns */
  string   name[10];
  name[0] = line.substr(76,2);
  name[1] = line.substr(13,1);
  name[2] = line.substr(13,2);
  name[3] = line.substr(12,1);
  name[4] = line.substr(12,2);

  for(int i=0;i<5;i++) {
    mass = GetMass(name[i]);
    if (mass > 0.0) {break;}
  }
  if (mass <= 0.0) {
    cout << "WARNING: unable to find atomic mass (element= '" << element
	 << "') during init, using unity mass!" << endl;
    cout << " line: " << line << endl;
    //	exit(-1);
    mass = 1.0;
  }

  x = ToDbl(line.substr(30,8));
  y = ToDbl(line.substr(38,8));
  z = ToDbl(line.substr(46,8));
  r = CVector3(x,y,z);
}

/*------------------------------------------------------------------------*/
/* Initializes from an XYZ file line created by XYZLine */
/* Requires:  line -- XYZ line from which to initialize */
/*            verbose -- flag indicating if output should be written to scrn */
/*------------------------------------------------------------------------*/
void CAtom::InitfromXYZ(string line, bool verbose) {
  int       len, nfields;
  double    x,y,z;
  string    copy,field[1000];

  /* eliminate any leading or trailing spaces */
  copy = line;
  Trim(copy);
  len = copy.length();

  /* split the line */
  nfields = Split(copy,field,1000);
  if (nfields != 9) {
    cout << "ERROR: (atom.InitfromXYZ) expected 7 fields on input line:" << endl;
    cout << copy << endl;
    exit(-1);
  }

  recordname = "ATOM";

  element = field[0];
  mass = GetMass(element);
  if (mass <= 0.0) {
    element.erase(1,element.length()-1);
    mass = GetMass(element);
    if (mass <= 0.0) {
      cout << "ERROR: unable to find atomic mass during init" << endl;
      cout << " element= '" << element << "'" << endl;
      cout << line << endl;
      exit(-1);
    }
  }

  x = ToDbl(field[1]);
  y = ToDbl(field[2]);
  z = ToDbl(field[3]);
  r = CVector3(x,y,z);

  atomname = field[4];
  atomno = ToInt(field[5]);

  resid = field[6];
  if (resid.length() == 0) {resid = 'z';}

  /* HACK, need to check for non-numerical residue numbers here? */
  resno = field[7];
  chainid = field[8];

  /* fill-in stuff not expected to be present */
  occ = 0.0;
  bfac = 0.0;
  segid = chainid;
}

/*------------------------------------------------------------------------*/
/* Overloads the << operator to dump atom coordinates to screen */
/*------------------------------------------------------------------------*/
ostream& operator<<(ostream& os, const CAtom& atom) {
  os << atom.r;
  return os;
}

/*------------------------------------------------------------------------*/
/* Returns a string containing a PDB entry for the atom */
/*------------------------------------------------------------------------*/
string CAtom::PDBLine() {
  char       cstring[80],name[5];
  string     line;

  if (atomname.length() == 4) {
    sprintf(name,"%s",atomname.c_str());
  } else if (atomname.length() < 4) {
    sprintf(name," %-3s",atomname.c_str());
  } else {
    cout << "UNEXPECTED ERROR: length of atomname > 4" << endl;
    exit(-1);
  }

  sprintf(cstring,
	  "%-6s%5d %4s %3s %1s%4s    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s",
	  recordname.c_str(),atomno,name,resid.c_str(),
	  chainid.c_str(),resno.c_str(),r.X(),r.Y(),r.Z(),
	  occ,bfac,segid.c_str(),element.c_str());
  line = cstring;
  return(line);
}

/*------------------------------------------------------------------------*/
/* Returns a string containing an XYZ entry for the atom */
/* Requires:  calphaflag -- if true, labels C_alpha's and C_beta's as CA,CB */
/*------------------------------------------------------------------------*/
string CAtom::XYZLine(const bool calphaflag) {
  char       cstring[120];
  string     line,name;

  name = element;
  if (calphaflag) {
    if ((atomname == "CA")||(atomname == "CB")) {
      name = atomname;
    }
  }

  sprintf(cstring,"%-2s  %10.4f  %10.4f  %10.4f  %s  %d",
	  name.c_str(),r.X(),r.Y(),r.Z(),atomname.c_str(),atomno);

  line = cstring;
  return(line);
}

/*--------------------------------------------------------------------*/
/* A display routine to dump atom information to screen */
/* Requires:  indent -- number of spaces to indent */
/*--------------------------------------------------------------------*/
void CAtom::Display(int indent) {
  char  str[100];
  char  spacing[100];
  string  xtra = "";

  for(int i=0;i<indent;i++) {spacing[i] = ' ';}
  spacing[indent] = '\0';

  sprintf(str,"%s%4s %4s %.2f %.2f %.2f",spacing,atomname.c_str(),
	  resid.c_str(), r[0], r[1], r[2]);
  if (repul) {xtra = "  REPULSIVE";}
  cout << str << xtra << endl;
  //  cout << str << "  mass= " << mass << endl;
}

