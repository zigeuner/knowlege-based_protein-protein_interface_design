// This system class contains all information about system coordinates
// and associated derivatives.  It also contains utility functions for 
// reading, writing and manipulating these variables.

#ifdef _USE_OLD_LIBRARIES
#include <math.h>
#else
#include <cmath>
#endif
#include <iostream>
#include <fstream>
#include "include/config.h"
#include "include/system.h"

using std::ifstream;
using std::ofstream;
using std::ostream;
using std::string;
using std::cout;
using std::ios;
using std::endl;
using std::cerr;

/*------------------------------------------------------------------------*/
/* Initialize the class from a pdb filename */
/* Requires:  rfile -- filename from which to read configuration */
/*            verbose -- flag indicating output to be written to screen */
/*------------------------------------------------------------------------*/
CSystem::CSystem(const char* rfile, bool verbose) {
  InitParameters();
  ReadInit(rfile, verbose);
}

/*------------------------------------------------------------------------*/
/* Initialize the portion of the class that can be read from a */
/* configuration file.  The remainder must remain as default */
/* Requires:  rfile -- filename from which to read configuration */
/*            verbose -- flag indicating output to be written to screen */
/*------------------------------------------------------------------------*/
void CSystem::ReadInit(const char* rfile, bool verbose) {
  string   filename = rfile, ending;
  ifstream config_file;

  /* Open the system configuration file */
  config_file.open(rfile, ios::in);
  if(!config_file) {
    cout << endl << " Error: cannot open file '" << rfile << 
      "'" << endl;
    exit(-1);
  }
  if (verbose) {
    printf("Reading system configuration from '%s'\n",rfile);
    fflush(stdout);
  }

  /* Read from configuration file */
  ending = FileEnding(filename);
  if (ending == "xyz") {
    ReadXYZ(config_file, verbose);
  } else if (ending == "pdb") {
    ReadPDB(config_file, verbose);
  } else {
    cout << "WARNING: (system.Init) could not determine input file type " 
	 << "for filename '" << filename << "'" << endl;
    cout << " it will assumed that it is a pdb file" << endl;
    ReadPDB(config_file, verbose);
  }

  config_file.close();

  /* setup the coordinate pointers to individual atoms */
  PopulateArrays();

  /* calculate the center of mass */
  CalcCOM();
}

/*------------------------------------------------------------------------*/
/* Load a set of chain coordinates from configuration file.  */
/* Requires:  file -- file to read from */
/*            verbose -- flag indicating output to be written to screen */
/*------------------------------------------------------------------------*/
void CSystem::ReadPDB(ifstream& file, bool verbose) {
  int      c, i, indx, start, end, natoms = 0, lineno = 0;
  char     buffer[1000];
  string   line, srchstr, srchstr2, current;
  CAtom*   atom;

  /* Initialize a temporary array of atoms */
  atom = new CAtom[MAXNATOMS];

  /* Read the coordinates */
  srchstr = "ATOM";
  srchstr2 = "HETATM";
  while (!file.eof()) {
    /* read the line */
    lineno += 1;
    file.getline(buffer,1000);
    line = buffer;

    indx = line.find(srchstr,0);
    if (indx < 0) {indx = line.find(srchstr2,0);}
    if ((indx >= 0)&&(indx <= 2)) {
      //      cout << lineno << ":  " << line << endl;
      atom[natoms].Init(line,verbose);
      //      cout << natoms << ":    " << atom[natoms].PDBLine() << endl;
      natoms += 1;
      if (natoms >= MAXNATOMS) {
	cout << endl << " Error: exceeded MAXNATOMS during pdb read " 
	     << natoms << endl;
	exit(-1);
      }
    }
  }

  /* count the number of chains (different segid's) */
  nchains = 1;
  current = atom[0].ChainID();
  for(i=0; i<natoms; i++) {
    if (atom[i].ChainID() != current) {
      //      cout << current << " -> " << atom[i].ChainID() << endl;
      nchains += 1;
      current = atom[i].ChainID();
    }
  }

  if (verbose) {
    cout << "found " << nchains << " chains " << endl;
  }

  /* Size the arrays and initialize the chains */
  SizeNChains(nchains);

  /* Sort the atoms into chains etc */
  /* moves sequentially through the atom list and makes a new chain */
  /* every time the chainID changes, will fail if they're mixed! */
  current = atom[0].ChainID();
  c = 0;
  start = 0;
  for(i=0; i<natoms; i++) {
    if ((atom[i].ChainID() != current)||(i == natoms - 1)) {
      end = i - 1;
      if (i == natoms - 1) {end = i;}
      if (verbose) {
	cout << "initializing chain " << current 
	     << " (" << c+1 << ") " << " with atoms " << start + 1
	     << " -> " << end + 1 << endl;
      }      
      chain[c].Init(atom,start,end,verbose);
      if (verbose) {cout << "sequence: " << chain[c].Sequence() << endl;}
      c += 1;
      current = atom[i].ChainID();
      start = i;
    }
  }

  /* delete temporary array of atoms, necessary? */
  delete[] atom; 
}

/*------------------------------------------------------------------------*/
/* Load a set of chain coordinates from an xyz configuration file.  */
/* Requires:  file -- file to read from */
/*            verbose -- flag indicating output to be written to screen */
/*------------------------------------------------------------------------*/
void CSystem::ReadXYZ(ifstream& file, bool verbose) {
  int      c, i, start, end, natoms = 0, lineno = 0;
  char     buffer[1000];
  string   line, srchstr, srchstr2, current;
  CAtom*   atom;

  /* Initialize a temporary array of atoms */
  atom = new CAtom[MAXNATOMS];

  /* Read the coordinates */
  while (!file.eof()) {
    /* read the line */
    lineno += 1;
    file.getline(buffer,1000);
    line = buffer;

    /* skip the first two lines (1st is natoms, 2nd is comment) */
    if (lineno < 3) {continue;}

    //    cout << lineno << ":  '" << line << "'" << endl;
    if (line.length() == 0) {continue;}
    atom[natoms].InitfromXYZ(line,verbose);
    //    cout << natoms << ":    " << atom[natoms].PDBLine() << endl;
    natoms += 1;
    if (natoms >= MAXNATOMS) {
      cout << endl << " Error: exceeded MAXNATOMS during pdb read " 
	   << natoms << endl;
      exit(-1);
    }
  }

  /* NOTE: from here down it's the same as ReadPDB */

  /* count the number of chains (different segid's) */
  nchains = 1;
  current = atom[0].ChainID();
  for(i=0; i<natoms; i++) {
    if (atom[i].ChainID() != current) {
      //      cout << current << " -> " << atom[i].ChainID() << endl;
      nchains += 1;
      current = atom[i].ChainID();
    }
  }

  if (verbose) {
    cout << "found " << nchains << " chains " << endl;
  }

  /* Size the arrays and initialize the chains */
  SizeNChains(nchains);

  /* Sort the atoms into chains etc */
  /* moves sequentially through the atom list and makes a new chain */
  /* every time the chainID changes, will fail if they're mixed! */
  current = atom[0].ChainID();
  c = 0;
  start = 0;
  for(i=0; i<natoms; i++) {
    if ((atom[i].ChainID() != current)||(i == natoms - 1)) {
      end = i - 1;
      if (i == natoms - 1) {end = i;}
      if (verbose) {
	cout << "initializing chain " << current 
	     << " (" << c+1 << ") " << " with atoms " << start + 1
	     << " -> " << end + 1 << endl;
      }      
      chain[c].Init(atom,start,end,verbose);
      if (verbose) {cout << "sequence: " << chain[c].Sequence() << endl;}
      c += 1;
      current = atom[i].ChainID();
      start = i;
    }
  }

  /* delete temporary array of atoms, necessary? */
  delete[] atom; 
}

/*------------------------------------------------------------------------*/
/* Wrapper for Pairs to allow output to screen */
/* Requires:  see main call to Pairs() */
/*            label -- string label to insert at end of output */
/*------------------------------------------------------------------------*/
void CSystem::Pairs(string chain1, string atomtype1, string chain2, 
		    string atomtype2, double cutoff, string label) {
  int       npairs, arraysize = 5000;
  string*   pairs;

  pairs = new string[arraysize];

  npairs = Pairs(chain1,atomtype1,chain2,atomtype2,cutoff,pairs,arraysize);

  for(int i=0;i<npairs;i++) {
    if (label.length() > 0) {
      pairs[i] += "  ";
      pairs[i] += label;
    }
    cout << pairs[i] << endl;
  }
}

/*------------------------------------------------------------------------*/
/* Wrapper for Pairs to allow output to file */
/* Requires:  see main call to Pairs() */
/*            label -- string label to insert at end of output */
/*------------------------------------------------------------------------*/
void CSystem::Pairs(string chain1, string atomtype1, string chain2, 
		    string atomtype2, double cutoff, string label, 
		    ofstream& outfile) {
  int       npairs, arraysize = 5000;
  string*   pairs;

  pairs = new string[arraysize];

  npairs = Pairs(chain1,atomtype1,chain2,atomtype2,cutoff,pairs,arraysize);

  for(int i=0;i<npairs;i++) {
    if (label.length() > 0) {
      pairs[i] += "  ";
      pairs[i] += label;
    }
    outfile << pairs[i] << endl;
  }
}

/*------------------------------------------------------------------------*/
/* Find pairs of atoms matching chain, residue and atom type */
/* specifications separated by less than a given distance  */
/* Requires:  chain1 -- chain 1 ID (eg L) */
/*            atomtype1 -- atom type 1 (eg CA) */
/*            chain2 -- chain 2 ID */
/*            atomtype1 -- atom type 2 */
/*            cutoff -- cutoff distance */
/*            pairs -- array of strings containing output */
/*------------------------------------------------------------------------*/
int CSystem::Pairs(string chain1, string atomtype1, string chain2, 
		   string atomtype2, double cutoff, string* pairs, 
		   int& arraysize) {
  int       npairs = 0;
  int 	    chain1no = ChainNo(chain1);
  int 	    chain2no = ChainNo(chain2);
  double    dist;
  char      cstr[100];
  string    symbol1, symbol2, resno1, resno2;
  CVector3  coord1, coord2;

  if (chain1no < 0) {
    cout << "ERROR: could not match stored chain for ID '" << 
      chain1 << "'" << endl;
    return(0);
  }
  if (chain2no < 0) {
    cout << "ERROR: could not match stored chain for ID '" << 
      chain2 << "'" << endl;
    return(0);
  }

  pairs[npairs] = "# list of atom pairs from chain " + chain1 + " atom type '" 
    + atomtype1 + "' interactions with chain " + chain2 + 
    " atom type '" + atomtype2 + "'";
  //  cout << pairs[npairs] << endl;
  npairs += 1;
  
  for(int i=0;i<chain[chain1no].NRes();i++) {
    if (! chain[chain1no].AtomRWType(i,atomtype1,coord1)) {
      //      cout << "WARNING: chain " << chain1 << " residue " << i 
      //   << " is missing atom type '" << atomtype1 << "'" << endl;
      continue;
    }
    symbol1 = chain[chain1no].ResSymbol(i);
    if (symbol1.length() == 0) {symbol1 = "z";}
    resno1 = chain[chain1no].ResNo(i);

    for(int j=0;j<chain[chain2no].NRes();j++) {
      if ((chain1no == chain2no)&&(i == j)) {continue;}

      if (! chain[chain2no].AtomRWType(j,atomtype2,coord2)) {
	//	cout << "WARNING: chain " << chain2 << " residue " << j 
	//	     << " is missing atom type '" << atomtype2 << "'" << endl;
	continue;
      } 
      symbol2 = chain[chain2no].ResSymbol(j);
      if (symbol2.length() == 0) {symbol2 = "z";}
      resno2 = chain[chain2no].ResNo(j);

      /* find distance */
      dist = coord1.Dist(coord2);

      /* write feedback if distance is within cutoff */
      if (dist < cutoff) {
	sprintf(cstr,"%-4s %s  %-4s %s  %.2f",resno1.c_str(),symbol1.c_str(),
		resno2.c_str(),symbol2.c_str(),dist);
	if (! npairs%100) {
	  cout << npairs << "  " << arraysize << endl;
	}
	pairs[npairs] = cstr;
	npairs += 1;

	if (npairs > arraysize - 1) {
	  cout << "DEBUG at " << __FILE__ << ":" << __LINE__ << endl;
	  cout << "ERROR: resize array routine doesn't work, " << 
	    "cannot resize beyond " << arraysize << endl;
	  exit(-1);
	  Resize(pairs,arraysize,arraysize+100);
	  arraysize += 100;
	}
      }
      
    }
  }

  return(npairs);
}

/*------------------------------------------------------------------------*/
/* Wrapper for SomePairs to allow output to file */
/* Requires:  see main call to SomePairs() */
/*            label -- string label to insert at end of output */
/*------------------------------------------------------------------------*/
void CSystem::SomePairs(string chain1, string atomtype1, string chain2, 
			string atomtype2, double cutoff, string label, 
			string* iresline, int nires, ofstream& outfile) {
  int       npairs, arraysize = 5000;
  string*   pairs;

  pairs = new string[arraysize];

  npairs = SomePairs(chain1,atomtype1,chain2,atomtype2,cutoff,iresline,nires,
		     pairs,arraysize);

  for(int i=0;i<npairs;i++) {
    if (label.length() > 0) {
      pairs[i] += "  ";
      pairs[i] += label;
    }
    outfile << pairs[i] << endl;
  }
}

/*------------------------------------------------------------------------*/
/* Find pairs of atoms matching chain, residue and atom type */
/* specifications separated by less than a given distance  */
/* Same as 'Pairs', but this routine only takes residues from a list */
/* Requires:  chain1 -- chain 1 ID (eg L) */
/*            atomtype1 -- atom type 1 (eg CA) */
/*            chain2 -- chain 2 ID */
/*            atomtype1 -- atom type 2 */
/*            cutoff -- cutoff distance */
/*            iresline -- array of strings containing interacting residues */
/*            nires -- size of residues array */
/*            pairs -- array of strings containing output */
/*            arraysize -- number of output lines */
/*------------------------------------------------------------------------*/
int CSystem::SomePairs(string chain1, string atomtype1, string chain2, 
		       string atomtype2, double cutoff, string* iresline, 
		       int nires, string* pairs, int& arraysize) {
  int       npairs = 0, nfields, rnum1, rnum2;
  int 	    chain1no = ChainNo(chain1);
  int 	    chain2no = ChainNo(chain2);
  double    dist;
  char      cstr[100];
  string    symbol1, symbol2, resno1, resno2, resno, field[10];
  CVector3  coord1, coord2;

  if (chain1no < 0) {
    cout << "ERROR: could not match stored chain for ID '" << 
      chain1 << "'" << endl;
    return(0);
  }
  if (chain2no < 0) {
    cout << "ERROR: could not match stored chain for ID '" << 
      chain2 << "'" << endl;
    return(0);
  }

  pairs[npairs] = "# list of atom pairs from chain " + chain1 + " atom type '" 
    + atomtype1 + "' interactions with chain " + chain2 + 
    " atom type '" + atomtype2 + "'";
  //  cout << pairs[npairs] << endl;
  npairs += 1;
  
  for(int i=0;i<nires;i++) {
    /* compare with interface residue specifications, skip if no match */
    nfields = Split(iresline[i],' ',field,10);
    //    cout << field[1] << "=?" << chain1 << endl;
    if (field[1] != chain1) {continue;}
    resno = field[2];
    if (!chain[chain1no].FindResWNum(resno,rnum1)) {
      cout << "WARNING: could not find residue '" << resno << "' in chain '"
	   << chain1 << "'" << endl;
      continue;
    }

    if (! chain[chain1no].AtomRWType(rnum1,atomtype1,coord1)) {
      //      cout << "WARNING: chain " << chain1 << " residue " << rnum1
      //   << " is missing atom type '" << atomtype1 << "'" << endl;
      continue;
    }

    symbol1 = chain[chain1no].ResSymbol(rnum1);
    if ((symbol1.length() == 0)||(symbol1[0] == ' ')) {symbol1 = "z";}
    resno1 = chain[chain1no].ResNo(rnum1);

    for(int j=0;j<nires;j++) {

      /* compare with interface residue specifications, skip if no match */
      nfields = Split(iresline[j],' ',field,10);
      if (field[1] != chain2) {continue;}
      resno = field[2];
      if (!chain[chain2no].FindResWNum(resno,rnum2)) {
      	cout << "WARNING: could not find residue '" << resno << "' in chain '"
      	     << chain2 << "'" << endl;
      	continue;
      } 

      if (! chain[chain2no].AtomRWType(rnum2,atomtype2,coord2)) {
	//	cout << "WARNING: chain " << chain2 << " residue " << rnum2
	//	     << " is missing atom type '" << atomtype2 << "'" << endl;
	continue;
      } 
      symbol2 = chain[chain2no].ResSymbol(rnum2);
      if ((symbol2.length() == 0)||(symbol2[0] == ' ')) {symbol2 = "z";}
      resno2 = chain[chain2no].ResNo(rnum2);

      /* find distance */
      dist = coord1.Dist(coord2);

      //      cout << iresline[i] << "---"
      //	   << chain1 << rnum1 << symbol1 << ": " << coord1;
      //      cout << iresline[j] << "---"
      //	   << chain2 << rnum2 << symbol2 << ": " << coord2;
      //      cout << "  dist=" << dist << endl;

      /* write feedback if distance is within cutoff */
      if (dist < cutoff) {
	sprintf(cstr,"%-4s %s  %-4s %s  %.2f",resno1.c_str(),symbol1.c_str(),
		resno2.c_str(),symbol2.c_str(),dist);
	if (! npairs%100) {
	  cout << npairs << "  " << arraysize << endl;
	}
	pairs[npairs] = cstr;
	npairs += 1;

	if (npairs > arraysize - 1) {
	  cout << "DEBUG at " << __FILE__ << ":" << __LINE__ << endl;
	  cout << "ERROR: resize array routine doesn't work, " << 
	    "cannot resize beyond " << arraysize << endl;
	  exit(-1);
	  Resize(pairs,arraysize,arraysize+100);
	  arraysize += 100;
	}
      }
      
    }
  }

  return(npairs);
}

/*------------------------------------------------------------------------*/
/* Wrapper for SomePairAngles to allow output to file */
/* Requires:  see main call to SomePairAngles() */
/*            label -- string label to insert at end of output */
/*------------------------------------------------------------------------*/
void CSystem::SomePairAngles(string chain1, string chain2, double cutoff, 
			     string label, string* iresline, int nires, 
			     ofstream& outfile) {
  int       npairs, arraysize = 5000;
  string*   pairs;

  pairs = new string[arraysize];

  npairs = SomePairAngles(chain1,chain2,cutoff,iresline,nires,pairs,arraysize);

  for(int i=0;i<npairs;i++) {
    if (label.length() > 0) {
      pairs[i] += "  ";
      pairs[i] += label;
    }
    outfile << pairs[i] << endl;
  }
}

/*------------------------------------------------------------------------*/
/* Find pairs of residues matching chain ID's with CA's separated by less */
/* than a given distance.  Only residues from the passed list are used.  */
/* For each pair, a set of angles (and other) properties are calculated */
/* and returned in string format */
/* Requires:  chain1 -- chain 1 ID (eg L) */
/*            chain2 -- chain 2 ID */
/*            cutoff -- cutoff distance */
/*            iresline -- array of strings containing interacting residues */
/*            nires -- size of residues array */
/*            pairs -- array of strings containing output */
/*            arraysize -- number of output lines */
/*------------------------------------------------------------------------*/
int CSystem::SomePairAngles(string chain1, string chain2, double cutoff, 
			    string* iresline, int nires, string* pairs, 
			    int& arraysize) {
  int       npairs = 0, nfields, rnum1, rnum2;
  int 	    chain1no = ChainNo(chain1);
  int 	    chain2no = ChainNo(chain2);
  double    dist, cbcbdist, torsion, colinear_angle, angle1, angle2;
  double    dASA1, dASA2;
  char      cstr[100];
  string    symbol1, symbol2, resno1, resno2, resno, field[10];
  CVector3  calpha1, calpha2, cbeta1, cbeta2, a, b, c;

  if (chain1no < 0) {
    cout << "ERROR: could not match stored chain for ID '" << 
      chain1 << "'" << endl;
    return(0);
  }
  if (chain2no < 0) {
    cout << "ERROR: could not match stored chain for ID '" << 
      chain2 << "'" << endl;
    return(0);
  }

  pairs[npairs] = "# list of pairs from chain '" + chain1 +
    "' interactions with chain '" + chain2 + "' of structure ";
  //  cout << pairs[npairs] << endl;
  npairs += 1;
  
  for(int i=0;i<nires;i++) {

    /* compare with interface residue specifications, skip if no match */
    nfields = Split(iresline[i],' ',field,10);
    //    cout << field[1] << "=?" << chain1 << endl;
    if (field[1] != chain1) {continue;}
    resno = field[2];
    if (!chain[chain1no].FindResWNum(resno,rnum1)) {
      cout << "WARNING: could not find residue '" << resno << "' in chain '"
	   << chain1 << "'" << endl;
      continue;
    }
    dASA1 = 0.0;
    if (nfields > 4) {dASA1 = ToDbl(field[4]);}

    /* Extract Calpha and Cbeta coordinates for first residue */
    if (! chain[chain1no].AtomRWType(rnum1,"CA",calpha1)) {
      //      cout << "WARNING: chain " << chain1 << " residue " << rnum1
      //   << " is missing atom type '" << atomtype1 << "'" << endl;
      continue;
    }
    if (! chain[chain1no].AtomRWType(rnum1,"CB",cbeta1)) {
      //      cout << "WARNING: chain " << chain1 << " residue " << rnum1
      //   << " is missing atom type '" << atomtype1 << "'" << endl;
      continue;
    }

    symbol1 = chain[chain1no].ResSymbol(rnum1);
    if ((symbol1.length() == 0)||(symbol1[0] == ' ')) {symbol1 = "z";}
    resno1 = chain[chain1no].ResNo(rnum1);

    for(int j=0;j<nires;j++) {

      /* compare with interface residue specifications, skip if no match */
      nfields = Split(iresline[j],' ',field,10);
      if (field[1] != chain2) {continue;}
      resno = field[2];
      if (!chain[chain2no].FindResWNum(resno,rnum2)) {
      	cout << "WARNING: could not find residue '" << resno << "' in chain '"
      	     << chain2 << "'" << endl;
      	continue;
      } 
      dASA2 = 0.0;
      if (nfields > 4) {dASA2 = ToDbl(field[4]);}

      /* Extract Calpha and Cbeta coordinates for second residue */
      if (! chain[chain2no].AtomRWType(rnum2,"CA",calpha2)) {
        //      cout << "WARNING: chain " << chain1 << " residue " << rnum1
        //   << " is missing atom type '" << atomtype1 << "'" << endl;
        continue;
      }
      if (! chain[chain2no].AtomRWType(rnum2,"CB",cbeta2)) {
        //      cout << "WARNING: chain " << chain1 << " residue " << rnum1
        //   << " is missing atom type '" << atomtype1 << "'" << endl;
        continue;
      } 

      symbol2 = chain[chain2no].ResSymbol(rnum2);
      if ((symbol2.length() == 0)||(symbol2[0] == ' ')) {symbol2 = "z";}
      resno2 = chain[chain2no].ResNo(rnum2);

      /* find distance */
      dist = calpha1.Dist(calpha2);

      //      cout << iresline[i] << "---"
      //	   << chain1 << rnum1 << symbol1 << ": " << calpha1;
      //      cout << iresline[j] << "---"
      //	   << chain2 << rnum2 << symbol2 << ": " << calpha2;
      //      cout << "  dist=" << dist << endl;

      /* write feedback if distance is within cutoff */
      if (dist < cutoff) {

	/* find CB-CB distance */
	cbcbdist = cbeta1.Dist(cbeta2);

	/* calculate torsion angle */
	a = cbeta1 - calpha1;
	b = calpha2 - cbeta1;
	c = cbeta2 - calpha2;
	torsion = a.TorsionAngle(b,c)*RAD2DEG;

	/* calculate colinearity angle (translate so cbeta's overlap) */
	c = calpha2 - (cbeta2 - cbeta1);
	colinear_angle = calpha1.Angle(cbeta1,c)*RAD2DEG;

	/* calculate angle btwn alpha1-beta1 and beta2 */
	angle1 = calpha1.Angle(cbeta1,cbeta2)*RAD2DEG;

	/* calculate angle btwn alpha1-beta1 and beta2 */
	angle2 = calpha2.Angle(cbeta2,cbeta1)*RAD2DEG;

	sprintf(cstr,"%-4s %s  %-4s %s   %s  %s   %s  %s   %s  %s  %s  %s",
		resno1.c_str(),symbol1.c_str(),
		resno2.c_str(),symbol2.c_str(),
		(Num(dist)).c_str(),(Num(cbcbdist)).c_str(),
		(Num(dASA1)).c_str(),(Num(dASA2)).c_str(),
		(Num(torsion)).c_str(), (Num(colinear_angle)).c_str(),
		(Num(angle1)).c_str(), (Num(angle2)).c_str());
	if (! npairs%100) {
	  cout << npairs << "  " << arraysize << endl;
	}
	//	cout << cstr << endl;
	//	cout << "torsion: " << Num(torsion) << endl;
	//	cout << "angle: " << Num(colinear_angle) << endl;
	//	cout << "angle1: " << Num(angle1) << endl;
	//	cout << "angle2: " << Num(angle2) << endl;
	pairs[npairs] = cstr;
	npairs += 1;

	if (npairs > arraysize - 1) {
	  cout << "DEBUG at " << __FILE__ << ":" << __LINE__ << endl;
	  cout << "ERROR: resize array routine doesn't work, " << 
	    "cannot resize beyond " << arraysize << endl;
	  exit(-1);
	  Resize(pairs,arraysize,arraysize+100);
	  arraysize += 100;
	}
      }
      
    }
  }

  return(npairs);
}

/*------------------------------------------------------------------------*/
/* Constructs a distance matrix using a subset definition */
/* The passed matrix class MUST ALREADY BE SIZED to the number */
/* of interacting residues in the forcefield's subset definition */
/* Requires:  atomtype1 -- atom type 1 (eg CA) */
/*            atomtype2 -- atom type 2 */
/*            subset -- system subset definition */
/*            mtx -- matrix of res-res distances */
/*            commentline -- string of residue identifiers for output */
/*------------------------------------------------------------------------*/
void CSystem::DistMtx(string atomtype1, string atomtype2, CSubset& subset, 
		      CMatrix& mtx, string& commentline) {
  int       rnum1, rnum2, nres, tempno, chain1no, chain2no;
  string    symbol1, symbol2, resno1, resno2, resno, chain1, chain2;
  CVector3  coord1, coord2;

  if (subset.TotalResidues() != mtx.Size()) {
    cout << "ERROR: passed matrix class must have size corresponding to "
	 << "number of interface residues" << endl;
    exit(-1);
  }

  nres = subset.TotalResidues();

  for(int i=0;i<nres;i++) {

    /* Get chain and residue number */
    tempno = subset.Num2ID(i,chain1,resno1);
    chain1no = ChainNo(chain1);
    if (!chain[chain1no].FindResWNum(resno1,rnum1)) {
      cout << "WARNING: could not find residue '" << resno << "' in chain '"
	   << chain1 << "'" << endl;
      continue;
    }

    if (! chain[chain1no].AtomRWType(rnum1,atomtype1,coord1)) {
      //      cout << "WARNING: chain " << chain1 << " residue " << rnum1
      //   << " is missing atom type '" << atomtype1 << "'" << endl;
      continue;
    }

    for(int j=i+1;j<nres;j++) {

      /* Get chain and residue number */
      tempno = subset.Num2ID(j,chain2,resno2);
      chain2no = ChainNo(chain2);
      if (!chain[chain2no].FindResWNum(resno2,rnum2)) {
      	cout << "WARNING: could not find residue '" << resno << "' in chain '"
      	     << chain2 << "'" << endl;
      	continue;
      } 

      if (! chain[chain2no].AtomRWType(rnum2,atomtype2,coord2)) {
	//	cout << "WARNING: chain " << chain2 << " residue " << rnum2
	//	     << " is missing atom type '" << atomtype2 << "'" << endl;
	continue;
      } 

      /* find and record distance */
      mtx.Set(i,j,coord1.Dist(coord2));

    }
  }

  /* make the string of residue definitions as a comment line */
  string tempid,tempresno;
  commentline = "# ";
  for(int i=0; i<subset.TotalResidues(); i++) {
    tempno = subset.Num2ID(i,tempid,tempresno);
    commentline += tempid + ":" + ResSymbol(tempid,tempresno);
    commentline += tempresno +  " ";
  }  
}

/*------------------------------------------------------------------------*/
/* Fill an array of CVector3 pointers correponding to coordinates of  */
/* a single atom type of specified residue of a given chain */
/* This is used get the coordinates of interacting atoms */
/* Returns the number of coordinates that it successfully found */
/* Requires:  chainid -- chain 1 ID (eg L) */
/*            atomtype -- atom type 1 (eg CA) */
/*            resnumbers -- array of residue numbers (strings) */
/*            ptarray -- array of CVector3 pointers, pointing to atom coords */
/*            arraysize -- number of atoms */
/*------------------------------------------------------------------------*/
int CSystem::GetCoords(string chainid, string atomtype, string* resnumbers,
		       CVector3** ptarray, const int arraysize) {
  int     chainno = ChainNo(chainid);
  int     rnum, ncoords = 0;
  string  resno;

  if (chainno < 0) {
    cout << "ERROR: could not match stored chain for ID '" << 
      chainid << "'" << endl;
    return(0);
  }

  /* Search for the residue then the atom with the correct labels */
  for(int i=0; i<arraysize; i++) {
    //    cout << resnumbers[i] << " in " << chainid << endl;
    /* find the residue in the chain */
    if (!chain[chainno].FindResWNum(resnumbers[i],rnum)) {
      cout << "WARNING: (system.GetCoords) could not find residue number '" 
	   << resnumbers[i] << "' in chain '" 
	   << chainid << "', skipping this residue" << endl;
      continue;
    }

    /* find the atom type in the residue */
    if (chain[chainno].GetAtomRWType(rnum,atomtype,ptarray[ncoords])) {
      ncoords += 1;
    } else {
      //      cout << "WARNING: could not find atomtype '" << atomtype << 
      //	"' in residue '" << resnumbers[i] << "' of chain '" << 
      //	chainid << "'" << endl;
    }
  }

  return(ncoords);
}

/*------------------------------------------------------------------------*/
/* Return the pointer to a specific residue */
/* Returns false if it could not find the residue */
/* Requires:  chainid -- chain ID (eg L) */
/*            resno -- residue number (string) */
/*            ptr -- pointer to residue class */
/*------------------------------------------------------------------------*/
bool CSystem::GetResidue(string& chainid, string& resno, CResidue*& ptr) {
  int     chainno = ChainNo(chainid);

  if (chainno < 0) {
    cout << "ERROR: could not match stored chain for ID '" << 
      chainid << "'" << endl;
    return(false);
  }

  /* Search chain for the residue */
  if (chain[chainno].FindResWNum(resno,ptr)) {
    return(true);
  } else {
    return(false);
  }
}

/*------------------------------------------------------------------------*/
/* Fill an array of CResidue pointers correponding to a list of residue */
/* "numbers" supplied. */
/* Returns the number of residues that it successfully found */
/* Requires:  chainid -- chain 1 ID (eg L) */
/*            resnumbers -- array of residue numbers (strings) */
/*            ptarray -- array of CResidue pointers */
/*            arraysize -- number of residues in input array */
/*------------------------------------------------------------------------*/
int CSystem::GetResidues(string chainid, string* resnumbers,
			 CResidue** ptarray, const int arraysize) {
  int     chainno = ChainNo(chainid);
  int     rnum, nres = 0;

  if (chainno < 0) {
    cout << "ERROR: could not match stored chain for ID '" << 
      chainid << "'" << endl;
    return(0);
  }

  /* Search for the residue then the atom with the correct labels */
  for(int i=0; i<arraysize; i++) {
    //    cout << resnumbers[i] << " in " << chainid << endl;
    /* find the residue in the chain */
    if (!chain[chainno].FindResWNum(resnumbers[i],rnum)) {
      cout << "WARNING: (system.GetResidues) could not find residue number '" 
	   << resnumbers[i] << "' in chain '" 
	   << chainid << "', skipping this residue" << endl;
      ptarray[i] = 0;
      continue;
    }

    /* get a pointer to the residue */
    ptarray[i] = chain[chainno].GetResiduePtr(rnum);
    nres += 1;
  }

  return(nres);
}

/*------------------------------------------------------------------------*/
/* Fill an array of integers correponding to a list of residue */
/* "numbers" supplied. */
/* Returns the number of residues that it successfully found */
/* Requires:  chainid -- chain 1 ID (eg L) */
/*            resnumbers -- array of residue numbers (strings) */
/*            ptarray -- array of CResidue pointers */
/*            arraysize -- number of residues in input array */
/*------------------------------------------------------------------------*/
int CSystem::GetResidueIndices(string chainid, string* resnumbers,
			       int* array, const int arraysize) {
  int     chainno = ChainNo(chainid);
  int     rnum, nres = 0;

  if (chainno < 0) {
    cout << "ERROR: could not match stored chain for ID '" << 
      chainid << "'" << endl;
    return(0);
  }

  /* Search for the residue then the atom with the correct labels */
  for(int i=0; i<arraysize; i++) {
    //    cout << resnumbers[i] << " in " << chainid << endl;
    /* find the residue in the chain */
    if (!chain[chainno].FindResWNum(resnumbers[i],rnum)) {
      cout << "WARNING: (system.GetResidues) could not find residue number '" 
	   << resnumbers[i] << "' in chain '" 
	   << chainid << "', skipping this residue" << endl;
      array[i] = 0;
      continue;
    }

    /* get a pointer to the residue */
    array[i] = rnum;
    nres += 1;
  }

  return(nres);
}

/*-------------------------------------------------------------------------*/
/* Copies one system structure into another, reducing the complexity of */
/* specified chains to only Calpha and Cbeta atoms */
/* Requires:  old -- old structure to copy from */
/*            chains -- comma-separated list of chain ID's to reduce */
/*            atomtypes -- array of atom types to copy */
/*            ntypes -- number of atom types to copy */
/*------------------------------------------------------------------------*/
void CSystem::CopyReduce(const CSystem& old, string chains, 
			 string* atomtypes, int ntypes) {
  bool    reduce;
  int     nreduced;
  string  chainlist[100];

  nreduced = Split(chains,',',chainlist,100);

  SizeNChains(old.nchains);
  iseed = old.iseed;
  temperature = old.temperature;

  for(int i=0; i<nchains; i++) {    
    reduce = false;
    //    cout << "copying chain '" << old.chain[i].ChainID() << "'" << endl;
    for(int j=0; j<nreduced; j++) {
      if (old.chain[i].ChainID() == chainlist[j]) {reduce = true; break;}
    }
    if (reduce) {
      chain[i].CopyReduce(old.chain[i],atomtypes,ntypes);
    } else {
      chain[i] = old.chain[i];
    }
  }

  /* setup the coordinate pointers to individual atoms */
  PopulateArrays();

  /* calculate the center of mass */
  CalcCOM();
}

/*-------------------------------------------------------------------------*/
/* Copies a subset of one system structure into another, optionally */
/* reducing the complexity to given atom types */
/* Requires:  old -- old structure to copy from */
/*            subset -- system subset definition */
/*            chains -- comma-separated list of chain ID's to reduce */
/*            atomtypes -- array of atom types to copy */
/*            ntypes -- number of atom types to copy */
/*------------------------------------------------------------------------*/
void CSystem::CopySubset(const CSystem& old, CSubset& subset, string chains,
			 string* atomtypes, int ntypes) {
  bool    reduce;
  int     nreduced, nnew, chainno = 0;
  string  chainlist[100];

  nreduced = Split(chains,',',chainlist,100);

  /* determine the new number of chains (may not copy some of them) */
  nnew = 0;
  for(int i=0; i<old.nchains; i++) {    
    if (subset.InSubset(old.chain[i].ChainID())) {nnew += 1;}
  }

  SizeNChains(nnew);
  iseed = old.iseed;
  temperature = old.temperature;

  /* copy relevant chains */
  for(int i=0; i<old.nchains; i++) {    
    if (subset.InSubset(old.chain[i].ChainID())) {

      /* determine if it should be reduced or just copied */
      reduce = false;
      for(int j=0; j<nreduced; j++) {
	if (old.chain[i].ChainID() == chainlist[j]) {reduce = true; break;}
      }

      if (reduce) {
	chain[chainno].CopySubset(old.chain[i],subset,atomtypes,ntypes);
      } else {
	chain[chainno].CopySubset(old.chain[i],subset,atomtypes,0);
      }
      chainno += 1;

    }
  }

  /* setup the coordinate pointers to individual atoms */
  PopulateArrays();

  /* calculate the center of mass */
  CalcCOM();
}

/*-------------------------------------------------------------------------*/
/* Copies a chain-based subset of one system structure into another */
/* Requires:  old -- old structure to copy from */
/*            subset -- system subset definition */
/*------------------------------------------------------------------------*/
void CSystem::CopySubset(const CSystem& old, CSubset& subset) {
  int     nnew, chainno = 0;

  /* determine the new number of chains (may not copy some of them) */
  nnew = 0;
  for(int i=0; i<old.nchains; i++) {    
    if (subset.InSubset(old.chain[i].ChainID())) {nnew += 1;}
  }

  SizeNChains(nnew);
  iseed = old.iseed;
  temperature = old.temperature;

  /* copy relevant chains */
  for(int i=0; i<old.nchains; i++) {    
    if (subset.InSubset(old.chain[i].ChainID())) {
      chain[chainno].CopySubset(old.chain[i],subset);
      chainno += 1;
    }
  }

  /* setup the coordinate pointers to individual atoms */
  PopulateArrays();

  /* calculate the center of mass */
  CalcCOM();
}

/*-------------------------------------------------------------------------*/
/* Creates a new system with sidechain centroids from the old system for */
/* those residues that are specified in the subset */
/* Requires:  old -- old structure to copy from */
/*            subset -- system subset definition */
/*------------------------------------------------------------------------*/
void CSystem::CopyCentroids(const CSystem& old, CSubset& subset) {
  int     nnew, chainno = 0;

  /* determine the new number of chains (may not copy some of them) */
  nnew = 0;
  for(int i=0; i<old.nchains; i++) {    
    if (subset.InSubset(old.chain[i].ChainID())) {nnew += 1;}
  }

  SizeNChains(nnew);
  iseed = old.iseed;
  temperature = old.temperature;

  /* copy relevant chains */
  for(int i=0; i<old.nchains; i++) {    
    if (subset.InSubset(old.chain[i].ChainID())) {
      chain[chainno].CopyCentroids(old.chain[i],subset);
      chainno += 1;
    }
  }

  /* setup the coordinate pointers to individual atoms */
  PopulateArrays();

  /* calculate the center of mass */
  CalcCOM();
}

/*-------------------------------------------------------------------------*/
/* Calculate squared distances between like residue atoms in two systems */
/* for each residue in a subset definition */
/* Requires:  two -- second system */
/*            subset -- system subset definition */
/*            atomtype -- atom type */
/*            distsq -- array of squared distances, MUST BE PRE-SIZED! */
/*------------------------------------------------------------------------*/
void CSystem::DistSq(CSystem& two, CSubset& subset, 
		     const string atomtype, double* distsq) {
  int       rnum, tempno, chainno;
  int       nres = subset.TotalResidues();
  string    resno, chainID;
  CVector3  coord1, coord2;

  for(int i=0;i<nres;i++) {
    distsq[i] = 0.0; 

    /* Get chain and residue number */
    tempno = subset.Num2ID(i,chainID,resno);
    //    cout << "calculating distance for " << chainID << ":" << resno << endl;
    chainno = ChainNo(chainID);
    if (!chain[chainno].FindResWNum(resno,rnum)) {
      cout << "WARNING: (system.DistSq) could not find residue '" 
	   << resno << "' in chain '" << chainID << "'" << endl;
      continue;
    }

    /* get coordinate of residue atom in FIRST system */
    if (! chain[chainno].AtomRWType(rnum,atomtype,coord1)) {
      //      cout << "WARNING: chain " << chainID << " residue " 
      //	   << resno << " (" << rnum << ")"
      //	   << " is missing atom type '" << atomtype << "'" << endl;
      continue;
    }

    /* get coordinate of residue atom in SECOND system */
    chainno = two.ChainNo(chainID);
    if (!two.chain[chainno].FindResWNum(resno,rnum)) {
      cout << "WARNING: (system.DistSq) could not find residue '" 
	   << resno << "' in chain '" << chainID << "' of sys2" << endl;
      continue;
    }
    if (! two.chain[chainno].AtomRWType(rnum,atomtype,coord2)) {
      //      cout << "WARNING: chain (sys 2) " << chainID << " residue " 
      //	   << resno << " (" << rnum << ")"
      //	   << " is missing atom type '" << atomtype << "'" << endl;
      continue;
    }

    /* calculate the squared distance */
    distsq[i] = coord1.DistSq(coord2);
    //    cout << i << " 1: " << coord1 << endl;
    //    cout << i << " 2: " << coord2 << endl;
    //    cout << i << "  " << distsq[i] << endl;
  }
}

/*-------------------------------------------------------------------------*/
/* Calculate distances between like residue atoms in two systems */
/* for each residue in a subset definition */
/* Requires:  two -- second system */
/*            subset -- system subset definition */
/*            atomtype -- atom type */
/*            dist -- array of distances, MUST BE PRE-SIZED! */
/*------------------------------------------------------------------------*/
void CSystem::Dist(CSystem& two, CSubset& subset, 
		   const string atomtype, double* dist) {
  int       nres = subset.TotalResidues();
  
  DistSq(two,subset,atomtype,dist);

  /* calculate the distance */
  for(int i=0;i<nres;i++) {
    dist[i] = sqrt(dist[i]);
  }
}

/*-------------------------------------------------------------------------*/
/* Calculate the RMSD (root mean square deviation) between two systems */
/* for each residue in a subset definition */
/* Requires:  two -- second system */
/*            subset -- system subset definition */
/*            atomtype -- atom type */
/*------------------------------------------------------------------------*/
double CSystem::RMSD(CSystem& two, CSubset& subset, const string atomtype) {
  int       nres = subset.TotalResidues();
  double    dist2sum = 0.0, rmsd;
  double*   distsq;

  distsq = new double[nres];
  DistSq(two,subset,atomtype,distsq);

  /* calculate the sum of squared distances */
  for(int i=0;i<nres;i++) {
    dist2sum += distsq[i];
    //    cout << i << "  " << distsq[i] << "  " << dist2sum << endl;
  }
  rmsd = sqrt(dist2sum/nres);
  delete[] distsq;
  return(rmsd);
}

/*-------------------------------------------------------------------------*/
/* Check the subset definition to make sure all the residues exist */
/* Requires:  subset -- system subset definition */
/*------------------------------------------------------------------------*/
bool CSystem::CheckSubset(CSubset& subset) {
  bool      problem = false;
  string    chainid,*residues;
  CResidue  *residue;

  for(int c=0; c<subset.NChainTypes(); c++) {
    chainid = subset.ChainID(c);
    residues = subset.GetResidues(chainid);
    for(int r=0; r<subset.NResidues(chainid); r++) {
      //      cout << chainid << residues[r] << endl;
      if (!GetResidue(chainid,residues[r],residue)) {
	cout << "ERROR: (system.CheckSubset) Could not find " 
	     << chainid << residues[r] << endl;
	problem = true;
      } else {
	//	cout << "found it" << endl;	
      }
    }
  }

  if (problem) {return(false);}
  return(true);
}

/*-------------------------------------------------------------------------*/
/* Change residue types for selected residues in system */
/* NOTE: for now, it just changes the types of all residues */
/* Requires:  newtype -- new 3-letter residue type code */
/*------------------------------------------------------------------------*/
void CSystem::Mutate(const string newtype) {
  int r,c,nres;

  for(c=0; c<nchains; c++) {  
    nres = chain[c].NRes();
    for(r=0; r<nres; r++) {
      chain[c].SetResType(r,newtype);
    }
  }
}

/*--------------------------------------------------------------------*/
/* Returns the symbol of a residue given the residue "number" */
/* and the chain ID. */
/* Requires:  chainid -- chain identifier */
/*            resno -- residue "number" */
/*--------------------------------------------------------------------*/
string CSystem::ResSymbol(const string chainid, const string resno) {
  int chainno = ChainNo(chainid);
  string symbol;

  if (chainno < 0) {
    cout << "ERROR: could not find chain ID '" 
	 << chainid << "' in system" << endl;
    exit(-1);
  }
  if (!chain[chainno].ResSymbolWNum(resno,symbol)) {
    cout << "ERROR: (system.ResSymbol) could not find residue number '" 
	 << resno << "' in chain '" << chainid << "' of system" << endl;
    exit(-1);
  }
  return(symbol);
}

/*-------------------------------------------------------------------------*/
/* Recalculate the Center of Mass (COM) of the individual chains and the */
/* system as a whole.  Calculation of the masses is included here. */
/*------------------------------------------------------------------------*/
CVector3 CSystem::CalcCOM() {
  double     chainmass[nchains];
  CVector3   vec(0.0), chaincom[nchains];

  mass = 0.0;
  for(int i=0; i<nchains; i++) {    
    //    cout << "calculating COM for chain " << chain[i].ChainID() << endl;
    chaincom[i] = chain[i].CalcCOM();
    chainmass[i] = chain[i].Mass();
    vec += chainmass[i]*chaincom[i];
    mass += chainmass[i];
  }
  if ((mass <= 0.0)&&(NAtoms() > 0)) {
    cout << "ERROR: (system.CalcCOM) system mass is zero, something failed" 
	 << endl;
    exit(-1);
  }
  com = vec/mass;
  return(com);
}

/*-------------------------------------------------------------------------*/
/* Calculate the Center of Mass (COM) of a set of individual chains */
/* using the stored COM's in the chain classes*/
/* Requires:  chainlist -- string array of chain ID's */
/*            nselected -- number of chains in list */
/*------------------------------------------------------------------------*/
CVector3 CSystem::COM(const string* chainlist, const int nselected) {
  double     totmass;
  CVector3   vec(0.0);

  totmass = 0.0;
  for(int i=0; i<nchains; i++) {    
    vec += chain[i].Mass()*chain[i].COM();
    totmass += chain[i].Mass();
  }
  mass = totmass;
  if (mass <= 0.0) {
    cout << "ERROR: (system.COM) system mass is zero, something failed" << endl;
    exit(-1);
  }
  com = vec/totmass;
  return(com);
}

/*------------------------------------------------------------------------*/
/* Write system to a pdb file */
/* Requires:  file -- filename for configuration */
/*------------------------------------------------------------------------*/
void CSystem::WritePDB(char* file) {
  ofstream pdbfile;

  /* Open the pdb file for output*/
  pdbfile.open(file, ios::out);
  if(!pdbfile) {
    cout << endl << " Error: cannot open file '" << file << 
      "'" << endl;
    exit(-1);
  }
  //  printf("Writing system configuration to pdb file '%s'\n",file);

  WritePDB(pdbfile);

  pdbfile.close();
}

/*------------------------------------------------------------------------*/
/* Write system to a pdb file with a comment as a REMARK record */
/* Requires:  file -- filename for configuration */
/*            comment -- string for the comment line */
/*------------------------------------------------------------------------*/
void CSystem::WritePDB(char* file, const string comment) {
  ofstream pdbfile;

  /* Open the pdb file for output*/
  pdbfile.open(file, ios::out);
  if(!pdbfile) {
    cout << endl << " Error: cannot open file '" << file << 
      "'" << endl;
    exit(-1);
  }
  //  printf("Writing system configuration to pdb file '%s'\n",file);

  pdbfile << "REMARK  " << comment << endl;
  WritePDB(pdbfile);

  pdbfile.close();
}

/*------------------------------------------------------------------------*/
/* Write system to a pdb file */
/* Requires:  pdbfile -- output stream */
/*------------------------------------------------------------------------*/
void CSystem::WritePDB(ofstream& pdbfile) {

  /* Write each chain to pdb formatted file */
  for(int c=0; c<nchains; c++) {  
    chain[c].WritePDB(pdbfile);
    pdbfile << "TER" << endl;
  }

  pdbfile << "END" << endl;
}

/*------------------------------------------------------------------------*/
/* Write system snapshot to a pdb file with a comment as a REMARK record */
/* Requires:  comment -- string for the comment line */
/*------------------------------------------------------------------------*/
void CSystem::WritePDBSnapshot(const string comment) {
  string   file="snapshot.pdb";
  ofstream pdbfile;

  /* Open the pdb file for output*/
  pdbfile.open(file.c_str(), ios::out);
  if(!pdbfile) {
    cout << endl << " Error: cannot open file '" << file << 
      "'" << endl;
    exit(-1);
  }
  //  printf("Writing system configuration to pdb file '%s'\n",file.c_str());

  pdbfile << "REMARK  iseed=" << Int2String(GetSeed()) << endl;
  pdbfile << "REMARK  " << comment << endl;
  WritePDB(pdbfile);

  pdbfile.close();
}

/*------------------------------------------------------------------------*/
/* Write system snapshot to a pdb file */
/*------------------------------------------------------------------------*/
void CSystem::WritePDBSnapshot() {
  string   file="snapshot.pdb";
  ofstream pdbfile;

  /* Open the pdb file for output*/
  pdbfile.open(file.c_str(), ios::out);
  if(!pdbfile) {
    cout << endl << " Error: cannot open file '" << file << 
      "'" << endl;
    exit(-1);
  }
  //  printf("Writing system configuration to pdb file '%s'\n",file.c_str());

  pdbfile << "REMARK  iseed=" << Int2String(GetSeed()) << endl;
  WritePDB(pdbfile);

  pdbfile.close();
}

/*------------------------------------------------------------------------*/
/* Write system to an xyz file */
/* Requires:  file -- filename for configuration */
/*            comment -- string for the comment line */
/*------------------------------------------------------------------------*/
void CSystem::WriteXYZ(char* file, string comment) {
  ofstream xyzfile;

  /* Open the xyz file for output*/
  xyzfile.open(file, ios::out);
  if(!xyzfile) {
    cout << endl << " Error: cannot open file '" << file << 
      "'" << endl;
    exit(-1);
  }
  //  printf("Writing system configuration to xyz file '%s'\n",file);

  xyzfile << NAtoms() << endl;
  xyzfile << comment << endl;

  /* Write each chain to xyz formatted file */
  for(int c=0; c<nchains; c++) {  
    chain[c].WriteXYZ(xyzfile);
  }

  xyzfile.close();
}

/*------------------------------------------------------------------------*/
/* Write system to an open xyz file as a single frame */
/* Requires:  xyzfile -- open file for configuration */
/*            comment -- string for the comment line */
/*------------------------------------------------------------------------*/
void CSystem::WriteXYZ(ofstream& xyzfile, string comment) {

  xyzfile << NAtoms() << endl;
  xyzfile << comment << endl;

  /* Write each chain to xyz formatted file */
  for(int c=0; c<nchains; c++) {  
    chain[c].WriteXYZ(xyzfile);
  }
}

/*------------------------------------------------------------------------*/
/* Write a pair of profit script files, one for fitting the antibody */
/* portion and one for the antigen portion */
/* Requires:  abdefn -- comma-separated list of antibody chain ID's */
/*            agdefn -- comma-separated list of antigen chain ID's */
/*------------------------------------------------------------------------*/
void CSystem::WriteProfitScripts(const string abdefn, const string agdefn) {
  int       nchains,chainno;
  string    filename,fields[100];
  ofstream  file;

  /* Open the antibody script file for output*/
  filename = "profit_antibody_script";
  file.open(filename.c_str(), ios::out);
  if(!file) {
    cout << endl << " Error: cannot open file '" << filename << 
      "'" << endl;
    exit(-1);
  }
  cout << "Creating profit script for antibody chains (" 
       << abdefn << "):" << endl;
  nchains = Split(abdefn,',',fields,100);

  /* write antibody script */
  file << "atom ca" << endl;
  for(int c=0; c<nchains; c++) {  
    if ((chainno = ChainNo(fields[c])) < 0) {
      cout << "ERROR: could not create profit script because chain ID '"
	   << fields[c] << "' could not be found" << endl;
      file.close();
      return;
    }
    for(int i=0; i<chain[chainno].NRes(); i++) {
      file << "zone " << fields[c] << chain[chainno].ResNo(i) 
	   << "-" << fields[c] << chain[chainno].ResNo(i);
      file << ":" << fields[c] << chain[chainno].ResNo(i) 
	   << "-" << fields[c] << chain[chainno].ResNo(i) << endl;
    }
  }
  file << "ignore" << endl;
  file << "fit" << endl;
  file << "write fittedAB.pdb" << endl;
  file << "quit" << endl;
  file.close();

  /* Open the antibody script file for output*/
  filename = "profit_antigen_script";
  file.open(filename.c_str(), ios::out);
  if(!file) {
    cout << endl << " Error: cannot open file '" << filename << 
      "'" << endl;
    exit(-1);
  }
  cout << "Creating profit script for antigen chain(s) (" 
       << agdefn << "):" << endl;
  nchains = Split(agdefn,',',fields,100);

  /* write antibody script */
  file << "atom ca" << endl;
  for(int c=0; c<nchains; c++) {  
    if ((chainno = ChainNo(fields[c])) < 0) {
      cout << "ERROR: could not create profit script because chain ID '"
	   << fields[c] << "' could not be found" << endl;
      file.close();
      return;
    }
    for(int i=0; i<chain[chainno].NRes(); i++) {
      file << "zone " << fields[c] << chain[chainno].ResNo(i) 
	   << "-" << fields[c] << chain[chainno].ResNo(i);
      file << ":" << fields[c] << chain[chainno].ResNo(i) 
	   << "-" << fields[c] << chain[chainno].ResNo(i) << endl;
    }
  }
  file << "ignore" << endl;
  file << "fit" << endl;
  file << "write fittedAG.pdb" << endl;
  file << "quit" << endl;
  file.close();
}

/*-------------------------------------------------------------------------*/
/* Generate a subset class using all the residues in the system */
/* Requires:  subset -- system subset */
/*------------------------------------------------------------------------*/
CSubset CSystem::GenSubset() {
  int       r, c, nres, nlines = 0, maxlines = 10000;
  string    line[maxlines], chainid;
  CSubset   subset;

  for(c=0; c<nchains; c++) {  
    chainid = chain[c].ChainID();
    nres = chain[c].NRes();
    for(r=0; r<nres; r++) {
      //      cout << chainid << ":" << chain[c].ResNo(r) << endl;
      line[nlines] = Int2String(nlines) + "  " + chainid 
	+ "  " + chain[c].ResNo(r);
      //      cout << line[nlines] << endl;
      nlines += 1;
    }
  }
  subset.Init(line,nlines,false);

  return subset;
}

/*-------------------------------------------------------------------------*/
/* Generate a subset class using all the residues in the specified */
/* chains of the system */
/* Requires:  subset -- system subset */
/*            chains -- comma-separated chain identifiers */
/*------------------------------------------------------------------------*/
CSubset CSystem::GenSubset(const string chains) {
  int       r, c, k, nres, nlines = 0, maxlines = 10000, nkeep;
  bool      keepit;
  string    line[maxlines], chainid, chainlist[100];
  CSubset   subset;

  nkeep = Split(chains,',',chainlist,100);

  for(c=0; c<nchains; c++) {  
    chainid = chain[c].ChainID();
    keepit = false;
    for(k=0; k<nkeep; k++) {  
      if (chainid == chainlist[k]) {keepit = true;}
    }
    if (!keepit) {continue;}
    nres = chain[c].NRes();
    for(r=0; r<nres; r++) {
      //      cout << chainid << ":" << chain[c].ResNo(r) << endl;
      line[nlines] = Int2String(nlines) + "  " + chainid 
	+ "  " + chain[c].ResNo(r);
      //      cout << line[nlines] << endl;
      nlines += 1;
    }
  }
  subset.Init(line,nlines,false);

  return subset;
}

/*--------------------------------------------------------------------*/
/* A display routine to dump system characteristics to screen */
/* Requires:  indent -- number of spaces to indent */
/*--------------------------------------------------------------------*/
void CSystem::Display(int indent) {
  int   c;
  char  spacing[100];

  for(int i=0;i<indent;i++) {spacing[i] = ' ';}
  spacing[indent] = '\0';

  cout << spacing << "System characteristics:" << endl;
  cout << spacing << "system temperature (K): " << Num(temperature) << endl;
  cout << spacing << "current random number seed: " << iseed << endl;
  cout << spacing << "number of chains: " << nchains << endl;

  for(c=0; c<nchains; c++) {  
    cout << spacing << " chain " << c << "  (" << chain[c].NAtoms() <<
      " atoms): " << endl ;
    chain[c].Display(indent+1);
  }
}
