// Functions for the subset class.  These routines calculate energies
// and/or forces for the system.

#include <iostream>
#include "include/config.h"
#include "include/subset.h"

using std::string;
using std::cout;
using std::ios;
using std::endl;

/*--------------------------------------------------------------------*/
/* Requires:  residuesfile -- name of file containing residues specs */
/*            verbose -- flag indicating output to be written to screen */
/*--------------------------------------------------------------------*/
CSubset::CSubset(const string residuesfile, bool verbose) {
  Init(residuesfile,verbose);
}

/*--------------------------------------------------------------------*/
/* Requires:  residuesfile -- name of file containing residues specs */
/*            chains -- comma-separated list of chain ID's */
/*            verbose -- flag indicating output to be written to screen */
/*--------------------------------------------------------------------*/
CSubset::CSubset(const string residuesfile, string chains, bool verbose) {
  Init(residuesfile,chains,verbose);
}

/*--------------------------------------------------------------------*/
/* Requires:  residues -- array of strings containing subset lines */
/*            nlines -- number of lines (strings) */
/*            verbose -- flag indicating output to be written to screen */
/* Line Format:*/
/* [number, ignored]  [chain ID]  [residue ID]  [4th+ column ignored] */
/*--------------------------------------------------------------------*/
CSubset::CSubset(string* residues, const int nlines, bool verbose) {
  Init(residues,nlines,verbose);
}

/*--------------------------------------------------------------------*/
/* Initializes the system subset definition from a file */
/* Requires:  residuesfile -- name of file containing residues specs */
/*            verbose -- flag indicating output to be written to screen */
/*--------------------------------------------------------------------*/
void CSubset::Init(const string residuesfile, bool verbose) {
  int        nlines;
  const int  maxlines = 10000;
  string     line[maxlines];

  if (verbose) {
    cout << " reading interacting residue definitions from file '" << 
  	residuesfile << "'" << endl;
  }
  nchaintypes = 0;
  nlines = ReadFile(residuesfile.c_str(),line,maxlines);
  if (nlines < 0) {exit(-1);}
  if (verbose) {
    cout << "  read " << nlines << " lines from file '" 
	 << residuesfile << "'" << endl;
  }

  Init(line,nlines,verbose);
}

/*--------------------------------------------------------------------*/
/* Initializes the system subset definition from a file */
/* Using only lines that correspond to the passed chain ID's */
/* Requires:  residuesfile -- name of file containing residues specs */
/*            chains -- comma-separated list of chain ID's */
/*            verbose -- flag indicating output to be written to screen */
/*--------------------------------------------------------------------*/
void CSubset::Init(const string residuesfile, string chains, bool verbose) {
  int        nlines, nkeep, nfields;
  const int  maxlines = 10000, arraysize=100;
  bool       keep;
  string     line[maxlines],chainlist[arraysize],field[arraysize];

  nkeep = Split(chains,',',chainlist,arraysize);

  if (verbose) {
    cout << " reading interacting residue definitions from file '" << 
  	residuesfile << "'" << endl;
  }
  nchaintypes = 0;
  nlines = ReadFile(residuesfile.c_str(),line,maxlines);
  if (nlines < 0) {exit(-1);}
  if (verbose) {
    cout << "  read " << nlines << " lines from file '" 
	 << residuesfile << "'" << endl;
  }

  /* blank the lines that don't contain chain ID's which must be kept */
  for(int i=0; i<nlines; i++) {
    nfields = Split(line[i],' ',field,arraysize);
    keep = false;
    for(int j=0; j<nkeep; j++) {
      if (field[1] == chainlist[j]) {keep = true; break;}
    }
    if (!keep) {line[i] = "";}
  }

  Init(line,nlines,verbose);
}

/*--------------------------------------------------------------------*/
/* Initializes the system subset definition from an array of lines */
/* Requires:  line -- array of strings containing residues specs */
/*            verbose -- flag indicating output to be written to screen */
/* Line Format:*/
/* [number, ignored]  [chain ID]  [residue ID]  [4th+ column ignored] */
/*--------------------------------------------------------------------*/
void CSubset::Init(string* line, const int nlines, bool verbose) {
  bool        newtype,skip;
  int         i,j, chainno, nfields, nkeep = 0;
  const int   arraysize = 100;
  string      chaintypes[arraysize], field[arraysize];

  /* count number of different chain types */
  nchaintypes = 0;
  for(i=0; i<nlines; i++) {
    if (line[i] == "") {continue;}   /* skip blank lines */
    if ((nfields = Split(line[i],' ',field,arraysize)) < 3) {
      cout << "ERROR: (subset.Init) line number " << i+1 
	   << "' contains less than the minimum 3 expected columns " << endl;
      exit(-1);
    }
    nkeep += 1;
    newtype = true;
    for(j=0; j<nchaintypes; j++) {
      if (field[1] == chaintypes[j]) {newtype = false; break;}
    }
    if (newtype) {chaintypes[nchaintypes] = field[1]; nchaintypes += 1;}
  }
  if (verbose) {
    cout << "  found " << nchaintypes << " interacting chain types " << endl;
  }

  /* size the arrays defining the interacting residues */
  intchain = new string[nchaintypes];
  for(i=0; i<nchaintypes; i++) {
    intchain[i] = chaintypes[i];
  }
  nresidues = new int[nchaintypes];
  intres = new string*[nchaintypes];
  totalresidues = 0;
  for(i=0; i<nchaintypes; i++) {
    nresidues[i] = 0;
    intres[i] = new string[nkeep];
  }

  /* build the arrays defining the interacting residues */
  maxnresidues = 0;
  for(i=0; i<nlines; i++) {
    if (line[i] == "") {continue;}   /* skip blank lines */
    if ((nfields = Split(line[i],' ',field,arraysize)) < 3) {
      cout << "ERROR: (subset.Init) line number " << i+1 
	   << "' contains less than the minimum 3 expected columns " << endl;
      exit(-1);
    }

    /* find the correct internal chain definition index */
    for(j=0; j<nchaintypes; j++) {
      if (field[1] == chaintypes[j]) {chainno = j; break;}
    }

    /* make sure this this chain:residue pair hasn't already been recorded */
    //    cout << "checking: " << line[i] << endl;
    skip = false;
    for(j=0; j<nresidues[chainno]; j++) {
      if (intres[chainno][j] == field[2]) {
	cout << "WARNING: (subset.init) residue " << field[1] << ":" 
	     << field[2] << " is already in subset, skipping" << endl;
	//	cout << line[i] << endl;
	skip = true;
	break;
      }
    }
    if (skip) {continue;}
    //    cout << "adding: " << line[i] << endl;

    /* add this residue number to the array for the appropriate chain type */
    intres[chainno][nresidues[chainno]] = field[2];
    //    cout << "added '" << intres[chainno][nresidues[chainno]] << 
    //      "' at position " << nresidues[chainno] << endl;
    nresidues[chainno] += 1;
    totalresidues += 1;
    if (maxnresidues < nresidues[chainno]) {maxnresidues = nresidues[chainno];}
  }

  if (verbose) {
    for(i=0; i<nchaintypes; i++) {
      cout << "  found " << nresidues[i] << " residues of chain type '" 
	   << intchain[i] << "'" << endl;
    }
  }

  /* set the interacting chain types in matrix form, all interactions on */
  /* default is all interactions on except intrachain */
  intmtx = new bool*[nchaintypes];
  for(i=0; i<nchaintypes; i++) {
    intmtx[i] = new bool[nchaintypes];    
    for(j=0; j<nchaintypes; j++) {
      if (i == j) {
	intmtx[i][j] = false;
      } else {
	intmtx[i][j] = true;
      }
    }
  }

  /* initialize the basic unknown information */
  nunknown = 0;
  chainunknown = new bool[nchaintypes];
  for(int i=0;i<nchaintypes;i++) {
    chainunknown[i] = false;
  }

  /* Initialize the arrays for the alternative sequence information */
  altseq = false;
  altres = new string*[nchaintypes];
  altresindx = new int*[nchaintypes];
  for(int i=0; i<nchaintypes; i++) {
    altres[i] = new string[nresidues[i]];
    altresindx[i] = new int[nresidues[i]];
    for(j=0; j<nresidues[i]; j++) {
      altres[i][j] = "Z";
      altresindx[i][j] = -1;
    }
  }
}

/*--------------------------------------------------------------------*/
/* Sets the unknown chain types.  These are the chains which will be */
/* assumed to have variable residue types, called 'ligand' in forcefield */
/* Will also set the chain-chain interactions to only known-unknown */
/* type interactions, if unkflag is true */
/* Requires:  identifiers -- comma-separated string of unknown ID's */
/*            unkflag -- if true, set eval to only known-unknown pairs */
/*--------------------------------------------------------------------*/
void CSubset::SetUnknown(const string identifiers, const bool unkflag) {
  int        index;
  const int  arraysize = 100;
  string     field[arraysize];

  if ((nunknown = Split(identifiers,',',field,arraysize)) <= 0) {
    cout << "WARNING: no chain identifiers for variable residue chains " <<
      "(SetUnknown) found in '" << identifiers << "'" << endl;
    //    return;
  }

  chainunknown = new bool[nchaintypes];
  unknownchain = new string[nunknown];
  for(int i=0;i<nunknown;i++) {
    unknownchain[i] = field[i];
  }

  for(int i=0;i<nchaintypes;i++) {
    chainunknown[i] = false;
    index = IndexInArray(intchain[i],unknownchain,nunknown);
    if (index >= 0) {chainunknown[i] = true;}
  }

  /* set the interacting chain types to be only known-unknown pairs */
  if (unkflag) {
    for(int i=0; i<nchaintypes; i++) {
      for(int j=0; j<nchaintypes; j++) {
	intmtx[i][j] = false;
	if ( ((chainunknown[i])&&(!chainunknown[j]))||
	     ((!chainunknown[i])&&(chainunknown[j])) ) {
	  intmtx[i][j] = true;	
	}
      }
    }
  }

}

/*--------------------------------------------------------------------*/
/* Translate a chain ID and residue number (1...N in subset) to a number */
/* Requires:  chainid -- ID of chain */
/*            num -- residue number in subset */
/*--------------------------------------------------------------------*/
int CSubset::ID2Num(const std::string chainid, const int num) {
  int resnum = 0;
  int index = ChainIndex(chainid);
  for(int i=0; i<index; i++) {
    resnum += nresidues[i];
  }
  return(resnum + num);
}

/*--------------------------------------------------------------------*/
/* Translate the number from ID2Num to a chain ID and residue number */
/* (1...N in subset) */
/* Requires:  num -- linear residue number */
/*            chainid -- ID of chain */
/*            resno -- real (PDB) residue "number" */
/*--------------------------------------------------------------------*/
int CSubset::Num2ID(const int num, std::string& chainid, std::string& resno) {
  int resnum = num;
  for(int i=0; i<nchaintypes; i++) {
    if (resnum < nresidues[i]) {
      chainid = intchain[i];
      resno = intres[i][resnum];
      break;
    }
    resnum -= nresidues[i];
  }
  return(resnum);
}

/*--------------------------------------------------------------------*/
/* Turn specific chain-chain interactions on/off */
/* Requires:  chainid1 -- ID of first chain */
/*            chainid2 -- ID of second chain */
/*            flag -- new value of interaction true/false flag */
/*--------------------------------------------------------------------*/
void CSubset::SetChainChainInt(const string chainid1, 
			       const string chainid2, bool flag) {
  int index1 = ChainIndex(chainid1);
  int index2 = ChainIndex(chainid2);
  if (index1 < 0) {
    cout << "WARNING: (subset.SetChainChainInt) " 
	 << "unable to find index for chain ID '" << chainid1 << "' " 
	 << "skipping specification" << endl;
    return;
  }
  if (index2 < 0) {
    cout << "WARNING: (subset.SetChainChainInt) " 
	 << "unable to find index for chain ID '" << chainid2 << "' " 
	 << "skipping specification" << endl;
    return;
  }
  intmtx[index1][index2] = flag;
  intmtx[index2][index1] = flag;
}

/*--------------------------------------------------------------------*/
/* Returns the stored alternative sequence of a specific chain */
/* Requires:  chainid -- string chain ID */
/*--------------------------------------------------------------------*/
string CSubset::ChainSeq(const string chainid) {
  int      index = ChainIndex(chainid);
  string   chainseq;

  for(int r=0;r<nresidues[index];r++) {
    chainseq += altres[index][r];
  }
  return(chainseq);
}

/*--------------------------------------------------------------------*/
/* A display routine to dump subset specifications to screen */
/* Requires:  indent -- number of spaces to indent */
/*--------------------------------------------------------------------*/
void CSubset::Display(int indent) {
  char    spacing[100];
  string  chainseq;

  for(int i=0;i<indent;i++) {spacing[i] = ' ';}
  spacing[indent] = '\0';

  cout << "Residues in the subset (" << totalresidues << ") " 
       << "(chain ID, residue number): " << endl;
  for(int c=0; c<nchaintypes; c++) {
    for(int r=0; r<nresidues[c]; r++) {
      cout << intchain[c] << "  " << intres[c][r] << endl;
    }
  }

  if (nunknown > 0) {
    cout << "variable (unknown) residue type chain ID's (" << nunknown << "): ";
    for(int i=0; i<nunknown; i++) {
      cout << unknownchain[i] << "(" << i << ") ";
    }
    cout << endl;
  }

  cout << "unknown flag for each chain ID: ";
  for(int i=0; i<nchaintypes; i++) {
    cout << chainunknown[i] << "(" << i << ") ";
  }
  cout << endl;

  cout << "interacting chain types: ";
  for(int i=0; i<nchaintypes; i++) {
    cout << intchain[i] << "(" << i << ") ";
  }
  cout << endl;

  cout << "chain-chain interaction matrix:" << endl;
  for(int i=0; i<nchaintypes; i++) {
    for(int j=0; j<nchaintypes; j++) {
      cout << intmtx[i][j] << " ";
    }
    cout << endl;
  }

  cout << "alternative sequence for each chain:" << endl;
  for(int i=0; i<nchaintypes; i++) {
    chainseq = ChainSeq(ChainID(i));
    cout << ChainID(i) << ": " << chainseq << endl;
  }
  
}
