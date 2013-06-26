/* This file provides routines for performing lookup table
   model energy and force (?) evaluations.
   The interaction table file format should be:
   [number of entries]
   [separation distance]  [potential value]
   ...
*/

#include <string>
#include <iostream>
#include <fstream>
#include "include/config.h"
#include "include/pairint.h"

using std::string;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::cout;
using std::endl;

/*------------------------------------------------------------------------*/
/* Initializes from a file */
/* Requires:  filename -- filename from which to read configuration */
/*            label -- label for the interaction */
/*            verbose -- flag indicating output to be written to screen */
/*            warn -- flag indicating if warnings should be given */
/*------------------------------------------------------------------------*/
void CPairInt::Init(string filename, string label, bool verbose, bool warn) {
  int      lineno = 0, nfields, nwarn = 0;
  bool     initialzero = true;
  double   lastsep, sep, value;
  char     buffer[1000];
  string   line, field[10];
  ifstream file;

  id = label;

  /* Get the number of non-comment lines in the interaction file */
  nentries = NLines(filename.c_str());
  if (nentries < 0) {
    cout << " ERROR: could not read interaction table file '" 
	 << filename << "'" << endl;
    exit(-1);
  }
  table = new double[nentries];
  //  cout << filename << "' with " << nentries << " lines" << endl;

  /* Open the interaction file */
  file.open(filename.c_str(), ios::in);
  if(!file) {
    cout << endl << " ERROR: cannot open interaction table file '" << 
      filename << "'" << endl;
    exit(-1);
  }

  while (!file.eof()) {
    /* read the line */
    file.getline(buffer,1000);
    line = buffer;

    /* clean line, skip if there is no content */
    StripComment(line);
    if (line.length() == 0) {continue;}

    /* clean and check line */
    nfields = Split(line,' ',field,10);
    if ((warn)&&(nfields > 2)) {
      if (verbose) {
	if (nwarn == 0) {
	  cout << "WARNING: line " << lineno + 1 << " in file '" << filename << 
	    "' contains more than two fields" << endl;
	}
	nwarn += 1;
      }
    } else if (nfields < 2) {
      cout << "ERROR: line " << lineno + 1 << " in file '" << filename << 
	"' contains more less than than two fields" << endl;
      exit(-1);
    }

    /* extract information */
    sep = ToDbl(field[0]);
    value = ToDbl(field[1]);

    /* skip the first zero value lines (undefined) */
    if ((initialzero)&&(value == 0.00)) {continue;}
    initialzero = false; 

    /* record information */
    if (lineno == 0) {
      firstsep = sep;
    } else if (lineno == 1) {  
      sepinc = sep - firstsep;
    } else {
      if (fabs(sep - lastsep - sepinc) > 1.0e-4) {
	cout << "ERROR: line " << lineno + 1 << " in file '" << filename << 
	  "' does not maintain separation increment " << endl;
	cout << " previous separation: " << lastsep << endl;
	cout << " current separation:  " << sep << endl;
	cout << " established separation increment:  " << sepinc << endl;
	cout << " new separation increment:  " << sep - lastsep - sepinc << endl;
	exit(-1);
      }
    }
    lastsep = sep;
    table[lineno] = value;
    lineno += 1;
  }
  finalsep = sep;

  /* set the number of entries recorded (can be less than initialized size) */
  nentries = lineno;

  if (nwarn > 0) {
    cout << " warning repeated " << nwarn << " times" << endl;
  }
  file.close();
}

/*------------------------------------------------------------------------*/
/* Checks a g(r) function for negative values and values above a given */
/* value.  Returns the number of both types of values and warns if desired */
/* Requires:  max -- maximum value above which value is flagged */
/*            nneg -- number of negative values */
/*            nover -- number of values over maximum */
/*            warn -- flag indicating if warnings should be given */
/*------------------------------------------------------------------------*/
void CPairInt::Check(const double max, int& nneg, int& nover, bool warn) {
  nneg = 0;
  nover = 0;

  for(int i=0; i<nentries; i++) {
    if (table[i] < 0.0) {nneg += 1;}
    if (table[i] > max) {nover += 1;}
  }  

  if (warn) {
    if (nneg > 0) {
      cout << "WARNING: pair interaction '" << id << "' has " << nneg 
	   << " values" << endl;
    }
    if (nover > 0) {
      cout << "WARNING: pair interaction '" << id << "' has " << nover
	   << " values greater than " << Num(max) << endl;
    }
  }
}

/*------------------------------------------------------------------------*/
/* A display routine for the interactions */
/* Requires:  indent -- number of spaces to indent */
/*------------------------------------------------------------------------*/
void CPairInt::Display(int indent) {
  int     i;
  char    spacing[100];
  double  inc;

  for(i=0;i<indent;i++) {spacing[i] = ' ';}
  spacing[indent] = '\0';

  cout << spacing << "Lookup table for interaction: " << id << endl;
  cout << spacing << "Nentries= " << nentries 
       << "  separation increment= " << sepinc << endl;

  for(int i=0; i<nentries; i++) {
    inc = firstsep + i*sepinc;
    cout << i << "  " << spacing << Num(inc) << "   " << Num(table[i]) << endl;
  }
}


