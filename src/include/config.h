#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <string>

#define MAXNATOMS       150000
#define MINEXP         -30.0
#define MAXEXP          30.0

const double calToj = 4.184;
const double kjmole_kb = 8.314e-3;             /* kJ/mole/K */
const double kcalmole_kb = kjmole_kb/calToj;   /* kcal/mole/K */
const double PI = 3.1415927;
const double TWOPI = 2.0*PI;
const double DEG2RAD = PI/180.0;
const double RAD2DEG = 180.0/PI;

/* Sets the 3-letter amino acid codes */
const int NAATypes = 24;
const std::string AAType[NAATypes] = 
  {"ALA",
   "ARG",
   "ASN",
   "ASP",
   "CYS",
   "GLN",
   "GLU",
   "GLY",
   "HIS",
   "HSD",
   "HSE",
   "ILE",
   "LEU",
   "LYS",
   "MET",
   "PHE",
   "PRO",
   "SER",
   "THR",
   "TRP",
   "TYR",
   "VAL",
   "UNK",
   "HET"};

/* Sets the 1-letter amino acid symbols */
const std::string AASymbol[NAATypes] = 
  {"A",
   "R",
   "N",
   "D",
   "C",
   "Q",
   "E",
   "G",
   "H",
   "H",
   "H",
   "I",
   "L",
   "K",
   "M",
   "F",
   "P",
   "S",
   "T",
   "W",
   "Y",
   "V",
   "U",
   "Z"};

/* Sets the number of amino acid interaction types */
const int NIntResTypes = 22;
const std::string AAIntType[NIntResTypes] = 
  {"A",
   "R",
   "N",
   "D",
   "C",
   "Q",
   "E",
   "G",
   "H",
   "I",
   "L",
   "K",
   "M",
   "F",
   "P",
   "S",
   "T",
   "W",
   "Y",
   "V",
   "U",
   "Z"};

/* Sets certain residue type indices */
const int AlanineResType = 0;
const int GlycineResType = 7;
const int UnknownResType = 20;  /* must be first non-standard amino acid! */

/* Sets the recognized elements */
const int NElements = 11;
const std::string Element[NElements] = 
  {"H",
   "C",
   "N",
   "O",
   "NA",
   "MG",
   "P",
   "S",
   "MN",
   "FE",
   "ZN"};

/* Sets the mass of the recognized elements */
const double ElementMass[NElements] = 
  {1.00794,     
   12.011,
   14.00674,
   15.994,
   22.989768,
   24.3050,
   30.973762,
   32.066,
   54.93805,
   55.847,
   65.39
  };

double gasdev(long& idum);
double ran2(long& idum);

void flaglist(char* str);

/* global variables (used for debugging only!), accessable anywhere */
extern bool debugflag;
extern int  iteration;

#endif



