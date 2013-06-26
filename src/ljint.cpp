/* This file provides routines for performing Lennard-Jones 
   energy and force (?) evaluations.  It's here for testing purposes.
*/

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "include/config.h"
#include "include/ljint.h"

using std::string;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::cout;
using std::endl;

/*------------------------------------------------------------------------*/
/* Initializes the full Lennard-Jones interaction */
/* Requires:  sig -- atomic radius */
/*            eps -- energy scaling parameter */
/*            label -- string label for interaction type */
/*------------------------------------------------------------------------*/
void CLJInt::Init(double sig, double eps, string label) {
  sigma = sig;
  epsilon = eps;
  id = label;
  locut = 0.1;
  hicut = 13.0;
  locut2 = locut*locut;
  hicut2 = hicut*hicut;

  /* calculate A and B */
  /* A = 4*eps*(sig^12), B = 4*eps*(sig^6) */
  B = 4.0*epsilon*pow(sigma,6)*kcalmole_kb;
  A = B*pow(sigma,6);
}

/*------------------------------------------------------------------------*/
/* Initializes the repulsive-only Lennard-Jones interaction */
/* Requires:  param -- repulsive coefficient */
/*            label -- string label for interaction type */
/*------------------------------------------------------------------------*/
void CLJInt::Init(double param, string label) {
  A = param;
  id = label;
  sigma = 0.0;
  epsilon = 0.0;
  locut = 0.1;
  hicut = 13.0;
  locut2 = locut*locut;
  hicut2 = hicut*hicut;
}

/*------------------------------------------------------------------------*/
/* Checks an initialization (not used) */
/*------------------------------------------------------------------------*/
void CLJInt::Check(const double max, int& nneg, int& nover, bool warn) {
  cout << "ERROR: cannot check LJ potential, routine is non-functional" << endl;
  exit(-1);
}

/*------------------------------------------------------------------------*/
/* A display routine for the interactions */
/* Requires:  indent -- number of spaces to indent */
/*------------------------------------------------------------------------*/
void CLJInt::Display(int indent) { 
  int     i;
  char    spacing[100];

  for(i=0;i<indent;i++) {spacing[i] = ' ';}
  spacing[indent] = '\0';

  cout << spacing << id;
}
