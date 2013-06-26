/* This file provides routines for performing hard-sphere repulsion
   energy and force (?) evaluations.  
*/

#include <string>
#include <iostream>
#include <fstream>
#include "include/config.h"
#include "include/hardsphereint.h"

using std::string;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::cout;
using std::endl;

/*------------------------------------------------------------------------*/
/* Initializes the full hard sphere interaction */
/* Requires:  eps -- energy scaling parameter */
/*            rad -- atomic radius */
/*            label -- string label for interaction type */
/*------------------------------------------------------------------------*/
void CHardSphere::Init(double eps, double rad, string label) {
  epsilon = eps;
  radius = rad;
  radius2 = rad*rad;
  id = label;
}

/*------------------------------------------------------------------------*/
/* Checks an initialization (not used) */
/*------------------------------------------------------------------------*/
void CHardSphere::Check(const double max, int& nneg, int& nover, bool warn) {
  cout << "ERROR: cannot check HS potential, routine is non-functional" << endl;
  exit(-1);
}

/*------------------------------------------------------------------------*/
/* A display routine for the interactions */
/* Requires:  indent -- number of spaces to indent */
/*------------------------------------------------------------------------*/
void CHardSphere::Display(int indent) { 
  int     i;
  char    spacing[100];

  for(i=0;i<indent;i++) {spacing[i] = ' ';}
  spacing[indent] = '\0';

  cout << spacing << id;
}
