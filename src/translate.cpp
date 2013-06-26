/* This file provides routines for performing sub-system translations */

#include <string>
#include <iostream>
#include <fstream>
#include "include/config.h"
#include "include/translate.h"

using std::string;
using std::cout;
using std::endl;

/*------------------------------------------------------------------------*/
/* Initializes from a file */
/* Requires: chains -- comma-separated list of chain ID's to move */
/*------------------------------------------------------------------------*/
void CTranslate::Init(const string chains) {
  string     temp[100];

  stepsize = 1.0;

  nmobile = Split(chains,',',temp,100);
  mobile = new string[nmobile];
  for(int i=0; i<nmobile; i++) {mobile[i] = temp[i];}
}

/*------------------------------------------------------------------------*/
/* Performs a rotation on part of the system */
/* Requires:  distance -- separation distance */
/*------------------------------------------------------------------------*/
bool CTranslate::Move(CSystem& sys, const CForcefield& ff) {
  int         c,a,chainno;
  CVector3**  coord;

  /* generate the displacements at random */
  dx = stepsize*(sys.GetRan() - 0.5);
  dy = stepsize*(sys.GetRan() - 0.5);
  dz = stepsize*(sys.GetRan() - 0.5);
  dispvec = CVector3(dx,dy,dz);

  /* translate the specified chains */
  for(c=0; c<nmobile; c++) {
    chainno = sys.ChainNo(mobile[c]);
    coord = sys.GetChainR(chainno);
    for(a=0; a<sys.NAtoms(chainno); a++) {
      *(coord[a]) += dispvec;
    }
  }

  return(true);
}

/*--------------------------------------------------------------------*/
/* A one-line output of the state of the class */
/*--------------------------------------------------------------------*/
string CTranslate::Disp() {
  char      cstring[80];
  string    line;

  sprintf(cstring,"TRANSLATE (stepsize = %.3f) (%.2f %.2f %.2f)",
	  stepsize,dx,dy,dz);
  line = cstring;
  return(line);
}

/*--------------------------------------------------------------------*/
/* A two-line output of the state of the class */
/* Requires:  string1 -- label for move type */
/*            string2 -- stepsize for move type */
/*--------------------------------------------------------------------*/
void CTranslate::Disp(string& string1, string& string2) {
  string1 = "TRANSLATE";
  string2 = Num(stepsize);
}

/*------------------------------------------------------------------------*/
/* A display routine for the interactions */
/*------------------------------------------------------------------------*/
void CTranslate::Display() { 
  cout << "translate move" << endl;
  cout << "current displacements (dx,dy,dz): " 
       << dx << " " << dy << " " << dz << endl;
}
