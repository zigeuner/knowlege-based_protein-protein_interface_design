/* This file provides routines for performing sub-system rotations */

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include "include/config.h"
#include "include/rotate.h"

using std::string;
using std::cout;
using std::endl;

/*------------------------------------------------------------------------*/
/* Initializes from a file */
/* Requires: chains -- comma-separated list of chain ID's to move */
/*------------------------------------------------------------------------*/
void CRotate::Init(const string chains) {
  string     temp[100];

  stepsize = 2.0;   // degrees
  stepsize *= DEG2RAD;

  nmobile = Split(chains,',',temp,100);
  mobile = new string[nmobile];
  for(int i=0; i<nmobile; i++) {mobile[i] = temp[i];}
}

/*------------------------------------------------------------------------*/
/* Performs a rotation on part of the system */
/* Requires:  distance -- separation distance */
/*------------------------------------------------------------------------*/
bool CRotate::Move(CSystem& sys, const CForcefield& ff) {
  int         c,a,chainno;
  double      costheta;
  CVector3    com;
  CVector3**  coord;

  /* translate the specified chains to place their COM at origin */
  com = sys.COM(mobile,nmobile);
  for(c=0; c<nmobile; c++) {
    chainno = sys.ChainNo(mobile[c]);
    if (chainno < 0) {
      cout << " Error: could not match chain ID '" << mobile[c] << "'" << endl;
      exit(-1);
    }
    coord = sys.GetChainR(chainno);
    for(a=0; a<sys.NAtoms(chainno); a++) {
      *(coord[a]) -= com;
    }
  }

  /* generate the phi and psi angles at random, sample theta from cos dist. */
  costheta = stepsize*(2.0*sys.GetRan() - 1.0);
  //  costheta -= round(costheta/2.0)*2.0;
  costheta -= rint(costheta/2.0)*2.0;
  theta = acos(costheta) - PI/2;
  //  cout << costheta << "  " << theta << endl;
  phi = stepsize*(2.0*sys.GetRan() - 1.0);
  psi = stepsize*(2.0*sys.GetRan() - 1.0);
  rotmtx.ConstructEulerMtx(phi,theta,psi);

  /* apply rotation matrix to system subset */
  for(c=0; c<nmobile; c++) {
    chainno = sys.ChainNo(mobile[c]);
    if (chainno < 0) {
      cout << " Error: could not match chain ID '" << mobile[c] << "'" << endl;
      exit(-1);
    }
    coord = sys.GetChainR(chainno);
    for(a=0; a<sys.NAtoms(chainno); a++) {
      *(coord[a]) = rotmtx * *(coord[a]);
    }
  }

  /* translate the specified chains to place their COM at original position */
  for(c=0; c<nmobile; c++) {
    chainno = sys.ChainNo(mobile[c]);
    if (chainno < 0) {
      cout << " Error: could not match chain ID '" << mobile[c] << "'" << endl;
      exit(-1);
    }
    coord = sys.GetChainR(chainno);
    for(a=0; a<sys.NAtoms(chainno); a++) {
      *(coord[a]) += com;
    }
  }

  return(true);
}

/*--------------------------------------------------------------------*/
/* A one-line output of the state of the class */
/*--------------------------------------------------------------------*/
string CRotate::Disp() {
  char      cstring[80];
  string    line;

  sprintf(cstring,"ROTATE (stepsize = %s degrees) (%s %s %s)",
	  Num(stepsize/DEG2RAD).c_str(),Num(phi/DEG2RAD).c_str(),
	  Num(theta/DEG2RAD).c_str(),Num(psi/DEG2RAD).c_str());
  line = cstring;
  return(line);
}

/*--------------------------------------------------------------------*/
/* A two-line output of the state of the class */
/* Requires:  string1 -- label for move type */
/*            string2 -- stepsize for move type */
/*--------------------------------------------------------------------*/
void CRotate::Disp(string& string1, string& string2) {
  string1 = "ROTATE";
  string2 = Num(stepsize/DEG2RAD);
}

/*------------------------------------------------------------------------*/
/* A display routine for the interactions */
/*------------------------------------------------------------------------*/
void CRotate::Display() { 
  cout << "rotate move" << endl;
  cout << "current angles (phi,theta,psi): " 
       << phi << " " << theta << " " << psi << endl;
  cout << "rotation matrix: " << endl << rotmtx;
}
