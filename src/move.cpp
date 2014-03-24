// Functions for the move class. 

#include "include/config.h"
#include "include/move.h"

using std::string;
using std::cout;
using std::ios;
using std::endl;

/*--------------------------------------------------------------------*/
/* Constructor, initializes the res-res interaction array */
/* Requires:  movelist -- array of strings specifying move types */
/*            listsize -- number of move types */
/*            chains -- comma-separated list of chain ID's to move */
/*            verbose -- flag indicating output to be written to screen */
/*--------------------------------------------------------------------*/
CMove::CMove(string* movelist, int listsize, string chains, bool verbose) {
  Init(movelist,listsize,chains,verbose);
}

/*--------------------------------------------------------------------*/
/* Initializes the array of move types */
/* Requires:  movelist -- array of strings specifying move types */
/*            listsize -- number of move types */
/*            chains -- comma-separated list of chain ID's to move */
/*            verbose -- flag indicating output to be written to screen */
/*--------------------------------------------------------------------*/
void CMove::Init(string* movelist, int listsize, string chains, bool verbose) {
  string     temp[100];

  /* Set the mobile chain(s) definition */
  nmobile = Split(chains,',',temp,100);
  mobile = new string[nmobile];
  for(int i=0; i<nmobile; i++) {mobile[i] = temp[i];}

  /* Allocate the 2D array of pointers */
  nmoves = listsize;
  movetype = new CMoveBase* [nmoves];

  if (verbose) {
    cout << " initializing " << nmoves << " move type(s) on chain(s): '" 
	 << chains << "'" << endl;
  }

  /* loop through the upper triangle array indices and initialize */
  for(int i=0; i<nmoves; i++) {

    /* Initialize each move type, picking the class according to specs */
    if (movelist[i] == "ROTATE") {
      if (verbose) {cout << " initializing ROTATE type" << endl;}
      movetype[i] = new CRotate(chains);
    } else if (movelist[i] == "TRANSLATE") {
      if (verbose) {cout << " initializing TRANSLATE type" << endl;}
      movetype[i] = new CTranslate(chains);
    } else {
      cout << "ERROR: Could not interpret move type '" << movelist[i] <<
	"'" << endl;
      exit(-1);
    }
  }

}

/*--------------------------------------------------------------------*/
/* Perform the move(s), getting the temperature from the system class */
/* Requires:  sys -- system */
/*            ff -- forcefield */
/*            verbose -- flag that controls verbosity */
/*--------------------------------------------------------------------*/
void CMove::Move(CSystem& sys, CForcefield& ff, const bool verbose) {
  Move(sys,ff,sys.T(),verbose);
}

/*--------------------------------------------------------------------*/
/* Perform the move(s) */
/* Requires:  sys -- system */
/*            ff -- forcefield */
/*            T -- temperature */
/*            verbose -- flag that controls verbosity */
/*--------------------------------------------------------------------*/
void CMove::Move(CSystem& sys, CForcefield& ff, const double T, 
		 const bool verbose) {
  bool       accept, overflow;
  int        moveno,a,c,chainno;
  double     initpot, finalpot, dpot, exponent, boltz;
  CVector3   coord[nmobile][sys.MaxNAtoms()];   // got to fix this (March 2014)
  //  CVector3   coord[1000][10000];
  CVector3** chaincoord = 0;

  /* copy relevant chain coordinates as backup */
  for(c=0; c<nmobile; c++) {
    chainno = sys.ChainNo(mobile[c]);
    chaincoord = sys.GetChainR(chainno);
    for(a=0; a<sys.NAtoms(chainno); a++) {
      coord[c][a] = *(chaincoord[a]);
    }
  }

  /* pick a move type at random */
  moveno = int(floor((nmoves - 1.0e-6)*sys.GetRan()));
  if (verbose) {
    cout << "Attempting move " << movetype[moveno]->ID() << endl;
  }

  /* evaluate the sub-system potential energy */
  initpot = ff.Potential(sys);
  //  initpot = forcefield.CurrentPotential();  // would be faster!

  /* make the move */
  //  sys.WritePDB("premove.pdb");
  movetype[moveno]->Move(sys,ff);
  //  sys.WritePDB("postmove.pdb");

  /* evaluate the sub-system potential energy */
  finalpot = ff.Potential(sys);
  dpot = finalpot - initpot;
  if (verbose) {
    cout << " dV = V_final - V_initial = " << dpot << " = " 
	 << finalpot << " - " << initpot << endl;
  }

  /* calculate Boltzmann factor and accept/reject move */
  accept = false;
  overflow = false;
  if (dpot < 0.0) {
    accept = true;
  } else {
    exponent = -1.0*dpot/(kcalmole_kb*T);
    if ((exponent < MINEXP)||(exponent > MAXEXP)) {
      if (verbose) {
	cout << "WARNING: will exceed range if number is exponentiated '" <<
	  exponent << "'" << endl;
	cout << " movetype " << movetype[moveno]->Disp() << endl;
	cout << " reducing stepsize by 2x and reseting coordinates" << endl;
      }
      movetype[moveno]->ResizeStep(0.50);
      movetype[moveno]->EnforceMinStep(0.001);
      overflow = true;
    } else {
      boltz = exp(exponent);
      if (boltz > sys.GetRan()) {accept = true;}
    }
  }

  if (!overflow) {
    /* modify the stepsize */
    if (accept) {
      movetype[moveno]->ResizeStep(1.05);
      movetype[moveno]->EnforceMaxStep(5.0);
    } else {
      movetype[moveno]->ResizeStep(0.95);
      movetype[moveno]->EnforceMinStep(0.01);
    }

    if (verbose) {
      if (accept) {
	cout << " accepted movetype " << movetype[moveno]->Disp() << endl;
      } else {
	cout << " rejected movetype " << movetype[moveno]->Disp() << endl;
      }
    }
  }

  /* restore relevant chain coordinates if rejected */
  if (!accept) {
    for(c=0; c<nmobile; c++) {
      chainno = sys.ChainNo(mobile[c]);
      chaincoord = sys.GetChainR(chainno);
      for(a=0; a<sys.NAtoms(chainno); a++) {
	*(chaincoord[a]) = coord[c][a];
      }
    }
    /* reset stored system potential energy in forcefield class */
    ff.SetCurrentPotential(initpot);
  }
  
}

/*--------------------------------------------------------------------*/
/* A one-line output of the state of the class */
/*--------------------------------------------------------------------*/
string CMove::Disp() {
  char      cstring[80];
  string    line;

  sprintf(cstring,"Moves state: ");
  for(int i=0; i<nmoves; i++) {
    sprintf(cstring,"%s %s",cstring,(movetype[i]->Disp()).c_str());
  }
  return(line = cstring);
}

/*--------------------------------------------------------------------*/
/* A two-line output of the state of the class */
/* Requires:  string1 -- labels for move types */
/*            string2 -- stepsizes for move types */
/*--------------------------------------------------------------------*/
void CMove::Disp(string& string1, string& string2) {
  string    str1 = "", str2 = "";
  string1 = "";
  string2 = "";
  for(int i=0; i<nmoves; i++) {
    movetype[i]->Disp(str1,str2);
    string1 += str1;
    string2 += str2;
    if (i < nmoves - 1) {
      string1 += "  ";
      string2 += "  ";
    }
  }
}

/*--------------------------------------------------------------------*/
/* A display routine to dump move parameters to screen */
/* Requires:  indent -- number of spaces to indent */
/*--------------------------------------------------------------------*/
void CMove::Display(int indent) {
  int   i;
  char  spacing[100];

  for(i=0;i<indent;i++) {spacing[i] = ' ';}
  spacing[indent] = '\0';

  cout << spacing << Disp() << endl;
  
}

