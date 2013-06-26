// A generic residue-residue type interaction class.  This class is 
// intended to identify the pair of residue types and calculate any
// appropriate interaction potentials.  I've just got it running with
// a generic LJ potential right now that is not residue-type specific.

#include "include/config.h"
#include "include/gentypeint.h"

using std::string;
using std::cout;
using std::ios;
using std::endl;

/*-----------------------------------------------------------------------------*/
/* Initializes the res-res based forcefield interaction analysis */
/* for a specific pair of atom types */
/* Requires:  atomtype1 -- string containing 1st atom type */
/*            atomtype2 -- string containing 2nd atom type */
/*            excludelist -- space separated list of residue types to exclude */
/*-----------------------------------------------------------------------------*/
void CGenTypeInt::Init(const string atomtype1, const string atomtype2, 
		       const string excludelist) {
  bool   ok;
  int    i,j,k,nexclude,exclude[100];
  string filename,pair,symbol,field[100];

  atype1 = atomtype1;
  atype2 = atomtype2;

  /* get the list of indices to exclude */
  nexclude = Split(excludelist,field,100);
  for(i=0; i<nexclude; i++) {
    symbol = AAType2Symbol(field[i]);
    exclude[i] = AAIntType2Index(symbol);
    if (exclude[i] < 0) {
      cout << "ERROR: could not interpret 3-letter residue type '" << field[i]
	   << "', '" << symbol << "', " << exclude[i] 
	   << " as residue to be excluded from interactions" << endl;
      exit(-1);
    }
  }  

  /* allocate the 2D array of pointers */
  arraydim = NIntResTypes;
  ff = new CFFBase** [arraydim];
  for(i=0; i<arraydim; i++) {
    ff[i] = new CFFBase* [arraydim];
  }    

  /* loop through the upper triangle array indices and initialize */
  for(i=0; i<arraydim; i++) {
    for(j=i; j<arraydim; j++) {

      ok = true;
      for(k=0; k<nexclude; k++) {
	if ((exclude[k] == i)||(exclude[k] == j)) {ok = false; break;}
      }
      if (ok) {
	//	ff[i][j] = new CLJInt(5.0,100.0,pair);   /* original LJ pots */

	/* sigma's 2x the radi used by Huang, Love, Mayo (2005) JCP v26 p1222 */
	if ((atomtype1 == "CB")&&(atomtype2 == "CB")) {
	  ff[i][j] = new CLJInt(4.3,100.0,pair); 
	  //	  ff[i][j] = new CLJInt(5.9,100.0,pair);  /* Val-Val */
	  //	  ff[i][j] = new CLJInt(6.5,100.0,pair);  /* Leu-Leu */
	} else {
	  ff[i][j] = new CLJInt(4.7,100.0,pair); 
	  //	  ff[i][j] = new CLJInt(6.1,100.0,pair);  /* Val-Val */
	  //	  ff[i][j] = new CLJInt(8.3,100.0,pair);  /* Leu-Leu */
	}

      } else {
	ff[i][j] = 0; 
      }
    }
  }

  /* Reflect the res-res evaluation classes into the lower triangle */
  for(i=0; i<arraydim; i++) {
    for(j=0; j<i; j++) {
      ff[i][j] = ff[j][i];
    }
  }

}

/*--------------------------------------------------------------------*/
/* Returns value of a residue-residue potential energy for a given */
/* pair of residues */
/* Returns false if there was an error during evaluation */
/* Requires:  residue1 -- pointer to the 1st residue */
/*            resindex1 -- integer residue type */
/*            residue2 -- pointer to the 2nd residue */
/*            resindex2 -- integer residue type */
/*            temperature -- system temperature in Kelvin (UNUSED) */
/*            potential -- returned potential */
/*--------------------------------------------------------------------*/
bool CGenTypeInt::Potential(CResidue* res1, const int resindex1, CResidue* res2, 
			    const int resindex2, const double& temperature, 
			    double& potential) {
  CVector3     sepvec,coord1,coord2;

  potential = 0.0;
  /* skip evaluation if this interaction pair isn't defined */
  if (ff[resindex1][resindex2] == 0) {return(true);}

  if (!(*res1).AtomRWType(atype1, coord1)) {
    return(true);   /* returns zero potential, no interaction */
    cout << " WARNING: could not find atom type '" << atype1 
	 << "' in passed residue with index '" << resindex1 << "'" << endl;
    (*res1).Display(2);
  }

  if (!(*res2).AtomRWType(atype2, coord2)) {
    return(true);   /* returns zero potential, no interaction */
    cout << " WARNING: could not find atom type '" << atype2 
	 << "' in passed residue with index '" << resindex2 << "'" << endl;
    (*res2).Display(2);
  }

  sepvec = coord1 - coord2;
  if ((ff[resindex1][resindex2])->Potential(sepvec,potential) < 0) {
    cout << "WARNING: potential evaluation failed" << endl;
    cout << coord1 << endl;
    cout << coord2 << endl;
    cout << "separation vector = " << sepvec << endl;
    return(false);
  }

  return (true);
}

/*--------------------------------------------------------------------*/
/* Returns value of a residue-residue potential energy due to only */
/* identities of the two amino acids. */
/* NOTE: this is zero for this type of forcefield */
/* Requires:  residue1 -- pointer to the 1st residue */
/*            resindex1 -- integer residue type */
/*            residue2 -- pointer to the 2nd residue */
/*            resindex2 -- integer residue type */
/*            temperature -- system temperature in Kelvin */
/*            potential -- returned potential */
/*--------------------------------------------------------------------*/
bool CGenTypeInt::ResResPot(CResidue* res1, const int resindex1, 
			    CResidue* res2, const int resindex2, 
			    const double& temperature, double& potential) {
  potential = 0.0;
  return(true);
}

/*--------------------------------------------------------------------*/
/* Returns value of the position-dependant parameter used to calculate */
/* the interaction.  This is the distance in most cases */
/* Requires:  residue1 -- pointer to the 1st residue */
/*            resindex1 -- integer residue type */
/*            residue2 -- pointer to the 2nd residue */
/*            resindex2 -- integer residue type */
/*            parameter -- returned parameter */
/*--------------------------------------------------------------------*/
bool CGenTypeInt::Param(CResidue* res1, const int resindex1, CResidue* res2, 
			const int resindex2, double& parameter) {
  CVector3     sepvec,coord1,coord2;

  parameter = 0.0;

  /* skip evaluation if this interaction pair isn't defined */
  if (ff[resindex1][resindex2] == 0) {return(true);}

  if (!(*res1).AtomRWType(atype1, coord1)) {
    return(true);   /* returns zero potential, no interaction */
    cout << " WARNING: could not find atom type '" << atype1 
	 << "' in passed residue with index '" << resindex1 << "'" << endl;
    (*res1).Display(2);
  }

  if (!(*res2).AtomRWType(atype2, coord2)) {
    return(true);   /* returns zero potential, no interaction */
    cout << " WARNING: could not find atom type '" << atype2 
	 << "' in passed residue with index '" << resindex2 << "'" << endl;
    (*res2).Display(2);
  }

  sepvec = coord1 - coord2;
  parameter = sepvec.Length();

  return (true);
}

/*--------------------------------------------------------------------*/
/* Returns value of a residue-residue potential energy for a given */
/* pair of residues */
/* Returns false if there was an error during evaluation */
/* Requires:  residue1 -- pointer to the 1st residue */
/*            resindex1 -- integer residue type */
/*            residue2 -- pointer to the 2nd residue */
/*            resindex2 -- integer residue type */
/*            prob -- returned array of calculated probabilities */
/*             prob[0] -- P_ZX */
/*             prob[1] -- g_ZX(r) */
/*             prob[2] -- Chi_ZX */
/*--------------------------------------------------------------------*/
int CGenTypeInt::PairProb(CResidue* res1, const int resindex1, CResidue* res2, 
			  const int resindex2, double*& prob) {
  cout << "ERROR: cannot calculate pair probabilities using " 
       << "general forcefield" << endl;
  exit(-1);
}

/*--------------------------------------------------------------------*/
/* Returns value of a residue-residue probability for a single */
/* residue-residue pair due ONLY to their identities. It's always unity */
/* for this forcefield type */
/* Requires:  residue1 -- pointer to the 1st residue */
/*            resindex1 -- integer residue type */
/*            residue2 -- pointer to the 2nd residue */
/*            resindex2 -- integer residue type */
/*--------------------------------------------------------------------*/
double CGenTypeInt::ResResPairProb(CResidue* res1, const int resindex1, 
				   CResidue* res2, const int resindex2) {
  return (1.0);
}

/*--------------------------------------------------------------------*/
/* A display routine to dump forcefield parameters to screen */
/* Requires:  indent -- number of spaces to indent */
/*--------------------------------------------------------------------*/
void CGenTypeInt::Display(int indent) {
  int   i,j;
  char  spacing[100];

  for(i=0;i<indent;i++) {spacing[i] = ' ';}
  spacing[indent] = '\0';

  /* Loop through the res type FF array */
  for(i=0; i<arraydim; i++) {
    cout << spacing << "  ";
    for(j=0; j<arraydim; j++) {
      ff[i][j]->Display(indent + 1);
      cout << "  ";
    }
    cout << endl;
  }
  
}
