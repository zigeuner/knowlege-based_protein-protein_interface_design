// A residue-residue type interaction class based on knowledge-based 
// potential evaluations.  This class is intended to identify the pair of 
// residue types and calculate any appropriate interaction potentials

#include "include/config.h"
#include "include/kbtypeint.h"

using std::string;
using std::cout;
using std::ifstream;
using std::ios;
using std::endl;

/*-----------------------------------------------------------------------------*/
/* Initializes the res-res based forcefield interaction analysis */
/* for a specific pair of atom types (eg CA-CA or CB-CB) */
/* Requires:  potentialsdir -- directory location of knowledge-based files */
/*            atomtype1 -- string containing 1st atom type */
/*            atomtype2 -- string containing 2nd atom type */
/*            refindex -- index of reference residue type in AAIntType */
/*            excludelist -- space separated list of residue types to exclude */
/*            verbose -- flag indicating output to be written to screen */
/*-----------------------------------------------------------------------------*/
void CKBTypeInt::Init(const string potentialsdir, const string
		      atomtype1, const string atomtype2, int refindex,
		      const string excludelist, bool verbose) {
  bool     ok;
  int      i,j,k,nexclude,exclude[100],nneg,nover;
  double   factor, max;
  string   basefilename,filename,pair,symbol,field[100],typei,typej;
  ifstream gofrfile;

  atype1 = atomtype1;
  atype2 = atomtype2;
  refresindex = refindex;

  /* set up the base filename */
  basefilename = potentialsdir;
  if (basefilename[basefilename.length() - 1] != '/') {basefilename += "/";}
  basefilename += "gofr_";
  if ((atype1 == "CA")&&(atype2 == "CA")) {
    basefilename += "caca_";
  } else if ((atype1 == "CB")&&(atype2 == "CB")) {
    basefilename += "cbcb_";
  } else {
    cout << "ERROR: not prepared to handle knowledge-based potential " 
	 << "initialization for atom types '" << atype1 << "' and '" 
	 << atype2 << endl;
    exit(-1);
  }

  /* get the list of indices to exclude */
  nexclude = Split(excludelist,field,100);
  for(i=0; i<nexclude; i++) {
    symbol = AAType2Symbol(field[i]);
    exclude[i] = AAIntType2Index(symbol);
    if (exclude[i] < 0) {
      cout << "error: could not interpret 3-letter residue type '" << field[i]
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

  string   oldrestype = "C", newrestype = "S";
  cout << "WARNING: substituting SER for CYS " << atomtype1 << "-" 
       << atomtype2 << " potentials!" << endl;

  /* loop through the upper triangle array indices and initialize */
  for(i=0; i<arraydim; i++) {
    for(j=i; j<arraydim; j++) {

      ok = true;
      for(k=0; k<nexclude; k++) {
	if ((exclude[k] == i)||(exclude[k] == j)) {ok = false; break;}
      }
      if (ok) {
	typei = AAIntType[i];
	typej = AAIntType[j];
	if (typei == oldrestype) {typei = newrestype;}
	if (typej == oldrestype) {typej = newrestype;}
        pair = typei + typej;	
	filename = basefilename + pair + ".dat";

	/* make sure the file exists by trying to open it */
	gofrfile.open(filename.c_str(), ios::in);
	if(!gofrfile.is_open()) {
	  pair = typej + typei;	
	  filename = basefilename + pair + ".dat";
	  gofrfile.close();
	  ifstream file2;
	  file2.open(filename.c_str(), ios::in);
	  if(!file2) {
	    cout << "ERROR: could not open '" << filename << "'" << endl;
	    cout << "ERROR: could not open either permutations of pair " 
		 << pair << endl;
	    exit(-1);
	  }
	  file2.close();
	}	
	gofrfile.close();

	ff[i][j] = new CPairInt(filename,pair,verbose,false);

	/* check the g(r) function */
	max = 4.0;
	(ff[i][j])->Check(max,nneg,nover,false);
	if (nneg > 0) {
	  cout << "ERROR: pair interaction '" << (ff[i][j])->ID()
	       << " with atom types '" << atype1 << "' and '" << atype2 
	       << "' has " << nneg << " values" << endl;
	  exit(-1);
	}
	if (nover > 0) {
	  cout << "WARNING: pair interaction '" << (ff[i][j])->ID()
	       << "' with atom types '" << atype1 << "' and '" << atype2 
	       << "' has " << nover << " values over " << Num(max) << endl;
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

  /* set up the base filename for residue type pair count matrix */
  basefilename = potentialsdir;
  if (basefilename[basefilename.length() - 1] != '/') {basefilename += "/";}
  basefilename += "pairmatrix_";
  if ((atype1 == "CA")&&(atype2 == "CA")) {
    basefilename += "caca";
  } else if ((atype1 == "CB")&&(atype2 == "CB")) {
    basefilename += "cbcb";
  } else {
    cout << "ERROR: not prepared to handle knowledge-based potential " 
	 << "initialization for atom types '" << atype1 << "' and '" 
	 << atype2 << endl;
    exit(-1);
  }

  /* read it into a temporary matrix */
  string   lastcomment;
  string** tempmtx;
  ReadMtxFile(basefilename.c_str(),tempmtx,NIntResTypes-1,lastcomment);

  /* reorder the entries to match internal residue type order */
  int      nfields,index = 0,internalindex[NIntResTypes];
  nfields = Split(lastcomment,field,100);
  for(int i=0; i<NIntResTypes; i++) {
    internalindex[i] = -1;
  }
  for(int i=0; i<nfields; i++) {
    if ((field[i] == "#")||(field[i] == "")) {continue;}
    internalindex[index] = AAIntType2Index(field[i]);
    index += 1;
  }

  /* copy the matrix into permanent storage */
  chimtx.Init(NIntResTypes);
  for(int i=0; i<NIntResTypes-1; i++) {
    for(int j=0; j<NIntResTypes-1; j++) {
      if ((internalindex[i] < 0)||(internalindex[j] < 0)) {continue;}
      chimtx.Set(i,j,ToDbl(tempmtx[internalindex[i]][internalindex[j]]));
    }
  }

  /* normalize the chi matrix */
  chimtx.Symmeterize();
  chimtx.TriangleNormalize();

  /* nullify the unknown residue type entries */
  for(int i=0; i<UnknownResType; i++) {
    chimtx.Set(i,UnknownResType,0.0);
    chimtx.Set(UnknownResType,i,0.0);
  }
  chimtx.Set(UnknownResType,UnknownResType,0.0);

  /* normalize by the random probability of picking a res-res pair */
  factor = (UnknownResType - 1)*(UnknownResType - 1) + (UnknownResType - 1);
  factor /= 2.0;
  for(int i=0; i<UnknownResType; i++) {
    for(int j=0; j<UnknownResType; j++) {
      chimtx.Multiply(i,j,factor);
    }
  }

  /* sum residue-specific entries to get the unknown residue type entries */
  for(int i=0; i<UnknownResType; i++) {
    for(int j=0; j<UnknownResType; j++) {
      chimtx.Add(i,UnknownResType,chimtx.Entry(i,j));
      chimtx.Add(UnknownResType,j,chimtx.Entry(i,j));
      chimtx.Add(UnknownResType,UnknownResType,chimtx.Entry(i,j));
    }
  }

  /* Normalize the unknown residue type entries by number of amino acids */
  /* This is done to give the unknown type average characteristics */
  chimtx.Divide(UnknownResType,UnknownResType,UnknownResType+1);
  chimtx.Divide(UnknownResType,UnknownResType,UnknownResType+1);
  for(int i=0; i<UnknownResType; i++) {
    chimtx.Divide(i,UnknownResType,UnknownResType+1);
    chimtx.Divide(UnknownResType,i,UnknownResType+1);
  }

  /* Check matrix for zeros for all pairs involving the natural 20 or U */
  for(int i=0; i<UnknownResType; i++) {
    for(int j=0; j<UnknownResType; j++) {
      if ((atype1 == "CB")&&(atype2 == "CB")&&
	  ((i == GlycineResType)||(j == GlycineResType))) {continue;}
      if (chimtx.Entry(i,j) <= 0.0) {
	cout << "ERROR: knowledge-based Chi matrix has a zero entry "
	     << "at indices " << i << "," << j << "  (" 
	     << AAIntType[i] << "," << AAIntType[j] << ")" << endl;
	cout << "  This must be corrected before potentials can be calculated." 
	     << endl;
	Display(0);
	exit(-1);
      }
    }
  }
  //  Display(0);
}

/*--------------------------------------------------------------------*/
/* Returns value of a residue-residue potential energy for a single */
/* center-center interaction type on a given pair of residues. */
/* ONLY the g(r) portion is included in this output, see ResResPot */
/* for the residue-residue type potential contribution.  */
/* The potential is relative to the potential if one of the residues */
/* has an alanine or other reference residue type: */
/* dE = E_resZ - E_resALA = kT*ln[(P_AX)/(P_ZX)] */
/* where: P_AX is the probability of residue types ALA and X */
/*        being at their current separation distance */
/* P_AX = g_AX(r) */
/* where:  g_AX(r) is the g(r) for the residue type pair */
/*         at the current separation distance */
/* This returns a potential with units of kcal/mol */
/* NOTE: the SECOND residue is one for which the identity is */
/* referenced to alanine or another reference state */
/* Requires:  residue1 -- pointer to the 1st residue */
/*            resindex1 -- integer residue type */
/*            residue2 -- pointer to the 2nd residue */
/*            resindex2 -- integer residue type */
/*            temperature -- system temperature in Kelvin */
/*            potential -- returned potential */
/*--------------------------------------------------------------------*/
bool CKBTypeInt::Potential(CResidue* res1, const int resindex1, CResidue* res2, 
			   const int resindex2, const double& temperature, 
			   double& potential) {
  int          success;
  double       ref_gofr = 0.0, other_gofr = 0.0, gofr_ratio, kT;
  CVector3     sepvec,coord1,coord2;

  potential = 0.0;
  kT = kcalmole_kb*temperature;

  /* skip evaluation if this interaction pair isn't defined */
  if (ff[resindex1][resindex2] == 0) {return(true);}

  if (!(*res1).AtomRWType(atype1, coord1)) {
    return(true);
    cout << " WARNING: could not find atom type '" << atype1 
	 << "' in passed residue with index '" << resindex1 << "'" << endl;
    (*res1).Display(2);
  }

  if (!(*res2).AtomRWType(atype2, coord2)) {
    return(true);
    cout << " WARNING: could not find atom type '" << atype2 
	 << "' in passed residue with index '" << resindex2 << "'" << endl;
    (*res2).Display(2);
  }

  sepvec = coord1 - coord2;

  success = (ff[resindex1][resindex2])->Potential(sepvec,other_gofr);
  //  cout << success << "    " << other_gofr 
  //       << "  distance: " << sepvec.Length() << endl;

  /* return large interaction potential for separations which are too small */
  if (success == 0) {potential = 100000; return(true);}  

  /* return zero interaction potential for large separations */
  if (success == 2) {potential = 0.0; return(true);}  

  /* finish calculating the potential if other g(r) is interpolable */
  if (success == 1) {
    success = (ff[resindex1][refresindex])->Potential(sepvec,ref_gofr);

    if (success != 1) {
      cout << "ERROR: g(r) interpolation failed on 2nd res-res pair type, " 
	   << "this is unexpected" << endl;
      exit(-1);
    }

    gofr_ratio = ref_gofr/other_gofr;
    potential = kT*log(gofr_ratio);
#ifdef DEBUG
    cout << "DEBUG at " << __FILE__ << ":" << __LINE__ << endl;
    cout << "ratios are: (" << AAIntType[resindex1] << "-" 
	 << AAIntType[refresindex] << ")/(" 
	 << AAIntType[resindex1] << "-" 
	 << AAIntType[resindex2] << ")" << endl;
    cout << "g(r) ratio: " << gofr_ratio << " @ " << sepvec.Length() << endl;
    cout << "kT = " << kT << endl;
    cout << "potential (kcal/mole): " << potential << endl;
#endif
  }

  if (success < 0) {
    cout << "WARNING: g(r) interpolation failed" << endl;
    cout << "atom types are: '" << atype1 << "' and '" << atype2 << "'" << endl;
    cout << coord1 << endl;
    cout << coord2 << endl;
    cout << "separation vector = " << sepvec << endl;
    cout << "separation distance = " << sepvec.Length() << endl;
    return(false);
  }

  return (true);
}

/*--------------------------------------------------------------------*/
/* Returns value of a residue-residue potential energy due to only */
/* identities of the two amino acids. */
/* has an alanine or other reference residue type: */
/* dE = E_resZ - E_resALA = kT*ln[(Chi_AX)/(Chi_ZX)] */
/* where:  Chi_AX is the fraction of residue type - residue type */
/*         pairs that are ALA and residue type X */
/* This returns a potential with units of kcal/mol */
/* NOTE: the SECOND residue is one for which the identity is */
/* referenced to alanine or another reference state */
/* Requires:  residue1 -- pointer to the 1st residue */
/*            resindex1 -- integer residue type */
/*            residue2 -- pointer to the 2nd residue */
/*            resindex2 -- integer residue type */
/*            temperature -- system temperature in Kelvin */
/*            potential -- returned potential */
/*--------------------------------------------------------------------*/
bool CKBTypeInt::ResResPot(CResidue* res1, const int resindex1, 
			   CResidue* res2, const int resindex2, 
			   const double& temperature, double& potential) {
  double       chiratio;
  double       kT = kcalmole_kb*temperature;

  potential = 0.0;

  /* skip evaluation if this interaction pair isn't defined */
  if (ff[resindex1][resindex2] == 0) {
    return(true);
    cout << " WARNING: residue-residue potential is not defined: " 
	 << resindex1 << "--" << resindex2 << endl;
  }

  chiratio = chimtx.Entry(resindex1,refresindex) /
    chimtx.Entry(resindex1,resindex2);
  potential = kT*log(chiratio);

  //    cout << "DEBUG at " << __FILE__ << ":" << __LINE__ << endl;
  //    cout << "ratios are: (" << AAIntType[resindex1] << "-" 
  //	 << AAIntType[refresindex] << ")/(" 
  //	 << AAIntType[resindex1] << "-" 
  //	 << AAIntType[resindex2] << ")" << endl;
  //    cout << "chi ratio: " << chiratio << endl;
  //    cout << "kT = " << kT << endl;
  //    cout << "potential (kcal/mole): " << potential << endl;

  return (true);
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
bool CKBTypeInt::Param(CResidue* res1, const int resindex1, CResidue* res2, 
		       const int resindex2, double& parameter) {
  CVector3     sepvec,coord1,coord2;

  parameter = 0.0;

  /* skip evaluation if this interaction pair isn't defined */
  if (ff[resindex1][resindex2] == 0) {return(true);}

  if (!(*res1).AtomRWType(atype1, coord1)) {
    return(true);
    cout << " WARNING: could not find atom type '" << atype1 
	 << "' in passed residue with index '" << resindex1 << "'" << endl;
    (*res1).Display(2);
  }

  if (!(*res2).AtomRWType(atype2, coord2)) {
    return(true);
    cout << " WARNING: could not find atom type '" << atype2 
	 << "' in passed residue with index '" << resindex2 << "'" << endl;
    (*res2).Display(2);
  }

  sepvec = coord1 - coord2;
  parameter = sepvec.Length();

  return (true);
}

/*--------------------------------------------------------------------*/
/* Returns value of a residue-residue probability for a single */
/* center-center interaction type on a given pair of residues. */
/* ONLY the g(r) portion is included in this output, see ResResPot */
/* for the residue-residue type potential contribution.  */
/* P_AX = Chi_AX*g_AX(r) */
/* where: P_ZX is the probability of residue types Z and X */
/*        being at their current separation distance */
/* where: Chi_ZX is the fraction of residue type - residue type */
/*         pairs that are Z and residue type X */
/* where:  g_ZX(r) is the g(r) for the residue type pair */
/*         at the current separation distance */
/* Requires:  residue1 -- pointer to the 1st residue */
/*            resindex1 -- integer residue type */
/*            residue2 -- pointer to the 2nd residue */
/*            resindex2 -- integer residue type */
/*            prob -- returned array of calculated probabilities */
/*             prob[0] -- g_ZX(r) */
/*             prob[1] -- separation distance */
/* Returns: -1 -- failed */
/*           0 -- under range, returned zero probability */
/*           1 -- ok, within range */
/*           2 -- above range, returned unity probability */
/*           3 -- no interaction defined, returned unity probability */
/*--------------------------------------------------------------------*/
int CKBTypeInt::PairProb(CResidue* res1, const int resindex1, CResidue* res2, 
			 const int resindex2, double*& prob) {
  int          success;
  double       gofr = 0.0;
  CVector3     sepvec,coord1,coord2;

  prob[0] = 1.0;

  /* skip evaluation if this interaction pair isn't defined */
  if (ff[resindex1][resindex2] == 0) {return(3);}

  if (!(*res1).AtomRWType(atype1, coord1)) {
    return(3);
    cout << " WARNING: could not find atom type '" << atype1 
	 << "' in passed residue with index '" << resindex1 << "'" << endl;
    (*res1).Display(2);
  }

  if (!(*res2).AtomRWType(atype2, coord2)) {
    return(3);
    cout << " WARNING: could not find atom type '" << atype2 
	 << "' in passed residue with index '" << resindex2 << "'" << endl;
    (*res2).Display(2);
  }

  sepvec = coord1 - coord2;
  prob[1] = sepvec.Length();

  success = (ff[resindex1][resindex2])->Potential(sepvec,gofr);
  //  cout << success << "    " << gofr 
  //       << "  distance: " << sepvec.Length() << endl;

  /* return zero probabilities for separations which are too small */
  if (success == 0) {
    prob[0] = 0.0;
    return(0);
  }  

  /* return unity probability if it's outside the tabulated max distance */
  if (success == 2) {return(2);}  

  /* finish calculating the probabilities if g(r) is interpolable */
  if (success == 1) {
    //    cout << "DEBUG at " << __FILE__ << ":" << __LINE__ << endl;
    //    cout << "g(r): " << gofr << endl;
    prob[0] = gofr;
  } else {   /* this can't occur, storing error code for future use */
    cout << "WARNING: g(r) interpolation failed" << endl;
    cout << "atom types are: '" << atype1 << "' and '" << atype2 << "'" << endl;
    cout << coord1 << endl;
    cout << coord2 << endl;
    cout << "separation vector = " << sepvec << endl;
    cout << "separation distance = " << sepvec.Length() << endl;
    return(-1);
  }

  return (1);
}

/*--------------------------------------------------------------------*/
/* Same as PairProb, but returns probabilities for a set of residue */
/* types.  The residue type of the... */
/* Returns value of a residue-residue probability for a single */
/* center-center interaction type on a given pair of residues. */
/* ONLY the g(r) portion is included in this output, see ResResPot */
/* for the residue-residue type potential contribution.  */
/* P_AX = Chi_AX*g_AX(r) */
/* where: P_ZX is the probability of residue types Z and X */
/*        being at their current separation distance */
/* where: Chi_ZX is the fraction of residue type - residue type */
/*         pairs that are Z and residue type X */
/* where:  g_ZX(r) is the g(r) for the residue type pair */
/*         at the current separation distance */
/* Requires:  residue1 -- pointer to the 1st residue */
/*            resindex1 -- integer residue type */
/*            residue2 -- pointer to the 2nd residue */
/*            resindex2 -- integer residue type */
/*            prob -- returned array of calculated probabilities */
/*             prob[0] -- separation distance */
/*             prob[N] -- g_ZX(r) for residue type 1...N */
/* Returns: -1 -- failed */
/*           0 -- under range, returned zero probability */
/*           1 -- ok, within range */
/*           2 -- above range, returned unity probability */
/*           3 -- no interaction defined, returned unity probability */
/*--------------------------------------------------------------------*/
/* WRITE THIS LATER IF NEEDED */

/*--------------------------------------------------------------------*/
/* Returns value of a residue-residue probability for a single */
/* residue-residue pair due ONLY to their identities. This probability */
/* Chi_ZX is taken directly from the chi matrix and is the probability */
/* of a given residue-residue pair having a given identity relatively */
/* to random chance.  So, Chi_ZX=1 means that there is no favoring of */
/* residue pairs of type Z-X relative to random chance */
/* Requires:  residue1 -- pointer to the 1st residue */
/*            resindex1 -- integer residue type */
/*            residue2 -- pointer to the 2nd residue */
/*            resindex2 -- integer residue type */
/*--------------------------------------------------------------------*/
double CKBTypeInt::ResResPairProb(CResidue* res1, const int resindex1, 
				  CResidue* res2, const int resindex2) {
  double   prob = 1.0;

  /* skip evaluation if this interaction pair isn't defined */
  if (ff[resindex1][resindex2] == 0) {
    return(prob);
    cout << " WARNING: residue-residue potential is not defined: " 
	 << resindex1 << "--" << resindex2 << endl;
  }

  prob = chimtx.Entry(resindex1,resindex2);
  return (prob);
}

/*--------------------------------------------------------------------*/
/* A display routine to dump forcefield parameters to screen */
/* Requires:  indent -- number of spaces to indent */
/*--------------------------------------------------------------------*/
void CKBTypeInt::Display(int indent) {
  int   i,j;
  char  spacing[100];

  for(i=0;i<indent;i++) {spacing[i] = ' ';}
  spacing[indent] = '\0';

  cout << spacing << "Knowledge-based potentials for pair: " << atype1 << "--"
       << atype2 << endl;
  cout << spacing << "residue-residue potential matrix:" << endl;

  /* Loop through the res type FF array */
  for(i=0; i<arraydim; i++) {
    cout << spacing << " ";
    for(j=0; j<arraydim; j++) {
      if (ff[i][j] == 0) {continue;}
      //      ff[i][j]->Display(indent + 1);
      cout << ff[i][j]->ID() << "  ";
    }
    cout << endl;
  }

  /* Display the chi matrix */
  cout << spacing << "Chi matrix [Nij/(sum_i sum_i>=j Nij)]" << endl;
  cout << spacing << "# ";
  for(i=0; i<NIntResTypes; i++) {
    cout << AAIntType[i] << " ";
  }
  cout << endl;
  chimtx.Display(indent);
  
}
