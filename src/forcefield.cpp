// Functions for the forcefield class.  These routines calculate energies
// and/or forces for the system.

#include "include/config.h"
#include "include/forcefield.h"

using std::string;
using std::cout;
using std::ios;
using std::endl;

/*--------------------------------------------------------------------*/
/* Constructor, initializes the res-res interaction array */
/* Requires:  forcefieldtype -- name defining type of forcefield */
/*            potentialsdir -- string containing directory name */
/*            refresidue -- reference residue string (3-letter) */
/*            verbose -- flag indicating output to be written to screen */
/*--------------------------------------------------------------------*/
CForcefield::CForcefield(const string forcefieldtype, const string potentialsdir,
			 const string refresidue, bool verbose) {
  fftype = forcefieldtype;
  Init(potentialsdir,refresidue,verbose);
}

/*--------------------------------------------------------------------*/
/* Sets the subset portion of the forcefield by copying from another */
/* subset class */
/* Requires:  oldsubset -- subset class to copy */
/*--------------------------------------------------------------------*/
void CForcefield::SetSubset(const CSubset& oldsubset) {
  intsubset = true;
  subset = oldsubset;
}

/*--------------------------------------------------------------------*/
/* Sets the subset portion of the forcefield by reading from a file */
/* Requires:  residuesfile -- name of file containing residues specs */
/*            verbose -- flag indicating output to be written to screen */
/*--------------------------------------------------------------------*/
void CForcefield::SetSubset(const string residuesfile, bool verbose) {
  /* Read definition of interacting residues if necessary */
  if (residuesfile != "NULL") {
    intsubset = true;
    subset.Init(residuesfile,verbose);
  }
}

/*--------------------------------------------------------------------*/
/* Initializes the res-res based forcefield interaction analysis */
/* Two "sides" or sets of chains are defined between which interactions */
/* will be evaluated */
/* Requires:  potentialsdir -- string containing directory name */
/*            refresidue -- reference residue string (3-letter) */
/*            verbose -- flag indicating output to be written to screen */
/*--------------------------------------------------------------------*/
void CForcefield::Init(const string potentialsdir, const string refresidue,
		       bool verbose) {
  int    i;
  string filename,pair;

  /* set some defaults, with -O3 they aren't set by default constructor */
  sysnrg = 0.0;
  subsetseq = false;
  evalchipot = false;
  chidistcutoff = 12.0; 
  intsubset = false; 
  useunk = false;

  /* determine which residue index to use as reference */
  refresindex = AAType2IntIndex(refresidue);
  if (refresindex == -1) {
    cout << "ERROR: unable to convert 3-letter code '" 
	 << refresidue << "' to interaction residue index" << endl;
    exit(-1);    
  }

  /* Allocate the array of pointers */
  npairtypes = 2;
  if (fftype.find("CA_ONLY",0) != string::npos) {npairtypes = 1;}
  intset = new CMtxIntBase* [npairtypes];

  /* loop through the residue pair interaction types and initialize */
  for(i=0; i<npairtypes; i++) {

    if (fftype.find("SIMPLE",0) != string::npos) {
      if (i == 0) {
        intset[i] = new CGenTypeInt("CA","CA","HET"); 
      } else if (i == 1) {
        intset[i] = new CGenTypeInt("CB","CB","GLY HET"); 
      } else {
        cout << "ERROR: unable to interpret numeric pairtype '" 
      	   << i << "'" << endl;
        exit(-1);
      } 
    } else if (fftype.find("KB",0) != string::npos) {
      if (i == 0) {
        intset[i] = new CKBTypeInt(potentialsdir,"CA","CA",
				   refresindex,"HET",verbose);
      } else if (i == 1) {
        intset[i] = new CKBTypeInt(potentialsdir,"CB","CB",
				   refresindex,"GLY HET",verbose); 
      } else {
        cout << "ERROR: unable to interpret numeric pairtype '" 
      	   << i << "'" << endl;
        exit(-1);
      } 
    } else {
      cout << "ERROR: unable to interpret forcefield type '" << fftype 
	   << "'" << endl;
    }
  }

  /* initialize a generic 2.0 Angstrom radius hard sphere potential */
  hspotential.Init(100.0,4.0,"HS");

}

/*--------------------------------------------------------------------*/
/* Returns value of repulsive potential energy */
/* Requires:  sys -- system */
/*--------------------------------------------------------------------*/
double CForcefield::RepulsivePot(CSystem& sys) {
  int      a1, a2, c1, c2, nchains;
  int      natoms1, natoms2;
  double   potential = 0.0, pot;
  string   chainid1, chainid2;
  bool*    interaction;
  CVector3 sepvec;
  CAtom    **atom1,**atom2;

  if (!intsubset) {
    cout << "not prepared for forcefield evaluations w/o subset" << endl;
    exit(-1);
  }

  nchains = subset.NChainTypes();
  for(c1=0; c1<nchains; c1++) {
    chainid1 = subset.ChainID(c1);
    natoms1 = sys.GetRepulsiveAtoms(chainid1,atom1);
    //    cout << "pointer in RepulPot: " << atom1 << endl;
    //    for(a1=0; a1<1; a1++) {
    //      cout << a1 << "  " << atom1[a1] << endl;
    //      (*(atom1[a1])).Display(2);
    //      cout << a1 << "  " << (*(atom1[a1])).R() << endl;
    //    }
    //    exit(0);

    interaction = subset.GetIntChains(chainid1);
    for(c2=c1; c2<nchains; c2++) {
      /* check to see if this chain-chain pair is interacting */
      if (!interaction[c2]) {continue;}

      /* get the interaction residues in chain 2 */
      chainid2 = subset.ChainID(c2);
      natoms2 = sys.GetRepulsiveAtoms(chainid2,atom2);

      for(a1=0; a1<natoms1; a1++) {
	for(a2=0; a2<natoms2; a2++) {
	  sepvec = (atom1[a1])->R() - (atom2[a2])->R();
	  hspotential.Potential(sepvec,pot);
	  //	  cout << c1 << ":" << a1 << "-" 
	  //	       << c2 << ":" << a2 << "  " << pot << endl;
	  potential += pot;
	}
      }
    }
  }

  return(potential);

}

/*--------------------------------------------------------------------*/
/* Returns value of total potential energy */
/* Requires:  sys -- system */
/*--------------------------------------------------------------------*/
double CForcefield::Potential(CSystem& sys) {
  int          i, r1, r2, c1, c2, nchains;
  int          nres1, nres2, resindex1, resindex2;
  int          maxnres = subset.MaxNResidues();
  bool         ligand1, ligand2, flip, success, evaldist;
  double       totalpotential = 0.0, resrespot, chipotential, potential;
  string       chainid1, chainid2;
  CVector3     coord1, coord2, sepvec;
  bool*        interaction;
  CResidue*    res1[maxnres];
  CResidue*    res2[maxnres];

  if (!intsubset) {
    cout << "not prepared for forcefield evaluations w/o subset" << endl;
    exit(-1);
  }
  
  nchains = subset.NChainTypes();
  for(c1=0; c1<nchains; c1++) {
    chainid1 = subset.ChainID(c1);
    nres1 = sys.GetResidues(chainid1,subset.GetResidues(chainid1),
			    res1,subset.NResidues(chainid1));

    interaction = subset.GetIntChains(chainid1);
    for(c2=c1; c2<nchains; c2++) {
      /* check to see if this chain-chain pair is interacting */
      if (!interaction[c2]) {continue;}

      /* get the interaction residues in chain 2 */
      chainid2 = subset.ChainID(c2);
      nres2 = sys.GetResidues(chainid2,subset.GetResidues(chainid2),
			      res2,subset.NResidues(chainid2));

      //      cout << " evaluating potential for chain '" << chainid1
      //	   << "' with chain '" << chainid2 << "'" << endl;

      /* Change residue types to ligand if chain(s) is flagged as unknown */
      ligand1 = false; ligand2 = false; 
      if (subset.ChkUnknown(chainid1)) {
	ligand1 = true; resindex1 = UnknownResType;
      }
      if (subset.ChkUnknown(chainid2)) {
	ligand2 = true; resindex2 = UnknownResType;
      }

      /* if just one residue type is unknown, want it to be 2nd */
      /* the unknown or ligand residue is type-switched for reference */
      flip = false;
      if ((ligand1)&&(!ligand2)) {flip = true;}

      /* keep track of the chi-derived res-res potential by ligand residue */
      chipotential = 0.0;

      /* loop over the residue pairs and calculate potentials */
      for(r1=0; r1<nres1; r1++) {
	if ((!useunk)||(!ligand1)) {resindex1 = res1[r1]->ResIndex();}
	if ((ligand1)&&(subsetseq)) {resindex1 = subset.AltResIndex(c1,r1);}

	for(r2=0; r2<nres2; r2++) {
	  if ((!useunk)||(!ligand2)) {resindex2 = res2[r2]->ResIndex();}
	  if ((ligand2)&&(subsetseq)) {resindex2 = subset.AltResIndex(c2,r2);}

	  //	  cout << chainid1 << ":" << r1 << " -- " << chainid2 
	  //	       << ":" << r2 << endl;
	  /* loop over residue pair interaction types */
	  resrespot = 0.0;
	  for(i=0; i<npairtypes; i++) {
	    //	    cout << chainid1 << resindex1 << "  " 
	    //		 << chainid2 << resindex2 << endl;
	    if (flip) {
	      success = (intset[i])->Potential(res2[r2],resindex2,res1[r1],
					       resindex1,sys.T(),potential);
	    } else {
	      success = (intset[i])->Potential(res1[r1],resindex1,res2[r2],
					       resindex2,sys.T(),potential);
	    }
	    //	    cout << "returned potential: " << potential << endl;
	    if (!success) {
	      cout << "WARNING: potential evaluation failed" << endl;
	      cout << chainid1 << ":" << r1 << " -- " << chainid2 
		   << ":" << r2 << endl;
	      potential = 0;
	      exit(-1);
	    }
	    //	  cout << chainid1 << ":" << r1 << " -- " << chainid2 
	    //	       << ":" << r2 << "  potential= " << potential << endl;
	    resrespot += potential;
	  }
	  totalpotential += resrespot;  /* add basic pair potential */

	  /* evaluate restype-restype potentials if necessary */
	  if (evalchipot) {
	    evaldist = true;
	    /* get Calpha-Calpha distance for restype-restype potential */
	    if (!((*res1[r1]).AtomRWType("CA", coord1))) {
	      //cout << "ERROR: could not find Calpha in residue" << endl;
	      //(*res1[r1]).Display(1);
	      //exit(-1);
	      evaldist = false;
	    }
	    if (!((*res2[r2]).AtomRWType("CA", coord2))) {
	      //cout << "ERROR: could not find Calpha in residue" << endl;
	      //(*res2[r2]).Display(1);
	      //exit(-1);
	      evaldist = false;
	    }
	    if (evaldist) {
	      sepvec = coord1 - coord2;
	      if (sepvec.Length() < chidistcutoff) {
		if (flip) {
		  success = (intset[0])->ResResPot(res2[r2],resindex2,res1[r1],
						   resindex1,sys.T(),potential);
		  resrespot = potential/nres2;
		} else {
		  success = (intset[0])->ResResPot(res1[r1],resindex1,res2[r2],
						   resindex2,sys.T(),potential);
		  resrespot = potential/nres1;
		}
		chipotential += resrespot;
	      }
	    }
	  }

	}  /* residue 2 loop */
      }  /* residue 1 loop */

      totalpotential += chipotential;  /* add chi-derived potential */

    }  /* chain 2 loop */
  }  /* chain 1 loop */

  /* handle the repulsive-center interactions */
  potential = RepulsivePot(sys);
  totalpotential += potential;  /* add additional repulsive potential */

  /* save the current potential and return */
  //  cout << "potential = " << totalpotential << endl;
  sysnrg = totalpotential;
  return (totalpotential);
}

/*--------------------------------------------------------------------*/
/* Returns value of total potential energy */
/* This version also fills a matrix of res-res interaction */
/* The passed matrix class MUST ALREADY BE SIZED to the number */
/* of interacting residues in the forcefield's subset definition */
/* Requires:  sys -- system */
/*            potmtx -- matrix of res-res potentials */
/*            commentline -- string of residue identifiers for output */
/*--------------------------------------------------------------------*/
double CForcefield::Potential(CSystem& sys, CMatrix& potmtx, 
			      string& commentline) {
  int          i, r1, r2, c1, c2, nchains, indx,jndx;
  int          nres1, nres2, resindex1, resindex2;
  int          maxnres = subset.MaxNResidues();
  bool         ligand1, ligand2, flip, success, evaldist;
  double       totalpotential = 0.0, resrespot, chipotential, potential;
  string       chainid1, chainid2;
  CVector3     coord1, coord2, sepvec;
  bool*        interaction;
  CResidue*    res1[maxnres];
  CResidue*    res2[maxnres];

  if (!intsubset) {
    cout << "not prepared for forcefield evaluations w/o subset" << endl;
    exit(-1);
  }

  if (subset.TotalResidues() != potmtx.Size()) {
    cout << "ERROR: passed matrix class must have size corresponding to "
	 << "number of interface residues" << endl;
    exit(-1);
  }

  nchains = subset.NChainTypes();
  for(c1=0; c1<nchains; c1++) {
    chainid1 = subset.ChainID(c1);
    nres1 = sys.GetResidues(chainid1,subset.GetResidues(chainid1),
			    res1,subset.NResidues(chainid1));

    interaction = subset.GetIntChains(chainid1);
    for(c2=c1; c2<nchains; c2++) {
      /* check to see if this chain-chain pair is interacting */
      if (!interaction[c2]) {continue;}

      /* get the interaction residues in chain 2 */
      chainid2 = subset.ChainID(c2);
      nres2 = sys.GetResidues(chainid2,subset.GetResidues(chainid2),
			      res2,subset.NResidues(chainid2));

      //      cout << " evaluating potential for chain '" << chainid1
      //	   << "' with chain '" << chainid2 << "'" << endl;

      /* Change residue types to ligand if chain(s) is flagged as unknown */
      ligand1 = false; ligand2 = false; 
      if (subset.ChkUnknown(chainid1)) {
	ligand1 = true; resindex1 = UnknownResType;
      }
      if (subset.ChkUnknown(chainid2)) {
	ligand2 = true; resindex2 = UnknownResType;
      }

      /* if just one residue type is unknown, want it to be 2nd */
      /* the unknown or ligand residue is type-switched for reference */
      flip = false;
      if ((ligand1)&&(!ligand2)) {flip = true;}

      /* keep track of the chi-derived res-res potential by ligand residue */
      chipotential = 0.0;

      /* loop over the residue pairs and calculate potentials */
      for(r1=0; r1<nres1; r1++) {
	if ((!useunk)||(!ligand1)) {resindex1 = res1[r1]->ResIndex();}
	if ((ligand1)&&(subsetseq)) {resindex1 = subset.AltResIndex(c1,r1);}
	indx = subset.ID2Num(chainid1,r1);

	if (resindex1 < 0) {
	  std::cout << "DEBUG at " << __FILE__ << ":" << __LINE__ << std::endl;
	  cout << subsetseq << endl;
	  cout << "ERROR: residue type is negative" << endl;
	  cout << "chain ID: " << chainid1 << endl;
	  cout << "residue number: " << res1[r1]->ResNo() << endl;
	  subset.Display(0);
	  exit(0);
	}

	for(r2=0; r2<nres2; r2++) {
	  if ((!useunk)||(!ligand2)) {resindex2 = res2[r2]->ResIndex();}
	  if ((ligand2)&&(subsetseq)) {resindex2 = subset.AltResIndex(c2,r2);}
	  jndx = subset.ID2Num(chainid2,r2);

	  if (resindex2 < 0) {
	    cout << "ERROR: residue type is negative" << endl;
	    cout << "chain ID: " << chainid2 << endl;
	    cout << "residue number: " << res2[r2]->ResNo() << endl;
	    exit(0);
	  }

	  //	  cout << chainid1 << ":" << r1 << " -- " << chainid2 
	  //	       << ":" << r2 << endl;
	  /* loop over residue pair interaction types */
	  resrespot = 0.0;
	  for(i=0; i<npairtypes; i++) {
	    //	    cout << chainid1 << resindex1 << "  " 
	    //		 << chainid2 << resindex2 << endl;
	    if (flip) {
	      success = (intset[i])->Potential(res2[r2],resindex2,res1[r1],
					       resindex1,sys.T(),potential);
	    } else {
	      success = (intset[i])->Potential(res1[r1],resindex1,res2[r2],
					       resindex2,sys.T(),potential);
	    }
	    //	    cout << "returned potential: " << potential << endl;
	    if (!success) {
	      cout << "WARNING: potential evaluation failed" << endl;
	      cout << chainid1 << ":" << r1 << " -- " << chainid2 
		   << ":" << r2 << endl;
	      potential = 0;
	      exit(-1);
	    }
	    //	  cout << chainid1 << ":" << r1 << " -- " << chainid2 
	    //	       << ":" << r2 << "  potential= " << potential << endl;
	    resrespot += potential;
	  }
	  totalpotential += resrespot;  /* add basic pair potential */
	  potmtx.Set(indx,jndx,resrespot); /* set basic pair potential */

	  /* evaluate restype-restype potentials if necessary */
	  if (evalchipot) {
	    evaldist = true;
	    /* get Calpha-Calpha distance for restype-restype potential */
	    if (!((*res1[r1]).AtomRWType("CA", coord1))) {
	      //cout << "ERROR: could not find Calpha in residue" << endl;
	      //(*res1[r1]).Display(1);
	      //exit(-1);
	      evaldist = false;
	    }
	    if (!((*res2[r2]).AtomRWType("CA", coord2))) {
	      //cout << "ERROR: could not find Calpha in residue" << endl;
	      //(*res2[r2]).Display(1);
	      //exit(-1);
	      evaldist = false;
	    }
	    if (evaldist) {
	      sepvec = coord1 - coord2;
	      if (sepvec.Length() < chidistcutoff) {
		if (flip) {
		  success = (intset[0])->ResResPot(res2[r2],resindex2,res1[r1],
						   resindex1,sys.T(),potential);
		  resrespot = potential/nres2;
		} else {
		  success = (intset[0])->ResResPot(res1[r1],resindex1,res2[r2],
						   resindex2,sys.T(),potential);
		  resrespot = potential/nres1;
		}
		potmtx.Add(indx,jndx,resrespot);  /* add chi-derived potential */
		chipotential += resrespot;
	      }
	    }
	  }

	}  /* residue 2 loop */
      }  /* residue 1 loop */

      totalpotential += chipotential;  /* add chi-derived potential */

    }
  }

  /* handle the repulsive-center interactions */
  potential = RepulsivePot(sys);
  cout << "Repulsive potential= " << Num(potential) << endl;
  totalpotential += potential;  /* add additional repulsive potential */

  /* make the string of residue definitions as a comment line */
  int tempno;
  string tempid,tempresno;
  commentline = "# ";
  for(int i=0; i<subset.TotalResidues(); i++) {
    tempno = subset.Num2ID(i,tempid,tempresno);
    commentline += tempid + ":" + sys.ResSymbol(tempid,tempresno);
    commentline += tempresno +  " ";
  }  

  sysnrg = totalpotential;
  return (totalpotential);
}

/*--------------------------------------------------------------------*/
/* Evaluates the probability of each amino acid type being at a each */
/* of the unknown residue positions */
/* Requires:  sys -- system */
/*            output -- array of strings containing output info */
/*            maxlines -- size of array of strings */
/*--------------------------------------------------------------------*/
int CForcefield::EvalAlphabet(CSystem& sys, string*& output, 
			      const int maxlines) {
  int          i, restype, r1, r2, c1, c2, nchains, nneighbors;
  int          nres1, nres2, resindex1, resindex2, nchains2, nlines, success;
  bool         neighbor;
  double       resprobsum,resgofrsum,distsum;
  double       resprob[NIntResTypes],resgofr[NIntResTypes];
  double*      prob;
  char         cstring[500];
  string       chainid1, chainid2, diststring;
  bool*        interaction;
  CResidue*    res1[subset.MaxNResidues()];
  CResidue*    res2[subset.MaxNResidues()];

  if (!intsubset) {
    cout << "subset must be defined before alphabets can be evaluated" << endl;
    exit(-1);
  }
  if (subset.NUnknown() == 0) {return(0);}

  prob = new double[4];

  nchains = subset.NUnknown();
  nchains2 = subset.NChainTypes();
  for(c1=0; c1<nchains; c1++) {
    chainid1 = subset.UnknownChainID(c1);
    nres1 = sys.GetResidues(chainid1,subset.GetResidues(chainid1),
			    res1,subset.NResidues(chainid1));

    /* loop over the residues in chain1 */
    nlines = 0;
    for(r1=0; r1<nres1; r1++) {
      nneighbors = 0;
      resindex1 = res1[r1]->ResIndex();
      diststring = "";

      /* zero the cumulative probabilities for each residue type */
      for(restype=0; restype<UnknownResType; restype++) {
	resprob[restype] = 1.0;
	resgofr[restype] = 1.0;
      }

      /* loop over the second chain ID's */
      interaction = subset.GetIntChains(chainid1);
      for(c2=0; c2<nchains2; c2++) {
	chainid2 = subset.ChainID(c2);
	/* check to see if this chain-chain pair is interacting */
	if (!interaction[c2]) {continue;}

	/* warn if the second chain type is unknown */
	if (subset.ChkUnknown(chainid2)) {
	  cout << "ERROR: both chains in interaction pair have unknown "
	       << "residue type -- not prepared for this." << endl;
	  exit(-1);
	}

	/* get the interaction residues in chain 2 */
	chainid2 = subset.ChainID(c2);
	nres2 = sys.GetResidues(chainid2,subset.GetResidues(chainid2),
			      res2,subset.NResidues(chainid2));

	//      cout << " evaluating pair probabilities for chain '" << chainid1
	//	   << "' with chain '" << chainid2 << "'" << endl;

	for(r2=0; r2<nres2; r2++) {
	  resindex2 = res2[r2]->ResIndex();

	  distsum = 0.0;
	  neighbor = false;
	  for(restype=0; restype<UnknownResType; restype++) {
		 
	    /* loop over residue pair interaction types */
	    for(i=0; i<npairtypes; i++) {
	    //	    for(i=0; i<1; i++) {
	      //	    cout << chainid1 << resindex1 << "  " 
	      //		 << chainid2 << resindex2 << endl;
	      success = (intset[i])->PairProb(res1[r1],restype,res2[r2],
					      resindex2,prob);
	      //	      if (prob[1] > 9.0) {continue;}

	      //#ifdef FEEDBACK	      
	      cout << chainid1 << ":" << res1[r1]->ResNo() 
		   << "(" << res1[r1]->ResType() 
		   << ")(" << restype << "=" 
		   << AAIntType[restype] << ") -- " << chainid2 << ":" 
		   << res2[r2]->ResNo() 
		   << "(" << res2[r2]->ResType() << ")"
		   << "  int: " << i << endl;
	      cout << "returned: g(r)=" << prob[0] 
		   << "  dist=" << prob[1]
		   << "  success=" << success << endl;
	      //#endif

	      if (success < 0) {
		cout << "WARNING: pair probability evaluation failed" << endl;
		cout << chainid1 << ":" << r1 << "(" << restype 
		     << ") -- " << chainid2 << ":" << r2 << endl;
		exit(-1);
	      }

	      /* record information if evaluation was successful */
	      if (success > 1) {continue;}
	      resgofr[restype] *= prob[0];
	      if (restype == 0) {distsum += prob[1];}
	      neighbor = true;
	    }

	    /* evaluate probability due only to pair identity */
	    if (neighbor) {
	      resprob[restype] = (intset[0])->ResResPairProb(res1[r1],restype,
							     res2[r2],resindex2);
	    //#ifdef FEEDBACK	      
	    cout << chainid1 << ":" << res1[r1]->ResNo() 
		 << "(" << res1[r1]->ResType() 
		 << ")(" << restype << "=" 
		 << AAIntType[restype] << ") -- " << chainid2 << ":" 
		 << res2[r2]->ResNo() 
		 << "(" << res2[r2]->ResType() << ")"
		 << "  final probabilities: " << endl;
	    cout << "Pij = chi_ij*g(r) = " << resprob[restype]*resgofr[restype]
		 << "  chi_ij=" << resprob[restype]
		 << "  g(r)=" << resgofr[restype]
		 << "  dist=" << prob[1] << endl;
	    //#endif
	    }

	  } /* residue1 type variation loop */
	  
	  /* calculate average distance to this residue2 and make string */
	  if (neighbor) {
	    nneighbors += 1;
	    sprintf(cstring,"%s%s(%s)@%s ",chainid2.c_str(),
		    (res2[r2]->ResNo()).c_str(),
		    (res2[r2]->ResType()).c_str(),
		    (Num(distsum/npairtypes)).c_str());
	    diststring += cstring;
	  }

	}  /* residue2 loop */
      }    /* chain2 loop */  

      /* calculate probabilities to find res type in residue1 position */
      resprobsum = 0.0; resgofrsum = 0.0;
      resprob[4] = 0.0;
      resgofr[4] = 0.0;
      for(restype=0; restype<UnknownResType; restype++) {
	resprobsum += resprob[restype];
	resgofrsum += resgofr[restype];
      }

      /* create the feedback information */
      sprintf(cstring,"%s%s(%s): ",chainid1.c_str(),
	      (res1[r1]->ResNo()).c_str(),(res1[r1]->ResType()).c_str());
      for(restype=0; restype<UnknownResType; restype++) {
	if ((resprobsum > 0.0)&&(resgofrsum > 0.0)) {
	  sprintf(cstring,"%s%s(%s:%s%s) ",cstring,(AAIntType[restype]).c_str(),
		  (Num(resgofr[restype]/resgofrsum * 100.0)).c_str(),
		  (Num(resprob[restype]/resprobsum * 100.0)).c_str(),"%");
	} else {
	  sprintf(cstring,"%s%s(%s:%s%s) ",cstring,(AAIntType[restype]).c_str(),
		  (Num(0.0)).c_str(),(Num(0.0)).c_str(),"%");
	}
      }
      if (diststring == "") {diststring = "no neighbors";}
      sprintf(cstring,"%s # %d distances: %s",cstring,
	      nneighbors,diststring.c_str());
      output[nlines] = cstring;
      nlines += 1;
      if (nlines >= maxlines) {
	cout << "ERROR: string array passed to forcefield.EvalAlphabet was "
	     << "not large enough" << endl;
	exit(-1);
      }
      cout << nneighbors << " " << cstring << endl;

    } /* residue1 loop */

  } /* chain1 loop */  
  exit(0);
  return (nlines);
}

/*--------------------------------------------------------------------*/
/* Optimize the unknown residue identities in each chain using the */
/* probability of each amino acid type being at each position. */
/* The alternative sequence is stored in the internal subset class */
/* Requires:  sys -- system */
/*--------------------------------------------------------------------*/
void CForcefield::OptimizeSeq(CSystem& sys) {
  int          i, restype, r1, r2, c1, c2, nchains, nneighbors, mostlikely;
  int          nres1, nres2, resindex1, resindex2, nchains2, success;
  bool         neighbor;
  bool*        interaction;
  double       resprobsum,resgofrsum,maxprob;
  double       resprob[NIntResTypes],resgofr[NIntResTypes],prob[UnknownResType];
  double*      p;
  string       chainid1, chainid2;
  CResidue*    res1[subset.MaxNResidues()];
  CResidue*    res2[subset.MaxNResidues()];

  if (!intsubset) {
    cout << "subset must be defined before alphabets can be evaluated" << endl;
    exit(-1);
  }
  if (subset.NUnknown() == 0) {
    cout << "ERROR: forcefield.OptimizeSeq called without any unknown "
	 << "residue type positions being defined in subset" << endl;
    exit(-1);  
  }

  p = new double[4];

  nchains = subset.NUnknown();
  nchains2 = subset.NChainTypes();
  for(c1=0; c1<nchains; c1++) {
    chainid1 = subset.UnknownChainID(c1);
    nres1 = sys.GetResidues(chainid1,subset.GetResidues(chainid1),
			    res1,subset.NResidues(chainid1));

    /* loop over the residues in chain1 */
    for(r1=0; r1<nres1; r1++) {
      nneighbors = 0;
      resindex1 = res1[r1]->ResIndex();

      /* zero the cumulative probabilities for each residue type */
      for(restype=0; restype<UnknownResType; restype++) {
	resprob[restype] = 1.0;
	resgofr[restype] = 1.0;
      }

      /* loop over the second chain ID's */
      interaction = subset.GetIntChains(chainid1);
      for(c2=0; c2<nchains2; c2++) {
	chainid2 = subset.ChainID(c2);
	/* check to see if this chain-chain pair is interacting */
	if (!interaction[c2]) {continue;}

	/* warn if the second chain type is unknown */
	if (subset.ChkUnknown(chainid2)) {
	  cout << "ERROR: both chains in interaction pair have unknown "
	       << "residue type -- not prepared for this." << endl;
	  exit(-1);
	}

	/* get the interaction residues in chain 2 */
	chainid2 = subset.ChainID(c2);
	nres2 = sys.GetResidues(chainid2,subset.GetResidues(chainid2),
			      res2,subset.NResidues(chainid2));

	//      cout << " evaluating pair probabilities for chain '" << chainid1
	//	   << "' with chain '" << chainid2 << "'" << endl;

	for(r2=0; r2<nres2; r2++) {
	  resindex2 = res2[r2]->ResIndex();

	  neighbor = false;
	  for(restype=0; restype<UnknownResType; restype++) {
		 
	    /* loop over residue pair interaction types */
	    for(i=0; i<npairtypes; i++) {
	      //	    cout << chainid1 << resindex1 << "  " 
	      //		 << chainid2 << resindex2 << endl;
	      success = (intset[i])->PairProb(res1[r1],restype,res2[r2],
					      resindex2,p);

#ifdef FEEDBACK	      
	      cout << chainid1 << ":" << res1[r1]->ResNo() 
		   << "(" << res1[r1]->ResType() 
		   << ")(" << restype << "=" 
		   << AAIntType[restype] << ") -- " << chainid2 << ":" 
		   << res2[r2]->ResNo() 
		   << "(" << res2[r2]->ResType() << ")"
		   << "  int: " << i << endl;
	      cout << "returned: g(r)=" << p[0] 
		   << "  dist=" << p[1]
		   << "  success=" << success << endl;
#endif

	      if (success < 0) {
		cout << "WARNING: pair probability evaluation failed" << endl;
		cout << chainid1 << ":" << r1 << "(" << restype 
		     << ") -- " << chainid2 << ":" << r2 << endl;
		exit(-1);
	      }

	      /* record information if evaluation was successful */
	      if (success > 1) {continue;}
	      resgofr[restype] *= p[0];
	      neighbor = true;
	    }

	    /* evaluate probability due only to pair identity */
	    if (neighbor) {
	      resprob[restype] = (intset[0])->ResResPairProb(res1[r1],restype,
							     res2[r2],resindex2);
#ifdef FEEDBACK	      
	    cout << chainid1 << ":" << res1[r1]->ResNo() 
		 << "(" << res1[r1]->ResType() 
		 << ")(" << restype << "=" 
		 << AAIntType[restype] << ") -- " << chainid2 << ":" 
		 << res2[r2]->ResNo() 
		 << "(" << res2[r2]->ResType() << ")"
		 << "  final probabilities: " << endl;
	    cout << "Pij = chi_ij*g(r) = " << resprob[restype]*resgofr[restype]
		 << "  chi_ij=" << resprob[restype]
		 << "  g(r)=" << resgofr[restype]
		 << "  dist=" << p[1] << endl;
#endif
	    }

	  } /* residue1 type variation loop */
	  
	}  /* residue2 loop */
      }    /* chain2 loop */  

      /* calculate probabilities to find res type in residue1 position */
      resprobsum = 0.0; resgofrsum = 0.0;
      resprob[4] = 0.0;
      resgofr[4] = 0.0;
      for(restype=0; restype<UnknownResType; restype++) {
	resprobsum += resprob[restype];
	resgofrsum += resgofr[restype];
      }

      /* find the residue type with highest probability */
      mostlikely = -1;
      maxprob = 0.0;
      for(restype=0; restype<UnknownResType; restype++) {
	if ((resprobsum > 0.0)&&(resgofrsum > 0.0)) {
	  //	  prob[restype] = resprob[restype]/resgofrsum;
	  prob[restype] = resgofr[restype]/resgofrsum;
	} else {
	  prob[restype] = 0.0;
	}
	if (prob[restype] > maxprob) {
	  mostlikely = restype;
	  maxprob = prob[restype];
	}
      }

      /* store the residue type information in the subset class */
      if (mostlikely >= 0) {
	subset.SetResType(chainid1,res1[r1]->ResNo(),mostlikely);
      } else {
	cout << "WARNING: Unable to find an optimal residue type for residue '"
	     << res1[r1]->ResNo() << "' in chain '" << chainid1 << "'" << endl;
      }

    } /* residue1 loop */
  } /* chain1 loop */  

}

/*--------------------------------------------------------------------*/
/* Return the alternative sequence for any chains specified as unknown */
/*--------------------------------------------------------------------*/
string CForcefield::AltSeq() {
  int     c,nchains;
  string  chainid,seq = "";

  if (!intsubset) {
    cout << "subset is undefined, cannot return alternative sequence" << endl;
    exit(-1);
  }
  if (subset.NUnknown() == 0) {
    cout << "ERROR: forcefield.AltSeq called without any unknown "
	 << "residue type positions being defined in subset" << endl;
    exit(-1);  
  }

  nchains = subset.NUnknown();
  for(c=0; c<nchains; c++) {
    chainid = subset.UnknownChainID(c);
    seq += subset.ChainSeq(chainid);
    if (c != nchains-1) {seq += "|";}
  }
  return(seq);
}

/*--------------------------------------------------------------------*/
/* Turn specific chain-chain interactions on/off in the subset class */
/* Requires:  chainid1 -- ID of first chain */
/*            chainid2 -- ID of second chain */
/*            flag -- new value of interaction true/false flag */
/*--------------------------------------------------------------------*/
void CForcefield::SetChainChainInt(const string chainid1, 
				   const string chainid2, bool flag) {
  subset.SetChainChainInt(chainid1, chainid2, flag);
}

/*--------------------------------------------------------------------*/
/* Compare the sequence stored in the internal subset class to that */
/* from the system class */
/* Requires:  sys -- system */
/*            verbose -- if true, use column format with residue no's */
/* NOTE: I haven't implemented the column format form yet */
/*--------------------------------------------------------------------*/
void CForcefield::CompareSeqs(CSystem& sys, const bool) {
  int          r, c, nres;
  int          nchains = subset.NChainTypes();
  int          maxres = subset.MaxNResidues();
  string       chainid, resno;
  string       chainseq1, chainseq2;
  CResidue*    res[maxres];

  if (!intsubset) {
    cout << "subset must be defined before sequences can be compared" << endl;
    exit(-1);
  }
  if (subset.NUnknown() == 0) {
    cout << "ERROR: forcefield.CompareSeqs called without any unknown "
	 << "residue type positions being defined in subset" << endl;
    exit(-1);  
  }

  for(c=0;c<nchains;c++) {
    chainseq1 = "";
    chainid = subset.ChainID(c);
    nres = sys.GetResidues(chainid,subset.GetResidues(chainid),res,
			   subset.NResidues(chainid));
    for(r=0;r<nres;r++) {
      chainseq1 += res[r]->ResSymbol();
    }
    chainseq2 = subset.ChainSeq(chainid);

    cout << chainid << "(WT): " << chainseq1 << endl;
    cout << chainid << ":     " << chainseq2 << endl;

    /* add *'s for identical residues */
    cout << "       ";
    for(r=0;r<nres;r++) {
      if (chainseq1.substr(r,1) == chainseq2.substr(r,1)) {
	cout << "*";
      } else {
	cout << " ";
      }
    }
    cout << endl;
  }
}

/*--------------------------------------------------------------------*/
/* Set the sequence stored in the internal subset class from a string */
/* Requires:  sys -- system */
/*            chainid -- chain identifier */
/*            chainseq -- chain sequence */
/*--------------------------------------------------------------------*/
void CForcefield::SetChainSeq(CSystem& sys, const string chainid, 
			      const string chainseq) {
  int          r, nres;
  int          maxres = subset.MaxNResidues(), resindx[maxres];
  string       resno, subseq = "";
  
  if (!intsubset) {
    cout << "subset must be defined before sequence can be set" << endl;
    exit(-1);
  }

  /* get a list of residue numbers in the subset */
  nres = sys.GetResidueIndices(chainid,subset.GetResidues(chainid),
			       resindx,subset.NResidues(chainid));
			       
  /* extract the relevant subsequence from the passed sequence */
  int seqlen = chainseq.length();
  for(r=0;r<nres;r++) {
    if (resindx[r] > seqlen) {
      cout << "ERROR: (forcefield.SetChainSeq) residue index exceeds length "
	   << "of chain '" << chainid << "'" << endl;
      exit(-1);
    }
    subseq += chainseq.substr(resindx[r],1);
  }
  subset.SetAltSeq(chainid,subseq);

}

/*--------------------------------------------------------------------*/
/* A display routine to dump forcefield parameters to screen */
/* Requires:  indent -- number of spaces to indent */
/*--------------------------------------------------------------------*/
void CForcefield::Display(int indent) {
  int   i;
  char  spacing[100];

  for(i=0;i<indent;i++) {spacing[i] = ' ';}
  spacing[indent] = '\0';

  cout << spacing << "Forcefield type: " << fftype << endl;
  cout << spacing << " with " << npairtypes << " interaction type(s)" << endl;
  if (useunk) {
    cout << spacing << " will use unknown residue type definitions" << endl;
  } else {
    cout << spacing << " NOT using unknown residue type definitions" << endl;
  }
  if (evalchipot) {
    cout << spacing << " will evaluate chi-derived res-res potentials" << endl;
    cout << spacing << "  for all residue pairs with separations within " 
	 << chidistcutoff << " Angstroms" << endl;
  } else {
    cout << spacing << " NOT evaluating chi-derived res-res potentials" << endl;
  }

  /* Loop through the residue pair interaction types */
  for(i=0; i<npairtypes; i++) {
    intset[i]->Display(indent + 1);
  }

  if (intsubset) {
    if (subsetseq) {
      cout << spacing << "Will use alternative sequence in subset" << endl;
    } else {
      cout << spacing << "Will NOT use alternative sequence in subset" << endl;
    }
    cout << spacing << "Interaction subset information: " << endl;
    subset.Display(indent + 1);
  }
  
}
