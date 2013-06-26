#ifndef __SUBSET_H
#define __SUBSET_H
// This subset class is responsible for containing system subset specifications
// I consider this a temporary structure, because it doesn't seem very general

#include <string>
#include <cstdlib>
#include "config.h"
#include "stringutils.h"
 
class CSubset {
 private:
  int           nchaintypes,maxnresidues,totalresidues,nunknown;
  bool*         chainunknown;  /* true if chain has unknown residue types */
  int*          nresidues;     /* Nresidues in each chain */
  bool**        intmtx;        /* true if two chain types interact */
  std::string*  unknownchain;  /* chain ID's for unknown residue type chains*/
  std::string*  intchain;      /* a list of interacting chain ID's */
  std::string** intres;        /* list of interacting residues for each chain */
  /* variables for specifying an alternative sequence */
  bool          altseq;        /* true if an alternative sequence is present */
  int**         altresindx;    /* integer residue type (alternative seq) */
  std::string** altres;        /* 1-letter residue type (alternative seq) */

 public:
  /* constructors */
  CSubset() 
    : nchaintypes(0), maxnresidues(0), totalresidues(0), 
      nunknown(0), chainunknown(0), nresidues(0), 
      intmtx(0), unknownchain(0), intchain(0), intres(0),
      altseq(false), altresindx(0), altres(0) {
    //    std::cout << "calling CSubset constructor" << std::endl;
  }
  CSubset(const std::string, bool);
  CSubset(const std::string, std::string, bool);
  CSubset(std::string*, const int, bool);

  /* Destructor */
  ~CSubset() {
    //    std::cout << "calling CSubset destructor" << std::endl;
    if (nunknown > 0) {delete[] chainunknown;}
    delete[] nresidues;
    delete[] intmtx;
    if (nunknown > 0) {delete[] unknownchain;}
    delete[] intchain;
    delete[] intres;
    delete[] altres;
    delete[] altresindx;
  }

  /* Functions for returning the state variables of the object*/
  int MaxNResidues() {return (maxnresidues);}
  int TotalResidues() {return (totalresidues);}
  int NChainTypes() {return (nchaintypes);}
  int NUnknown() {return (nunknown);}
  int NResidues(const std::string chainid) {
    return nresidues[ChainIndex(chainid)];
  }
  bool ChkUnknown(const int index) {return (chainunknown[index]);}
  bool ChkUnknown(const std::string chainid) {
    int index = ChainIndex(chainid);
    if ((index < 0)||(index >= nchaintypes)) {
      std::cout << "subset.h: unexpected, unable to get find chain ID '" 
		<< chainid << "' in unknown chain array" << std::endl;
      exit(-1);
    }
    return (chainunknown[index]);
  }
  std::string ChainID(const int index) {return (intchain[index]);}
  std::string UnknownChainID(const int index) {return (unknownchain[index]);}
  int AltResIndex(const int c, const int r) {return (altresindx[c][r]);}
  std::string AltRes(const int c, const int r) {return (altres[c][r]);}
  std::string*& GetResidues(const std::string chainid) {
    return (intres[ChainIndex(chainid)]);
  }
  bool*& GetIntChains(const std::string chainid) {
    return (intmtx[ChainIndex(chainid)]);
  }
  bool InSubset(std::string chainid) {
    for(int i=0; i<nchaintypes; i++) {
      if (intchain[i] == chainid) {return(true);}
    }
    return(false);
  }
  bool InSubset(std::string chainid, std::string resno) {
    if (!InSubset(chainid)) {return(false);}
    int chainno = ChainIndex(chainid);
    for(int i=0; i<nresidues[chainno]; i++) {
      if (intres[chainno][i] == resno) {return(true);}
    }
    return(false);
  }
  void SetResType(const std::string chainid, const std::string resno, 
		  const int restypeindex) {
    int chainindex = ChainIndex(chainid);
    int resindex = ResIndex(chainindex,resno);
    altresindx[chainindex][resindex] = restypeindex;
    altres[chainindex][resindex] = AAIntType[restypeindex];
  }
  void SetAltSeq(const std::string chainid, const std::string seq) {
    int chainindex = ChainIndex(chainid);
    for(unsigned int r=0;r<seq.length();r++) {
      altresindx[chainindex][r] = AAIntType2Index(seq.substr(r,1));
      altres[chainindex][r] = seq.substr(r,1);
    }
  }

  /* Functions declared in .cpp file */
  void Init(const std::string, bool);
  void Init(std::string*, const int, bool);
  void Init(const std::string, std::string, bool);
  void SetUnknown(const std::string, const bool);
  void SetChainChainInt(const std::string, const std::string, bool);
  int ID2Num(const std::string, const int);
  int Num2ID(const int, std::string&, std::string&);
  std::string ChainSeq(const std::string);
  void Display(int);

  /* Operator overloads */
  CSubset& operator=(const CSubset& old) {
    //    std::cout << "calling CSubset =" << std::endl;
    if (this != &old) {  /* avoid old = old problems */

      nchaintypes = old.nchaintypes;
      nunknown = old.nunknown;
      chainunknown = new bool[nchaintypes];
      if (nunknown > 0) {
	unknownchain = new std::string[nunknown];
      }

      intmtx = new bool*[nchaintypes];
      nresidues = new int[nchaintypes];
      intres = new std::string*[nchaintypes];
      intchain = new std::string[nchaintypes];
      for(int i=0; i<nchaintypes; i++) {
	nresidues[i] = old.nresidues[i];
	intmtx[i] = new bool[nchaintypes];
	intres[i] = new std::string[nresidues[i]];
      }

      altres = new std::string*[nchaintypes];
      altresindx = new int*[nchaintypes];
      for(int i=0; i<nchaintypes; i++) {
	altres[i] = new std::string[old.nresidues[i]];
	altresindx[i] = new int[old.nresidues[i]];
      }

      /* Copy in the old values */
      Copy(old);
    }
    return *this;
  }

  friend CSubset operator-(CSubset& subset1, CSubset& subset2) {
    CSubset  ss;

    ss.nchaintypes = subset1.nchaintypes;
    ss.nunknown = subset1.nunknown;
    ss.chainunknown = new bool[ss.nchaintypes];

    for(int i=0; i<ss.nunknown; i++) {    
      ss.chainunknown[i] = subset1.chainunknown[i];
    }
    if (ss.nunknown > 0) {
      ss.unknownchain = new std::string[ss.nunknown];
      for(int i=0; i<ss.nunknown; i++) {    
	ss.unknownchain[i] = subset1.unknownchain[i];
      }
    }

    ss.intmtx = new bool*[ss.nchaintypes];
    ss.intchain = new std::string[ss.nchaintypes];
    for(int i=0; i<ss.nchaintypes; i++) {
      ss.intchain[i] = subset1.intchain[i];
      ss.intmtx[i] = new bool[ss.nchaintypes];
      for(int j=0; j<ss.nchaintypes; j++) {
	ss.intmtx[i][j] = subset1.intmtx[i][j];
      }
    }

    /* determine if any residues should be removed from subset, size + copy */
    ss.maxnresidues = 0;
    ss.totalresidues = 0;
    ss.altseq = subset1.altseq;
    ss.nresidues = new int[ss.nchaintypes];
    ss.intres = new std::string*[ss.nchaintypes];
    ss.altres = new std::string*[ss.nchaintypes];
    ss.altresindx = new int*[ss.nchaintypes];
    for(int i=0; i<ss.nchaintypes; i++) {
      ss.nresidues[i] = 0;
      for(int j=0; j<subset1.nresidues[i]; j++) {
	if (!(subset2.InSubset(subset1.ChainID(i),subset1.intres[i][j]))) {
	  ss.nresidues[i] += 1;
	} 
      }

      /* size the number of residues */
      ss.intres[i] = new std::string[ss.nresidues[i]];
      ss.altres[i] = new std::string[ss.nresidues[i]];
      ss.altresindx[i] = new int[ss.nresidues[i]];
      ss.totalresidues += ss.nresidues[i];
      if (ss.nresidues[i] > ss.maxnresidues) {ss.maxnresidues = ss.nresidues[i];}
      int count = 0;

      /* copy the interacting residues */
      for(int j=0; j<subset1.nresidues[i]; j++) {
	if (!subset2.InSubset(subset1.ChainID(i),subset1.intres[i][j])) {
	  ss.intres[i][count] = subset1.intres[i][j];
	  ss.altres[i][count] = subset1.altres[i][j];
	  ss.altresindx[i][count] = subset1.altresindx[i][j];
	  count += 1;
	}
      }
    }

    return ss;
  }

  /* Inline functions */

  /*------------------------------------------------------------------------*/
  /* Copies one structure into another, without initializing the first */
  /* IMPORTANT: allocation sizes are assumed to match! */
  /* Requires:  old -- old structure to copy from */
  /*------------------------------------------------------------------------*/
  void Copy(const CSubset& old) {
    int i,j;

    //    std::cout << "calling CSubset copy" << std::endl;
    nchaintypes = old.nchaintypes;
    nunknown = old.nunknown;
    maxnresidues = old.maxnresidues;
    totalresidues = old.totalresidues;

    /* Copy unknown-chain quantities */
    for(i=0; i<nchaintypes; i++) {    
      chainunknown[i] = old.chainunknown[i];
    }
    if (nunknown > 0) {
      for(i=0; i<nunknown; i++) {    
	unknownchain[i] = old.unknownchain[i];
      }
    }

    /* Copy normal quantities */
    for(i=0; i<nchaintypes; i++) {
      nresidues[i] = old.nresidues[i];
      intchain[i] = old.intchain[i];
      for(j=0; j<nchaintypes; j++) {
	intmtx[i][j] = old.intmtx[i][j];
      }
      for(int j=0; j<nresidues[i]; j++) {
	intres[i][j] = old.intres[i][j];
      }
    }

    /* Copy alternative sequence information */
    altseq = old.altseq;
    for(i=0; i<nchaintypes; i++) {
      for(j=0; j<nresidues[i]; j++) {
	altres[i][j] = old.altres[i][j];
	altresindx[i][j] = old.altresindx[i][j];
      }
    }
  }

  /*------------------------------------------------------------------------*/
  /* Returns the internal chain index given a chain ID */
  /* Requires:  chainid -- string chain ID */
  /*------------------------------------------------------------------------*/
  int ChainIndex(const std::string chainid) {
    for(int i=0; i<nchaintypes; i++) {
      if (intchain[i] == chainid) {return(i);}
    }
    return(-1);
    std::cout << "ERROR: (subset.ChainIndex) could not find chain ID '" << 
      chainid << "' in the set of specified interacting residues " << std::endl;
    exit(-1);
  }

  /*------------------------------------------------------------------------*/
  /* Returns the internal residue index given the internal chain index */
  /* and residue number */
  /* Requires:  chainid -- chain internal index */
  /*            resno -- string residue number */
  /*------------------------------------------------------------------------*/
  int ResIndex(const int chainid, const std::string resno) {
    for(int j=0; j<nresidues[chainid]; j++) {
      if (intres[chainid][j] == resno) {return(j);}
    }
    std::cout << "ERROR: (subset.ChainIndex) could not find residue number '" << 
      resno << "' in the set of specified interacting residues " << std::endl;
    exit(-1);
  }

};

#endif
