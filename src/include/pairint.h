// This file provides the core routines for performing lookup table
// model energy and force evaluations.  Routines here are inline
// to allow the compiler to put them into the loop in forcefield.

#include <cmath>
#include <string>
#include <iostream>
#include <stdlib.h>
#include "config.h"
#include "ffbase.h"
#include "stringutils.h"

class CPairInt : public CFFBase {
 private:
  int          nentries;
  double       firstsep, finalsep, sepinc;
  std::string  id;
  double*      table;

 public:
  /* Constructors */
  CPairInt() 
    : nentries(0), firstsep(0.0), sepinc(0.0), id(""), table(0) {};
  CPairInt(std::string filename, std::string label, bool verbose, bool warn) {
    Init(filename,label,verbose,warn);
  }
  
  /* Destructor */
  ~CPairInt() { 
    /* cout << " CPairInt destructor.." << endl; */ 
    delete[] table;
  }

  /* Functions declared in .cpp file */
  void Init(std::string, std::string, bool, bool);
  void Check(const double, int&, int&, bool);
  void Display(int);

  std::string ID() {return(id);}

  /*------------------------------------------------------------------------*/
  /* Energy routine */
  /* Performs linear interpolation between two table entries */
  /* Returns false if there was an error during evaluation */
  /* Since the interpolated numbers are g(r)'s, this does not return */
  /* a potential directly, but must be converted by kbtypeint */
  /* Requires:  sepvec -- Rij separation vector */
  /*            potential -- returned potential */
  /* Returns: -1 -- failed (can't currently happen) */
  /*           0 -- under range, returns zero */
  /*           1 -- ok, within range */
  /*           2 -- above range, returned unity */
  /*------------------------------------------------------------------------*/
  int Potential(const CVector3& sepvec, double& potential) {
    int     bin;
    double  rawbin;

    potential = 0.0;        /* default for small r */
    double r = sepvec.Length();
    //    cout << "pairint.h:Potential separation distance: " << r << std::endl;
    if (r > finalsep) {potential = 1.0; return(2);}
    if (r >= firstsep) {
      rawbin = (r - firstsep)/sepinc;
      bin = int(rawbin);
      if (bin == nentries) {potential = table[bin];}
      potential = table[bin] + 
	(rawbin - bin)*(table[bin+1] - table[bin]);
      //      std::cout << "distance= " << r << " bin1= " << bin << " bin2= " 
      //		<< bin + 1 << " maxbin= " << nentries 
      //		<< " g(r)= " << potential << std::endl;
      //      std::cout << "btwn " << table[bin] << " and " 
      //                << table[bin+1] << std::endl;
      //      std::cout << "rawbin " << rawbin << "  bin " << bin << std::endl;
      //      std::cout << "sepinc " << sepinc << std::endl;
      //      Display(0);
      return(1);
    }
    return(0);
  }

  /*------------------------------------------------------------------------*/
  /* Force routine */
  /* Requires:  sepvec -- Rij separation vector */
  /*            potential -- returned potential */
  /*            force -- returned force vector */
  /*------------------------------------------------------------------------*/
  bool Force(const CVector3& sepvec, double potential, CVector3& force) {
    std::cout << "Force routine not yet ready" << std::endl;
    exit(-1);    
  }
};

