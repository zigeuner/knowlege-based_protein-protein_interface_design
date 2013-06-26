#ifndef __LANDSCAPE_H
#define __LANDSCAPE_H
// The landscape class.  The purpose of this class is to
// facilitate the analysis of a potential energy landscape by recording
// order parameters and potentials.
// 
// At the moment, it writes the minimum energy configuration for each
// RMSD bin to file.  The idea is to tabulate the energy landscape for
// this single order parameter and get a set of representative structures

#include <string>
#include "config.h"
#include "stringutils.h"
#include "system.h"

class CLandscape {
 private:
  /* information for the pseudo-histogram */
  int      nbins;
  double   minbin,maxbin,binsize;
  double   *minnrg;

 public:
  /* constructors */
  CLandscape() : nbins(0), minbin(0.0), maxbin(0.0) {};
  CLandscape(const double, const double, const int);

  /* Destructor */
  ~CLandscape() {
    delete[] minnrg;
  };

  /* Functions declared in the .ccp file */
  void Process(CSystem&, double, double, double, std::string);
  void Write(std::string);
  void Display(const int);

  /* Locally declared functions */
  int Bin(double param) {
    int bin = int((param - minbin)/binsize);
    if (bin >= nbins) {
      //      std::cout << "WARNING: value exceeds binned parameter range (" 
      //		<< minbin << "->" << maxbin << "): " << param 
      //		<< ", returning maximum bin " << std::endl;
      return(nbins-1);
    }
    return(bin);
  }

};

#endif

