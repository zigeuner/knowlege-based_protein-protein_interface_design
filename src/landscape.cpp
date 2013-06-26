// Functions for the landscape class.  The purpose of this class is to
// facilitate the analysis of a potential energy landscape by recording
// order parameters and potentials.
// 
// At the moment, it writes the minimum energy configuration for each
// RMSD bin to file.  The idea is to tabulate the energy landscape for
// this single order parameter and get a set of representative structures

#include <iostream>
#include "include/landscape.h"

using std::string;
using std::ofstream;
using std::cout;
using std::ios;
using std::endl;

/*--------------------------------------------------------------------*/
/* Requires:  minbin -- minimum parameter value to bin */
/*            maxbin -- maximum parameter value to bin */
/*            numbins -- number of bins */
/*--------------------------------------------------------------------*/
CLandscape::CLandscape(const double min, const double max, const int numbins) {
  minbin = min;
  maxbin = max;
  nbins = numbins;
  binsize = (maxbin-minbin)/nbins;

  /* allocate memory for the table and initialize */
  minnrg = new double[nbins];
  for(int b=0;b<nbins;b++) {
    minnrg[b] = 1e100;
  }
}

/*--------------------------------------------------------------------*/
/* Compare parameter and configuration energy to table value for that */
/* parameter range */
/* Requires:  sys -- system class */
/*            param -- parameter value (same as bin variable, RMSD) */
/*            potnrg -- potential energy */
/*            maxdump -- maximum param value for which structure dumped */
/*            filebase -- base name for output file */
/*--------------------------------------------------------------------*/
void CLandscape::Process(CSystem& sys, double param, double potnrg, 
			 double maxdump, string filebase) {
  int      bin;
  char     file[50];
  ofstream pdbfile;

  /* determine if this is a new minimum energy for this parameter range */
  bin = Bin(param);
  if (potnrg >= minnrg[bin]) {return;}
  minnrg[bin] = potnrg;

  /* skip dump if param value is too high */
  if (param > maxdump) {return;}

  /* open structure file for output*/
  sprintf(file,"%s_bin%03d.pdb",filebase.c_str(),bin);
  pdbfile.open(file, ios::out);
  if(!pdbfile) {
    cout << endl << " Error: cannot open file '" << file << "'" << endl;
    exit(-1);
  }  

  /* write system configuration to file */
  pdbfile << "REMARK  RMSD= " << Num(param) << "  PotNrg= " 
	  << Num(potnrg) << endl;
  sys.WritePDB(pdbfile);
  pdbfile.close();
}

/*--------------------------------------------------------------------*/
/* Write the landscape information to a file for graphing */
/* Requires:  filename -- output file name */
/*--------------------------------------------------------------------*/
void CLandscape::Write(string filename) {
  int      minimumbin;
  double   min,dist;
  char     file[100];
  ofstream output;

  /* open file for output*/
  sprintf(file,"%s",filename.c_str());
  output.open(file, ios::out);
  if(!output) {
    cout << endl << " Error: cannot open file '" << file << "'" << endl;
    exit(-1);
  }  

  /* find the minimum value */
  min = 1000;
  for(int b=0;b<nbins;b++) {
    if (minnrg[b] < min) {
      min = minnrg[b];
      minimumbin = b;
    }
  }
  dist = minbin + minimumbin*binsize + binsize/2.0;

  /* write header */
  output << "# landscape information" << endl;
  output << "# global minimum energy found in bin centered at: " 
	 << Num(dist) << " with energy: " << Num(min) << endl;
  output << "# columns are: bin number, min bin value, max bin value, " 
	 << "minimum energy" << endl;

  /* write information */
  for(int b=0;b<nbins;b++) {
    min = minbin + b*binsize;
    if (minnrg[b] == 1e100) {
      output << b << "  " << Num(min) << "  " 
	     << Num(min+binsize) << "  NA" << endl;
    } else {
      output << b << "  " << Num(min) << "  " << Num(min+binsize) 
	     << "  " << Num(minnrg[b]) << endl;
    }
  }
  output.close();
}

/*--------------------------------------------------------------------*/
/* A display routine to dump landscape information to screen */
/* Requires:  indent -- number of spaces to indent */
/*--------------------------------------------------------------------*/
void CLandscape::Display(const int indent) {
  char    spacing[100];
  string  chainseq;

  for(int i=0;i<indent;i++) {spacing[i] = ' ';}
  spacing[indent] = '\0';

  cout << "landscape information: " << endl;
  for(int i=0; i<20; i++) {
  }
  
}
