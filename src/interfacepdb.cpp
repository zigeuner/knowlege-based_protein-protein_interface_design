/*
   This driver reads in an interface specification file and a pdb file,
   reduces the pdb file to just residues on the interfaces and writes
   the new system .pdb file.
*/

#include <cstdio>
#include <ctime>
#include <unistd.h>
#include <dirent.h>
#include <iostream>
#include "include/system.h"

using std::cout;
using std::endl;
using std::ios;
using std::string;
using std::ofstream;
using std::ifstream;

void parseflags(int, char**, long&, double&, int&, string&, string&,
                string&, string&, string&, int&, bool&);

int main (int argc, char *argv[]) { 
  int        niter, nabchains, fftype;
  bool       caonly;
  double     temperature;
  long       startseed;
  char       file[50];
  string     ffname, commentline;
  string     atomtype1, atomtype2, chain1, chain2, line, comment, abchain[100];  
  string     potentialsdir, abdefn, agdefn, pdbfile, core, residuesfile;
  ofstream   output, statsfile, movie;
  CSystem    sys;

  /* set the core name of the output files */
  core = "poses";

  /* read run parameters from command line */
  parseflags(argc, argv, startseed, temperature, niter, pdbfile, abdefn, agdefn, 
	     residuesfile, potentialsdir, fftype, caonly);

  /* check input */
  nabchains = Split(abdefn,',',abchain,100);
  if (nabchains > 2) {
    cout << "WARNING: what's up, there's more than 2 antibody chains?" << endl;
    exit(-1);
  }

  /* make sure the system pdb file exists, if so initialize from it */
  sprintf(file,"%s",pdbfile.c_str());
  ifstream   test;    
  test.open(file, ios::in);
  if(!test) {
    cout << "ERROR: cannot open system configuration file '" << file << endl;
    exit(-1);
  } else {
    test.close();
  }
  cout << "reading system configuration from file: " << pdbfile << endl;
  CSystem fullsys(file,0);

  /* define the interface residues and check them for existence */
  CSubset interface(residuesfile,true);
  fullsys.CheckSubset(interface);

  /* make the reduced system configuration */
  string atomtypes[3] = {"CA","CB","MN"};
  //  sys.CopyReduce(fullsys,agdefn,atomtypes,3);
  //  sys.CopySubset(fullsys,interface,agdefn,atomtypes,3);
  sys.CopySubset(fullsys,interface,"",atomtypes,3);
  //  sys.CopySubset(fullsys,interface,"",atomtypes,3);
  sys.WritePDB("interface.pdb");

  cout << "Done" << endl;

  return 0;
}

void parseflags(int argc, char* argv[], long& startseed, double&
		temperature, int& nsteps, string& inputfile, string&
		abdefn, string& agdefn, string& residuesfile, string&
		potentialsdir, int& fftype, bool& caonly) {
  startseed = 0;
  temperature = 300.0;
  nsteps = 1000;
  inputfile = "NULL";
  abdefn = "NULL";
  agdefn = "NULL";
  residuesfile = "NULL";
  potentialsdir = ".";
  fftype = 0; /* 0= knowledge-based FF with normal residue types */
              /* 1= knowledge-based FF with unknown type Ab residues */
              /* 2= LJ potentials */
  caonly = false;   /* if true, use only CA-CA pairs (no CB-CB) */

  if( argc >= 2 ) {
    for(int i=1; i<argc; i++) {
      if( *(argv[i]) == '-' ) {
        switch( *(argv[i]+1) ) {
	case 'd' : 
	  potentialsdir = argv[i+1];
	  i++;
	  break;
	case 'b' :
	  abdefn = argv[i+1];
	  i++;
	  break;
	case 'g' :
	  agdefn = argv[i+1];
	  i++;
	  break;
	case 'r' :
	  residuesfile = argv[i+1];
	  i++;
	  break;
	case 's' :
	  startseed = ToInt(argv[i+1]);
	  i++;
	  break;
	case 'n' :
	  nsteps = ToInt(argv[i+1]);
	  i++;
	  break;
	case 'T' :
	  temperature = ToDbl(argv[i+1]);
	  i++;
	  break;
	case 'u' :
	  fftype = 1;
	  break;
	case 'l' :
	  fftype = 2;
	  break;
	case 'a' :
	  caonly = true;
	  break;
	default :
	  cout << " Unknown flag " << argv[i] << endl << endl;
	  flaglist(argv[0]);
	}
      } else {
	inputfile = argv[i];	
      }
    }
  }
  if (inputfile == "NULL") {
    cout << "ERROR: must specify -r [interface file] and input filename" << endl;
    flaglist(argv[0]);
    exit(-1);
  }

  /* Some feedback for the user */
  //  cout << "Antibody is composed of chains: '" << abdefn << "'" << endl;
  //  cout << "Antigen is composed of chains: '" << agdefn << "'" << endl;

}

void flaglist(char* str) {
  cout << endl 
       << " Usage: " << str
       << " -b [chain ID's] -g [chain ID's] [-fdrsnTula] "
       <<" [system .pdb or .xyz filename]"
       << endl;
  cout << " Supported flags: " << endl;
  cout << "  d - ";
  cout << "specify directory from which to read res-res potentials\n";
  cout << "  b - ";
  cout << "specify definition of antibody, example: [-b A,B]\n";
  cout << "  g - ";
  cout << "specify definition of antigen, example: [-g C,D]\n";
  cout << "  r - ";
  cout << "specify filename containing interacting residues\n";
  cout << "  s - ";
  cout << "specify starting random number seed\n";
  cout << "  n - ";
  cout << "specify number of MC steps in the run\n";
  cout << "  T - ";
  cout << "specify starting temperature in Kelvin\n";
  cout << "  u - ";
  cout << "use unknown residue type for antibody residues\n";
  cout << "  l - ";
  cout << "use generic LJ potentials \n";
  cout << "  a - ";
  cout << "use only Calpha-Calpha interaction type\n";
  cout << endl;

  exit(1);
}
