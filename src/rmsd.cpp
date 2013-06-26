/*
   This driver reads in two pdb files and calculates an RMSD
*/

#include <cstdio>
#include <ctime>
#include <unistd.h>
#include <dirent.h>
#include <iostream>
#include "include/system.h"
#include "include/forcefield.h"
#include "include/move.h"
#include "include/landscape.h"

using std::cout;
using std::endl;
using std::ios;
using std::string;
using std::ofstream;
using std::ifstream;

void parseflags(int, char**, string&, string&, string&);

int main (int argc, char *argv[]) { 
  double     rmsd;
  string     file1,file2,residuesfile;

  /* read run parameters from command line */
  parseflags(argc, argv, file1, file2, residuesfile);

  /* make sure the first system pdb file exists, if so initialize from it */
  ifstream   test;    
  test.open(file1.c_str(), ios::in);
  if(!test) {
    cout << "ERROR: cannot open ref system configuration file '" 
	 << file1 << endl;
    exit(-1);
  } else {
    test.close();
  }
  cout << "reading ref system configuration from file: " << file1 << endl;
  CSystem refsys(file1.c_str(),0);

  /* make sure the second system pdb file exists, if so initialize from it */
  test.open(file2.c_str(), ios::in);
  if(!test) {
    cout << "ERROR: cannot open system configuration file '" << file2 << endl;
    exit(-1);
  } else {
    test.close();
  }
  cout << "reading system configuration from file: " << file2 << endl;
  CSystem sys(file2.c_str(),0);

  /* define the interface residues */
  CSubset interface(residuesfile,0);

  rmsd = sys.RMSD(refsys,interface,"CA");
  cout << Num(rmsd) << endl;
  cout << "Done" << endl;

  return 0;
}

void parseflags(int argc, char* argv[], string& file1, string& file2,
		string& residuesfile) {
  residuesfile = "NULL";
  file1 = "NULL";
  file2 = "NULL";

  if( argc >= 2 ) {
    for(int i=1; i<argc; i++) {
      if( *(argv[i]) == '-' ) {
        switch( *(argv[i]+1) ) {
	case 'r' :
	  residuesfile = argv[i+1];
	  i++;
	  break;
	default :
	  cout << " Unknown flag " << argv[i] << endl << endl;
	  flaglist(argv[0]);
	}
      } else {
	if (i == argc-1) {
	  file2 = argv[i];
	} else {
	  file1 = argv[i];
	}
      }
    }
  }
  if ((file1 == "NULL")||(file2 == "NULL")||(residuesfile == "NULL")) {
    cout << "ERROR: must specify residues defn and input filenames" << endl;
    flaglist(argv[0]);
    exit(-1);
  }

}

void flaglist(char* str) {
  cout << endl 
       << " Usage: " << str
       << " -r [interface file] "
       <<" [system .pdb or .xyz filename] [system .pdb or .xyz filename]"
       << endl;
  cout << " Supported flags: " << endl;
  cout << "  r - ";
  cout << "specify filename containing interacting residues\n";
  cout << endl;

  exit(1);
}
