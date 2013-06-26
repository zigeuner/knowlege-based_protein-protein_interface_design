/*
   Blah blah blah
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

void parseflags(int, char**, char*);

int main (int argc, char *argv[]) {
  char       pdbfile[50], file[50];
  string     atomtype1, atomtype2, chain1, chain2, label;  
  ofstream   file1, file2;

  /* Read run parameters from command line */
  parseflags(argc, argv, pdbfile);

  cout << "reading from file: " << pdbfile << endl;

  /* initialize the system configuration */
  CSystem sys(pdbfile,1);

  /* Open the ca-ca file for output*/
  sprintf(file,"pairs_ca-ca.txt");
  file1.open(file, ios::out);
  if(!pdbfile) {
    cout << endl << " Error: cannot open file '" << file << 
      "'" << endl;
    exit(-1);
  }

  label = pdbfile + chain1 + chain2;

  /* get Calpha-Calpha pairs on VL-VH interface */
  atomtype1 = "CA";
  atomtype2 = "CA";
  chain1 = "L";
  chain2 = "H";
  sys.Pairs(chain1,atomtype1,chain2,atomtype2,10.0,label,file1);
  file1.close();

  /* Open the cb-cb file for output*/
  sprintf(file,"pairs_cb-cb.txt");
  file2.open(file, ios::out);
  if(!pdbfile) {
    cout << endl << " Error: cannot open file '" << file << 
      "'" << endl;
    exit(-1);
  }

  /* get Cbeta-Cbeta pairs on VL-VH interface */
  atomtype1 = "CB";
  atomtype2 = "CB";
  chain1 = "L";
  chain2 = "H";
  sys.Pairs(chain1,atomtype1,chain2,atomtype2,10.0,label,file2);
  file2.close();

  return 0;
}


void parseflags(int argc, char* argv[], char* pdbfile) {

  strcpy(pdbfile, "NULL");

  if( argc >= 2 ) {
    for(int i=1; i<argc; i++) {
      if( *(argv[i]) == '-' ) {
        switch( *(argv[i]+1) ) {
	case 'p' : 
	  strcpy(pdbfile, argv[i+1]);
	  break;
	default :
	  cout << " Unknown flag " << argv[i] << endl << endl;
	  flaglist(argv[0]);
	}
      } else {
	strcpy(pdbfile, argv[i]);	
      }
    }
  }

}

void flaglist(char* str) {
  cout << endl 
    << " Usage: " << str
      << " [-dcxp]"
	<< " [-d directory name]"
	  << " [-p parameter file]"
	    << endl << endl;
  cout << " Supported flags: " << endl << endl;
  cout << "  d - ";
  cout << "dummy\n";
  cout << endl;

  exit(1);
}
