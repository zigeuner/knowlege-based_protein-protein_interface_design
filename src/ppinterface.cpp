/*
   This driver reads in a file containing a list of structures and chain
   ID's that make up an interface.  Keskin's notation is used to facilitate
   reading their clustering analysis files.
   Example line from input file:
   1axdAB
   This indicates the interface between chains A and B in PDB file 1AXD
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

void parseflags(int, char**, string&, string&, string&, string&, double&);
bool getspec(string&, string, string&, ifstream&);

int main (int argc, char *argv[]) { 
  int        nires;
  char       pdbfile[50], file[50];
  double     cutoff;
  string     inputfile, singlespec, core;
  string     atomtype1, atomtype2, chain1, chain2, line, interface, pdb;  
  string     interfacedir, pdbdir, iresline[1000];
  ifstream   interfacefile, iresfile;
  ofstream   file1, file2, cacafile, cbcbfile, czczfile, anglesfile;

  /* Read run parameters from command line */
  parseflags(argc, argv, inputfile, singlespec, interfacedir, pdbdir, cutoff);
  if (cutoff == 0.0) {cutoff = 12.0;}

  /* Open the interface list file if necessary */
  if (inputfile != "NULL") {
    core = "pairs";
    interfacefile.open(inputfile.c_str(), ios::in);
    if(!interfacefile) {
      cout << endl << " Error: cannot open file '" << inputfile << 
	"'" << endl;
      exit(-1);
    }
    cout << "reading interface list from file: " << inputfile << endl;
  } else {
    core = singlespec.substr(0,6);
  }

  /* Open the ca-ca file for output*/
  sprintf(file,"%s_ca-ca.txt",core.c_str());
  cacafile.open(file, ios::out);
  if(!cacafile) {
    cout << endl << " Error: cannot open file '" << file << 
      "'" << endl;
    exit(-1);
  }

  /* Open the cb-cb file for output*/
  sprintf(file,"%s_cb-cb.txt",core.c_str());
  cbcbfile.open(file, ios::out);
  if(!cbcbfile) {
    cout << endl << " Error: cannot open file '" << file << 
      "'" << endl;
    exit(-1);
  }

  /* Open the cz-cz file for output*/
  sprintf(file,"%s_cz-cz.txt",core.c_str());
  czczfile.open(file, ios::out);
  if(!czczfile) {
    cout << endl << " Error: cannot open file '" << file << 
      "'" << endl;
    exit(-1);
  }

  /* Open the angles etc file for output*/
  sprintf(file,"%s_angles.txt",core.c_str());
  anglesfile.open(file, ios::out);
  if(!anglesfile) {
    cout << endl << " Error: cannot open file '" << file << 
      "'" << endl;
    exit(-1);
  }

  /* Read the interface specification lines and process */
  while (getspec(interface, inputfile, singlespec, interfacefile)) {
    /* parse the interface specification */
    if (interface.length() != 6) {
      cout << "unexpected interface specification '" 
	   << interface << "' skipping" << endl;
      continue;
    }
    cout << "found interface specification '" << interface << "'" << endl;
    pdb = interface.substr(0,4);
    chain1 = interface.substr(4,1);
    chain2 = interface.substr(5,1);

    /* open and read the file that specifies the interface residues */
    sprintf(file,"%s/%s_interface",interfacedir.c_str(),interface.c_str());
    cout << "reading interface residue list for '" << interface << 
      "' from file: " << file << endl;
    if ((nires = ReadFile(file,iresline,1000)) < 0) {
      continue;
    }
    //    cout << "DEBUG at " << __FILE__ << ":" << __LINE__ << endl;
    //    for(int iter=1; iter<nires; iter++) {
    //      cout << iresline[iter] << endl;
    //    }
    //    cout << "DEBUG at " << __FILE__ << ":" << __LINE__ << endl;
    //    exit(0);

    /* read the appropriate pdb file and analyze interface */
    sprintf(pdbfile,"%s/%s.pdb",pdbdir.c_str(),pdb.c_str());
    ifstream   test;    
    test.open(pdbfile, ios::in);
    if(!test) {
      cout << "ERROR: cannot open file '" << pdbfile << 
	"' skipping this interface analysis" << endl;
      continue;
    } else {
      test.close();
    }
    cout << " reading from file: " << pdbfile << endl;
    CSystem sys(pdbfile,0);
    //    CSystem sys(pdbfile,1);
    //    sys.WritePDB("junk.pdb");
    //    exit(0);

    /* analyze Calpha-Calpha pairs on interface */
    cout << "  analyzing Calpha-Calpha pairs on " << chain1 << "-" 
	 << chain2 << " interface of " << pdb << ".pdb ..." << endl;
    cacafile << "# Calpha-Calpha pairs on " << chain1 << "-" 
	     << chain2 << " interface of " << pdb << ".pdb" << endl;
    atomtype1 = "CA";
    atomtype2 = "CA";
    sys.SomePairs(chain1,atomtype1,chain2,atomtype2,cutoff,interface,
		  iresline,nires,cacafile);

    /* analyze Cbeta-Cbeta pairs on interface */
    cout << "  analyzing Cbeta-Cbeta pairs on " << chain1 << "-" 
	 << chain2 << " interface of " << pdb << ".pdb ..." << endl;
    cbcbfile << "# Cbeta-Cbeta pairs on " << chain1 << "-" 
	     << chain2 << " interface of " << pdb << ".pdb" << endl;
    atomtype1 = "CB";
    atomtype2 = "CB";
    sys.SomePairs(chain1,atomtype1,chain2,atomtype2,cutoff,interface,
		  iresline,nires,cbcbfile);

    /* analyze end-carbon pairs on interface */
    cout << "  analyzing Cmax-Cmax pairs on " << chain1 << "-" 
	 << chain2 << " interface of " << pdb << ".pdb ..." << endl;
    czczfile << "# Cmax-Cmax pairs on " << chain1 << "-" 
	     << chain2 << " interface of " << pdb << ".pdb" << endl;
    atomtype1 = "MAXC";
    atomtype2 = "MAXC";
    sys.SomePairs(chain1,atomtype1,chain2,atomtype2,cutoff,interface,
		  iresline,nires,czczfile);

    /* analyze angles for residue pairs on interface */
    cout << "  analyzing angles for residues pairs on " << chain1 << "-" 
	 << chain2 << " interface of " << pdb << ".pdb ..." << endl;
    anglesfile << "# Angle information for residue pairs on " << chain1 << "-" 
	       << chain2 << " interface of " << pdb << ".pdb" << endl;
    anglesfile << "# columns are: res1 num, res1 type, res2 num, res2 type, "
	       << "ca-ca dist, cb-cb dist, dASA1, dASA2, " << endl;
    anglesfile << "# torsion angle, colinearity angle, angle1, angle2" << endl;
    sys.SomePairAngles(chain1,chain2,cutoff,interface,iresline,nires,anglesfile);
  }

  cout << "Done analyzing interfaces" << endl;

  interfacefile.close();
  cacafile.close();
  cbcbfile.close();
  czczfile.close();
  anglesfile.close();

  return 0;
}

/* provide the next interface specification */
bool getspec(string& interface, string inputfile, string& singlespec, 
	       ifstream& interfacefile) {
  char    buffer[1000];
  string  line;

  if (inputfile == "NULL") {
    if (singlespec == "NULL") {
      cout << "ERROR: cannot get interface specification, none specified" 
	   << endl;
      cout << " Usage: [-fdc] [interface specification (6 characters: 1PDBXX)]"
	   << endl;
      exit(-1);
    } else if (singlespec == "DONE") {
      return(false);
    } else {
      interface = singlespec;
      singlespec = "DONE";
      return(true);
    }
  } else {
    while (!interfacefile.eof()) {
      interfacefile.getline(buffer,1000);
      line = buffer;
      interface = FirstWord(line);
      return(true);
    }
    return(false);
  }
}


void parseflags(int argc, char* argv[], string& inputfile, string& interface, 
		string& interfacedir, string& pdbdir, double& cutoff) {
  inputfile = "NULL";
  interface = "NULL";
  interfacedir = ".";
  pdbdir = ".";
  cutoff = 0.0;

  if( argc >= 2 ) {
    for(int i=1; i<argc; i++) {
      if( *(argv[i]) == '-' ) {
        switch( *(argv[i]+1) ) {
	case 'f' : 
	  inputfile = argv[i+1];
	  i++;
	  break;
	case 'd' : 
	  interfacedir = argv[i+1];
	  i++;
	  break;
	case 'p' : 
	  pdbdir = argv[i+1];
	  i++;
	  break;
	case 'c' :
	  cutoff = ToDbl(argv[i+1]);
	  i++;
	  break;
	default :
	  cout << " Unknown flag " << argv[i] << endl << endl;
	  flaglist(argv[0]);
	}
      } else {
	interface = argv[i];	
	if (interface.length() != 6) {
	  cout << "Inappropriate interface specification '" 
	       << interface << "', length must be 6" << endl;
	  exit(-1);
	}
      }
    }
  }
  if ((inputfile != "NULL")&&(interface != "NULL")) {
    cout << "ERROR: both input filename and single interface was specified, " <<
      "only one can be specified" << endl;
    flaglist(argv[0]);
    exit(-1);
  }

}

void flaglist(char* str) {
  cout << endl 
       << " Usage: " << str
       << " [-fdc] [interface specification (6 characters: 1PDBXX)]"
       << endl;
  cout << " Supported flags: " << endl;
  cout << "  f - ";
  cout << "specify filename from which to read interface specifications\n";
  cout << "  d - ";
  cout << "specify directory from which to read interface residue lists\n";
  cout << "interface file names must have format: [interface_spec]_interface\n";
  cout << "interface residue list file should have format:\n";
  cout << "[number (not read)]  [chain ID] [residue no] [residue] [dASA]\n";
  cout << "Example: 10  A 94 W -90.29\n";
  cout << "  c - ";
  cout << "specify cutoff in Angstroms\n";
  cout << "  p - ";
  cout << "specify directory containing pdb files\n";
  cout << endl;

  exit(1);
}
