/*
   This driver reads in a structure file and evaluates the single-point energy
*/

#include <cstdio>
#include <ctime>
#include <unistd.h>
#include <dirent.h>
#include <iostream>
#include "include/system.h"
#include "include/forcefield.h"
#include "include/move.h"

using std::cout;
using std::endl;
using std::ios;
using std::string;
using std::ofstream;
using std::ifstream;

void parseflags(int, char**, long&, double&, double&, int&, int&,
                int&, int&, string&, string&, string&, string&,
                string&, string&, int&, bool&, bool&, bool&, bool&);

int main (int argc, char *argv[]) { 
  int        nequil, ncool, ndecor, nruns, maxseqs=100, nseqs;
  int        n, nabchains, fftype;
  bool       caonly, seqopt, usechi, unkref;
  double     potential, temperature, highT, pot1; 
  long       startseed;
  char       file[50];
  string     ffname, refresidue, commentline, xtra, field[100];
  string     atomtype1, atomtype2, chain1, chain2, line, comment, abchain[100];  
  string     potentialsdir, pdbfile, core, residuesfile, altseq;
  string     abdefn, agdefn, centroids;
  string     header[maxseqs], seq[maxseqs];
  ofstream   output, statsfile, movie;

  /* set the core name of the output files */
  core = "poses";

  /* read run parameters from command line */
  parseflags(argc, argv, startseed, temperature, highT, nruns, nequil,
	     ncool, ndecor, pdbfile, abdefn, agdefn, residuesfile,
	     potentialsdir, centroids, fftype, caonly, seqopt, usechi,
	     unkref);

  /* check input */
  nabchains = Split(abdefn,',',abchain,100);
  nseqs = 0;
  if ((nabchains == 1)&&(abdefn.length() > 1)) {
    nseqs = ReadFastaFile(abdefn.c_str(),header,seq,maxseqs);
    abdefn = "";
    nabchains = 0;
    for(int i=0; i<nseqs; i++) {
      n = Split(header[i],field,100);
      for(int j=0; j<n; j++) {
	if (field[j] == "chain") {abchain[i] = field[j+1]; nabchains +=1; break;}
      }
      //      cout << header[i] << endl << seq[i] << endl;
    }
    abdefn = Join(abchain,',',nabchains);
  }
  if (nabchains > 2) {
    cout << "WARNING: what's up, there's more than 2 antibody chains?" << endl;
    //    exit(-1);
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
  CSystem inputsys(file,0);
  cout << "Setting system temperature to: " << Num(temperature) << endl;
  inputsys.SetT(temperature);
  if (startseed != 0) {
    cout << "Setting random number seed to: " << startseed << endl;
    inputsys.SetSeed(startseed);
  }
  //  inputsys.Display(1);
  //  inputsys.WritePDB("junk.pdb");
  //  exit(0);

  /* start from a system only containing the antibody and antigen chains */
  CSubset  coresubset = inputsys.GenSubset(agdefn+","+abdefn);
  //  coresubset.Display(0);
  CSystem  fullsys;
  fullsys.CopySubset(inputsys,coresubset);

  /* define the interface residues */
  CSubset interface(residuesfile,agdefn+","+abdefn,0);
  interface.SetUnknown(abdefn,false);
  //  interface.Display(0);
  CMatrix potmtx(interface.TotalResidues());
  CMatrix distmtx(interface.TotalResidues());

  /* initialize the forcefield */
  ffname = "KB";
  if (fftype == 2) {ffname = "SIMPLE";}
  refresidue = "ALA";
  if (unkref) {refresidue = "UNK";}
  if (caonly) {ffname += "_CA_ONLY";}
  cout << "Initializing forcefield of type '" << ffname << "' ... " << endl;
  CForcefield forcefield(ffname,potentialsdir,refresidue,1);
  if (fftype == 2) {
    forcefield.SetUnknownInt(false);
  } else if (fftype == 1) {
    forcefield.SetUnknownInt(true);
    cout << "Treating chain(s) ";
    for(int i=0; i<nabchains; i++) {
      cout << abchain[i];
      if (i != nabchains-1) {cout << ",";}
    }
    cout << " as possessing unknown residue types" << endl;
  } else {
    forcefield.SetUnknownInt(false);
  }
  forcefield.SetSubset(interface);
  if (nabchains >= 2) {
    for(int i=0; i<nabchains; i++) {    
      for(int j=i+1; j<nabchains; j++) {    
	forcefield.SetChainChainInt(abchain[i],abchain[j],false);
      }
    }
  }
  if (seqopt) {
    cout << "Optimizing unknown residue types for each configuration" << endl;
    forcefield.UseSubsetSeq(true);
  }
  if (nseqs > 0) {
    forcefield.UseSubsetSeq(true);
    //    cout << "chain sequences: " << forcefield.AltSeq() << endl;
    for(int i=0; i<nabchains; i++) {
      cout << "setting chain '" << abchain[i] << "' to seq: " << seq[i] << endl;
      forcefield.SetChainSeq(fullsys,abchain[i],seq[i]);
    }
    //    cout << "chain sequences: " << forcefield.AltSeq() << endl;
  }
  forcefield.EvalChiPot(usechi);
  forcefield.SetChiDistCutoff(12.0);
  //  forcefield.Display(0); 
  //  exit(0);

  /* make the reduced system configuration */
  CSystem  interfacesys;
  string atomtypes[3] = {"CA","CB","MN"};
  //  interfacesys.CopyReduce(fullsys,agdefn,atomtypes,3);
  interfacesys.CopySubset(fullsys,interface,agdefn,atomtypes,3);
  //  interfacesys.CopySubset(fullsys,interface,"",atomtypes,3);
  interfacesys.WritePDB("starting.pdb");

  /* make a system of repulsive centroids if any were specified */
  CSystem  centroidsys;
  if (centroids != "NULL") {
    int     ncenters,MAXcenters=1000;
    string  centers[MAXcenters], lines[MAXcenters];
    ncenters = Split(centroids,',',centers,MAXcenters);
    for(int i=0; i<ncenters; i++) {  
      lines[i] = Int2String(i+1) + "  " + centers[i].substr(0,1) + "  " 
	+ centers[i].substr(1,centers[i].length()-1) ;
      //      cout << lines[i] << endl;
    }
    CSubset centroidsubset(lines,ncenters,0);
    centroidsys.CopyCentroids(fullsys,centroidsubset);
    centroidsys.SetAllRepulsive();
  }
  if (centroidsys.NAtoms() == 0) {
    cout << "Adding no additional repulsive centroids to system" << endl;
  } else {
    cout << "Adding " << centroidsys.NAtoms() 
	 << " repulsive sidechain centroids to system" << endl;    
  }

  /* make the normal repulsive centers (first as separate system) */
  /* want to add all CA's not on interface on both side AND */
  /* all CA's on antibody side */
  CSystem  repulsivesys;
  CSubset fullsubset = fullsys.GenSubset();
  CSubset repulsivesubset = fullsubset - interface;
  string repulsive_atomtypes[1] = {"CA"};
  repulsivesys.CopySubset(fullsys,repulsivesubset,agdefn+","+abdefn,
			  repulsive_atomtypes,1);
  if (repulsivesys.NAtoms() == 0) {
    cout << "Added no additional repulsive centers to system" << endl;
  } else {
    cout << "Added " << repulsivesys.NAtoms() 
	 << " repulsive centers to system at " << repulsive_atomtypes[0]
	 << " atomic centers not on interface" << endl;    
  }
  repulsivesys.SetAllRepulsive();

  /* make the new merged system configuration */
  CSystem  sys;
  if (centroids != "NULL") {
    CSystem  temp = interfacesys + repulsivesys;
    sys = centroidsys + temp; 
    //    sys = temp + centroidsys; 
  } else {
    //    sys = repulsivesys + interfacesys;
    sys = interfacesys + repulsivesys;
  }
  sys.SetSomeRepulsive(abdefn,repulsive_atomtypes,1);

  /* make initial res-res potential and distance matrices */
  pot1 = forcefield.Potential(sys,potmtx,commentline);
  potmtx.Write("resresmtx",commentline);
  sys.DistMtx("CA","CA",interface,distmtx,commentline);
  distmtx.Write("distmtx",commentline);

  //  forcefield.Display(0);

  potential = forcefield.Potential(sys);
  cout << "full system potential= " << potential << endl;
  if (fabs(pot1 - potential) > 1.0e-5) {
    cout << "WARNING: results from potential matrix and normal potential "
	 << "routines do not match " << endl;
    cout << " potential matrix routine nrg = " << Num(pot1) << endl;
    cout << " normal potential routine nrg = " << Num(potential) << endl;
  }
  cout << "Done" << endl;

  return 0;
}

void parseflags(int argc, char* argv[], long& startseed, double&
		temperature, double& highT, int& nruns, int& nequil,
		int& ncool, int& ndecor, string& inputfile, string&
		abdefn, string& agdefn, string& residuesfile, string&
		potentialsdir, string& centroids, int& fftype, bool&
		caonly, bool& seqopt, bool& usechi, bool& unkref) {
  bool     simann = false, ncool_set = false, ndecor_set = false;

  startseed = 0;
  temperature = 300.0;
  highT = 0.0;  /* default, implies no simulated annealing */
  nruns = 1;
  nequil = 1000;
  ncool = 0;
  ndecor = 0;
  inputfile = "NULL";
  abdefn = "NULL";
  agdefn = "NULL";
  centroids = "NULL";
  residuesfile = "NULL";
  potentialsdir = ".";
  fftype = 0; /* 0= knowledge-based FF with normal residue types */
              /* 1= knowledge-based FF with unknown type Ab residues */
              /* 2= LJ potentials */
  caonly = false;  /* if true, use only CA-CA pairs (no CB-CB) */
  seqopt = false;  /* if true, optimize the unknown residue types at each step */
  usechi = false;  /* if true, add chi-derived res-res potentials */
  unkref = false;  /* if true, use unknown residue type as reference state */

  if( argc >= 2 ) {
    for(int i=1; i<argc; i++) {
      if( *(argv[i]) == '-' ) {
        switch( *(argv[i]+1) ) {
	case 'D' : 
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
	  nequil = ToInt(argv[i+1]);
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
	case 'o' :
	  seqopt = true;
	  fftype = 1;
	  break;
	case 'h' :
	  highT = ToDbl(argv[i+1]);
	  simann = true;
	  i++;
	  break;
	case 'c' :
	  ncool = ToInt(argv[i+1]);
	  ncool_set = true;
	  i++;
	  break;
	case 'd' :
	  ndecor = ToInt(argv[i+1]);
	  ndecor_set = true;
	  i++;
	  break;
	case 'N' :
	  nruns = ToInt(argv[i+1]);
	  i++;
	  break;
	case 'x' :
	  usechi = true;
	  break;
	case 'C' :
	  centroids = argv[i+1];
	  i++;
	  break;
	case 'U' :
	  unkref = true;
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
  if ((inputfile == "NULL")||(abdefn == "NULL")||(agdefn == "NULL")) {
    cout << "ERROR: must specify Ab defn, Ag defn and input filename" << endl;
    flaglist(argv[0]);
    exit(-1);
  }

  /* set ncool and ndecor based on nequil if highT was specified */
  //  if ((simann)&&(!ncool_set)) {ncool = int(nequil/10);}
  //  if ((simann)&&(!ndecor_set)) {ndecor = int(nequil/10);}
  if ((simann)&&(!ncool_set)) {ncool = nequil;}  
  if ((simann)&&(!ndecor_set)) {ndecor = nequil;}

  /* Some feedback for the user */
  if (highT > 0.0) {
    cout << "Simulated annealing run:" << endl;
    cout << "  " << ndecor << " steps at decorrelation (high) temperature: " 
	 << highT << "K" << endl;
    cout << "  " << ncool 
	 << " steps between decorrelation (high) temperature and equilibrium " 
	 << "temperature (dT=" << highT - temperature << "K)" << endl;
    cout << "  " << nequil << " steps at equilibrium temperature: " 
	 << temperature << "K" << endl;
  } else {
    cout << "Normal MC run at temperature: " << temperature << "K" << endl;
  }
  cout << "Antibody is composed of chains: '" << abdefn << "'" << endl;
  cout << "Antigen is composed of chains: '" << agdefn << "'" << endl;

}

void flaglist(char* str) {
  cout << endl 
       << " Usage: " << str
       << " -b [chain ID's] -g [chain ID's] [-fdrsnTulao] "
       <<" [system .pdb or .xyz filename]"
       << endl;
  cout << " Supported flags: " << endl;
  cout << "  D - ";
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
  cout << "specify number of MC steps at the run temperatue [-n X]\n";
  cout << "  T - ";
  cout << "specify starting temperature in Kelvin\n";
  cout << "  u - ";
  cout << "use unknown residue type for antibody residues\n";
  cout << "  l - ";
  cout << "use generic LJ potentials with radi from Huang, Love and Mayo \n";
  cout << "  a - ";
  cout << "use only Calpha-Calpha interaction type\n";
  cout << "  o - ";
  cout << "optimize unknown residue types at each step (also implies -u flag)\n";
  cout << "  h - ";
  cout << "specify high temperature and use simulated annealing (-h X)\n";
  cout << "  c - ";
  cout << "specify number of cooling steps in the run [-c X]\n";
  cout << "  d - ";
  cout << "specify number of high T decorrelation steps in the run [-D X]\n";
  cout << "  N - ";
  cout << "specify number runs, default is 1 [-N X]\n";
  cout << "  x - ";
  cout << "evaluate chi-derived res-res potentials\n";
  cout << "  C - ";
  cout << "place repulsive sidechain Centroids [Ex: -C L45,A343]\n";
  cout << "  U - ";
  cout << "use Unknown residue type as reference state\n";
  cout << endl;

  exit(1);
}
