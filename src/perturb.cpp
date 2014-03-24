/*
   This driver reads in a pdb file and makes rigid body perturbations to the
   orientation.  It is similar to orient.cpp, but does the following things
   differently:
   1) does not reduce the complexity of the complex, all atoms are retained
   2) outputs a selected number of structures in separate directories
   3) no movie or statistics files are output
   4) only single-temperature MC moves are used, no simulated annealing
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

void parseflags(int, char**, long&, double&, double&, double&, int&,
                int&, string&, string&, string&, string&, string&,
                string&, string&, int&, int&, int&, bool&, bool&,
                bool&, bool&);

int main (int argc, char *argv[]) { 
  int        n, nequil, nbins, lastiter, nseqs;
  const int  maxseqs=100;
  int        frameinterval, nabchains, fftype, resetfreq;
  bool       dumpmovie, caonly, seqopt, usechi;
  double     rmsd, maxrmsd, potential, temperature, binsize;
  long       startseed, seed;
  char       file[50];
  string     ffname, commentline, xtra, field[100];
  string     atomtype1, atomtype2, chain1, chain2, line, comment, abchain[100];  
  string     potentialsdir, pdbfile, core, residuesfile, altseq;
  string     abdefn, agdefn, centroids, referencefile;
  string     header[maxseqs], seq[maxseqs];
  string     refresidue; 
  ofstream   output, statsfile, movie;

  /* read run parameters from command line */
  parseflags(argc, argv, startseed, temperature, maxrmsd, binsize,
	     nbins, nequil, pdbfile, abdefn, agdefn, residuesfile,
	     potentialsdir, centroids, referencefile, fftype,
	     frameinterval, resetfreq, caonly, seqopt, dumpmovie,
	     usechi);

  /* check input */
  nabchains = Split(abdefn,',',abchain,100);
  nseqs = 0;
  if ((nabchains == 1)&&(abdefn.length() > 1)) {
    nseqs = ReadFastaFile(abdefn.c_str(),header,seq,maxseqs);
    if (nseqs < 0) {exit(-1);}
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
    exit(-1);
  }

  /* setup recording for landscape information */
  if (nbins <= 0) {nbins = int(maxrmsd*5);}
  if (binsize != 0.2) {nbins = int(maxrmsd*(1.0/binsize));}
  CLandscape landscape(0.0,maxrmsd,nbins);

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
    //    cout << "retrieved seed= " << inputsys.GetSeed() << endl;
  }
  //  inputsys.Display(1);
  //  exit(0);
  //  inputsys.WritePDB("junk.pdb");
  //  exit(0);

  /* start from a system only containing the antibody and antigen chains */
  CSubset  coresubset = inputsys.GenSubset(agdefn+","+abdefn);
  CSystem  fullsys;
  fullsys.CopySubset(inputsys,coresubset);

  /* define the interface residues */
  CSubset interface(residuesfile,0);
  interface.SetUnknown(abdefn,false);
  //  interface.Display(0);

  /* initialize the forcefield */
  ffname = "KB";
  refresidue = "ALA";  // could also be "UNK"
  if (fftype == 2) {ffname = "SIMPLE";}
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
  if (nabchains == 2) {
    forcefield.SetChainChainInt(abchain[0],abchain[1],false);
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

  /* initialize the system moves */
  cout << "Initializing move types... " << endl;
  string movetypes[2] = {"ROTATE","TRANSLATE"};
  string movelist,stepsizes;
  CMove mcmoves(movetypes,2,agdefn,1);
  mcmoves.Disp(movelist,stepsizes);
  statsfile << "# columns are: step, RMSD, system potential energy (kcal/mol), " 
	    << "step sizes for move type(s) (" << movelist << "), seq" << endl;

  /* make a system of repulsive centroids if any were specified */
  CSystem  centroidsys;
  if (centroids != "NULL") {
    if (centroids == "AGINTERFACE") {/* centroids on all AG interface residues */
      CSubset centroidsubset(residuesfile,agdefn,0);
      centroidsys.CopyCentroids(fullsys,centroidsubset);
    } else if (centroids == "INTERFACE") {/* centroids on interface residues */
      CSubset centroidsubset(residuesfile,0);
      centroidsys.CopyCentroids(fullsys,centroidsubset);
    } else {
      int         ncenters;
      const int   MAXcenters=1000;
      string  centers[MAXcenters], lines[MAXcenters];
      ncenters = Split(centroids,',',centers,MAXcenters);
      for(int i=0; i<ncenters; i++) {  
	lines[i] = Int2String(i+1) + "  " + centers[i].substr(0,1) + "  " 
	  + centers[i].substr(1,centers[i].length()-1) ;
	//      cout << lines[i] << endl;
      }
      CSubset centroidsubset(lines,ncenters,0);
      centroidsys.CopyCentroids(fullsys,centroidsubset);
    }
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
    CSystem  temp = fullsys + repulsivesys;
    sys = centroidsys + temp; 
    //    sys = temp + centroidsys; 
  } else {
    //    sys = repulsivesys + fullsys;
    sys = fullsys + repulsivesys;
  }
  sys.SetSomeRepulsive(abdefn,repulsive_atomtypes,1);
  //  sys.WritePDB("starting.pdb");

  /* read in a reference system configuration if required */
  CSystem refsys;
  if (referencefile != "NULL") {
    sprintf(file,"%s",referencefile.c_str());
    test.open(file, ios::in);
    if(!test) {
      cout << "ERROR: cannot open reference system file '" << file << endl;
      exit(-1);
    } else {
      test.close();
    }
    cout << "reading reference system from file: " << pdbfile << endl;
    refsys.ReadInit(file,false);
  } else {
    refsys = sys;
  }

  /* make a copy of the system configuration */
  CSystem syscopy(sys);

  /* optimize the initial sequence if required */
  if (seqopt) {
    cout << "DEBUG at " << __FILE__ << ":" << __LINE__ << endl;
    forcefield.OptimizeSeq(sys);
    forcefield.CompareSeqs(sys,false);
    altseq = forcefield.AltSeq();
    cout << altseq << endl;
    cout << "DEBUG at " << __FILE__ << ":" << __LINE__ << endl;
    //    exit(0);
  }

  /* perform the equilibrium moves on the system */
  cout << "Beginning Equilibrium Monte Carlo moves, simulation length = " 
       << nequil << endl;
  for(int iter=0; iter<nequil; iter++) {
    
    /* optimize the sequence if required */
    if (seqopt) {
      //      cout << "DEBUG at " << __FILE__ << ":" << __LINE__ << endl;
      forcefield.OptimizeSeq(sys);
      //      forcefield.CompareSeqs(sys,false);
      altseq = forcefield.AltSeq();
      //      cout << "DEBUG at " << __FILE__ << ":" << __LINE__ << endl;
      //      exit(0);
    }

    /* move the system */
    //    mcmoves.Move(sys,forcefield,temperature,true);
    mcmoves.Move(sys,forcefield,temperature,false);

    /* derive information */
    potential = forcefield.CurrentPotential();
    rmsd = sys.RMSD(refsys,interface,"CA");
    landscape.Process(sys,rmsd,potential,20.0,"minnrg");

    /* make sure the molecules haven't separated */
    if (rmsd > maxrmsd) {
      cout << "WARNING: RMSD > " << Num(maxrmsd) << " Ang, iter = " 
	   << iter << ", molecules have flown apart!, "
	   << "restarting with initial configuration" << endl;
      seed = sys.GetSeed();
      sys = syscopy;  
      sys.SetSeed(seed);
      mcmoves.Reset();
    }

    lastiter = iter;
  } /* equilibrium iteration loop */

  /* dump cumulative landscape information to file */
  landscape.Write("landscape.txt");

  cout << "Done" << endl;

  return 0;
}

void parseflags(int argc, char* argv[], long& startseed, double&
		temperature, double& maxrmsd, double& binsize, int&
		nbins, int& nequil, string& inputfile, string& abdefn,
		string& agdefn, string& residuesfile, string&
		potentialsdir, string& centroids, string&
		referencefile, int& fftype, int& iframedump, int&
		resetfreq, bool& caonly, bool& seqopt, bool& movie,
		bool& usechi) {

  startseed = 0;
  temperature = 300.0;
  binsize = 0.2;
  maxrmsd = 2.0;  /* maximum rmsd before reseting to previous config */
  nbins = 0;      /* implies default behavior */
  nequil = 1000;
  resetfreq = 0;
  inputfile = "NULL";
  abdefn = "NULL";
  agdefn = "NULL";
  centroids = "NULL";
  residuesfile = "NULL";
  referencefile = "NULL";
  potentialsdir = ".";
  fftype = 0; /* 0= knowledge-based FF with normal residue types */
              /* 1= knowledge-based FF with unknown type Ab residues */
              /* 2= LJ potentials */
  caonly = false;  /* if true, use only CA-CA pairs (no CB-CB) */
  seqopt = false;  /* if true, optimize the unknown residue types at each step */
  usechi = false;  /* if true, add chi-derived res-res potentials */

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
	case 'N' :
	  nbins = ToInt(argv[i+1]);
	  i++;
	  break;
	case 'R' :
	  resetfreq = ToInt(argv[i+1]);
	  i++;
	  break;
	case 'M' :
	  maxrmsd = ToDbl(argv[i+1]);
	  i++;
	  break;
	case 'x' :
	  usechi = true;
	  break;
	case 'C' :
	  centroids = argv[i+1];
	  i++;
	  break;
	case 'w' :
	  referencefile = argv[i+1];
	  i++;
	  break;
	case 'B' :
	  binsize = ToDbl(argv[i+1]);
	  i++;
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

  /* Some feedback for the user */
  cout << "Number of bins: " << nbins << endl;
  cout << "Normal MC run at temperature: " << temperature << "K" << endl;
  cout << "  " << " maximum allowed RMSD from initial structure: " 
	 << maxrmsd << " Angstroms" << endl;
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
  cout << "use generic LJ potentials with radii from Huang, Love and Mayo \n";
  cout << "  a - ";
  cout << "use only Calpha-Calpha interaction type\n";
  cout << "  o - ";
  cout << "optimize unknown residue types at each step (also implies -u flag)\n";
  cout << "  N - ";
  cout << "specify number of bins, default is maxRMSD/0.2 [-N X]\n";
  cout << "  M - ";
  cout << "specify Maximum allowed RMSD vs initial config (default 50) [-M X]\n";
  cout << "  x - ";
  cout << "evaluate chi-derived res-res potentials\n";
  cout << "  C - ";
  cout << "place repulsive sidechain Centroids [Ex: -C L45,A343]\n";
  cout << "  w - ";
  cout << "specify Wildtype system as reference for RMSD's [-w filename]\n";
  cout << "  B - ";
  cout << "specify bin size for landscape structure dumps [-B X]\n";
  cout << endl;

  exit(1);
}
