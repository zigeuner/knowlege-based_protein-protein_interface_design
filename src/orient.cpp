/*
   This driver reads in a pdb file and attempts to reorient some of the chains
   to interact with a selected region of the remainder of the chains.  It is
   used to provide poses for de novo design of antibodies.
*/

#include <cstdio>
#include <ctime>
#include <unistd.h>
#include <dirent.h>
#include <iostream>
#include <cstdlib>
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

void dumpmovieframe(ofstream&, int, CSystem&, CForcefield&, string&);
void dumpmovieframe(ofstream&, int, CSystem&, double, double, string&);

void parseflags(int, char**, long&, double&, double&, double&,
                double&, int&, int&, int&, int&, string&, string&,
                string&, string&, string&, string&, string&, int&,
                int&, int&, bool&, bool&, bool&, bool&, bool&);

int main (int argc, char *argv[]) { 
  int        n, nequil, ncool, ndecor, nruns, lastiter, nseqs;
  const int  maxseqs=100;
  int        frameinterval, nabchains, maxlines, fftype, resetfreq;
  bool       dumpstats, dumpmovie, caonly, seqopt, usechi, unkref;
  double     rmsd, maxrmsd, potential, temperature, highT, moveT;
  double     minpot = 1e50, binsize;
  long       startseed, seed;
  char       file[50], cstring[80];
  string     ffname, refresidue, commentline, xtra, field[100];
  string     atomtype1, atomtype2, chain1, chain2, line, comment, abchain[100];  
  string     potentialsdir, pdbfile, core, residuesfile, altseq;
  string     abdefn, agdefn, centroids, referencefile;
  string     *outline, header[maxseqs], seq[maxseqs];
  ofstream   output, statsfile, movie;

  /* set the core name of the output files */
  core = "poses";

  /* read run parameters from command line */
  parseflags(argc, argv, startseed, temperature, highT, maxrmsd,
	     binsize, nruns, nequil, ncool, ndecor, pdbfile, abdefn,
	     agdefn, residuesfile, potentialsdir, centroids,
	     referencefile, fftype, frameinterval, resetfreq, caonly,
	     seqopt, dumpmovie, usechi, unkref);

  /* check input */
  nabchains = Split(abdefn,',',abchain,100);
  nseqs = 0;
  if ((nabchains == 1)&&(abdefn.length() > 1)) {  /* read sequence */
    nseqs = ReadFastaFile(abdefn.c_str(),header,seq,maxseqs);
    if (nseqs < 0) {
      cout << "ERROR: could not read sequences from '" << abdefn << "'" << endl;
      exit(-1);
    }
    abdefn = "";
    nabchains = 0;
    /* read headers to get chain ID from "chain X" string */
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

  /* open the XXXXX file for output*/
  sprintf(file,"%s.txt",core.c_str());
  output.open(file, ios::out);
  if(!output) {
    cout << endl << " Error: cannot open file '" << file << 
      "'" << endl;
    exit(-1);
  }

  /* open the tracking file for output*/
  sprintf(file,"runstats.txt");
  statsfile.open(file, ios::out);
  if(!statsfile) {
    cout << endl << " Error: cannot open file '" << file << 
      "'" << endl;
    exit(-1);
  }
  statsfile << "# run statistics file" << endl;
  statsfile << "# initial configuration is from '" << pdbfile << "'" << endl;
  statsfile << "# Antibody is composed of chains: '" << abdefn << "'" << endl;
  statsfile << "# Antigen is composed of chains: '" << agdefn << "'" << endl;
  statsfile << "# number of runs= " << nruns << endl;
  if (highT > 0.0) {
    statsfile << "# Simulated annealing run:" << endl;
    statsfile << "#  " << ndecor 
	      << " steps at decorrelation (high) temperature: " 
	      << highT << "K" << endl;
    statsfile << "#  " << ncool << " steps between decorrelation (high) "
	      << "temperature and equilibrium " 
	      << "temperature (dT=" << highT - temperature << "K)" << endl;
    statsfile << "#  " << nequil << " steps at equilibrium temperature: " 
	      << temperature << "K" << endl;
  } else {
    statsfile << "# Normal MC run at temperature: " << temperature 
	      << "K" << endl;
  }

  /* open an xyz movie file for output*/
  sprintf(file,"movie.xyz");
  movie.open(file, ios::out);
  if(!output) {
    cout << endl << " Error: cannot open file '" << file << 
      "'" << endl;
    exit(-1);
  }

  /* setup recording for landscape information */
  CLandscape landscape(0.0,maxrmsd,int(maxrmsd*(1.0/binsize)));

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

  /* initialize the system moves */
  cout << "Initializing move types... " << endl;
  string movetypes[2] = {"ROTATE","TRANSLATE"};
  string movelist,stepsizes;
  CMove mcmoves(movetypes,2,agdefn,1);
  mcmoves.Disp(movelist,stepsizes);
  statsfile << "# columns are: step, RMSD, system potential energy (kcal/mol), " 
	    << "step sizes for move type(s) (" << movelist << "), seq" << endl;

  /* make the reduced system configuration */
  CSystem  interfacesys;
  string atomtypes[3] = {"CA","CB","MN"};
  //  interfacesys.CopyReduce(fullsys,agdefn,atomtypes,3);
  interfacesys.CopySubset(fullsys,interface,agdefn,atomtypes,3);
  //  interfacesys.CopySubset(fullsys,interface,"",atomtypes,3);

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
      int        ncenters;
      const int  MAXcenters=1000;
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
    CSystem  temp = interfacesys + repulsivesys;
    sys = centroidsys + temp; 
    //    sys = temp + centroidsys; 
  } else {
    //    sys = repulsivesys + interfacesys;
    sys = interfacesys + repulsivesys;
  }
  sys.SetSomeRepulsive(abdefn,repulsive_atomtypes,1);
  sys.WritePDB("starting.pdb");

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

  /* create profit scripts */
  sys.WriteProfitScripts(abdefn,agdefn);

  /* make two copies of the system configuration */
  CSystem syscopy(sys);
  CSystem intermediate(sys);

  /* allocate array for alphabet output */
  maxlines = interface.TotalResidues();
  outline = new string[interface.TotalResidues()];

  //  cout << "DEBUG at " << __FILE__ << ":" << __LINE__ << endl;
  //  forcefield.EvalAlphabet(sys, outline, maxlines);
  //  exit(0);

  /* make initial res-res potential and distance matrices */
  forcefield.Potential(sys,potmtx,commentline);
  potmtx.Write("resresmtx_initial",commentline);
  sys.DistMtx("CA","CA",interface,distmtx,commentline);
  distmtx.Write("distmtx_initial",commentline);

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

  /* write initial line to statistics file */
  rmsd = sys.RMSD(refsys,interface,"CA");
  potential = forcefield.Potential(sys);
  mcmoves.Disp(movelist,stepsizes);
  statsfile << 0 << "  " << Num(rmsd) << "  " << Num(potential) 
	    << "  " << stepsizes; 
  if (seqopt) {
    statsfile << "   " << altseq << endl;
  } else {
    statsfile << endl;
  }

  /* set up run */
  //  frameinterval = 1;  /* set in parseflags */
  xtra = "initial config";
  dumpmovieframe(movie,0,sys,forcefield,xtra);

  /* loop over runs */
  for(int runno=0; runno<nruns; runno++) {  
    int start = runno*(ndecor + ncool + nequil) + 1;

    /* feedback */
    if (nruns > 1) {
      statsfile << "# starting run " << runno << endl;
    }
    cout << "---- Run " << runno << endl;

    if ((resetfreq > 0)&&(runno > 0)&&(!(runno%resetfreq))) {
      cout << "Reseting to initial configuration" << endl;      
      statsfile << "# --LABEL-- " << start-1 
		<< ",-5 " << "\"RESET\"" << endl;
      sys = syscopy;
    }

    /* save a copy of the run start configuration for restarts */
    CSystem runstartconfig(sys);

    /* perform the simulated annealing moves */
    if (highT > 0.0) {
      moveT = highT;
      cout << "Beginning Simulated Annealing moves, simulation length = " 
           << ncool+ndecor << endl;
      cout << " decorrelating at " << highT << "K for " 
    	 << ndecor << " steps..." << endl;
      if (ndecor > 0) {
	statsfile << "# starting " << ndecor << " decorrelation steps" << endl;
	statsfile << "# --VERTICAL_BAR-- " << "col1=" << start-1 << endl;
	statsfile << "# --LABEL-- " << start-1 
		  << ",ymin " << "\"decor\"" << endl;
      }
      int subiter = 0;
      for(int iter=start; iter<ncool+ndecor+start; iter++) {
	subiter += 1;
        dumpstats = false;

	if (subiter == ndecor+1) {
	  cout << " cooling to " << temperature << "K in " 
	       << ncool << " steps..." << endl;
	  statsfile << "# starting " << ncool << " cooling steps" << endl;
	  statsfile << "# --VERTICAL_BAR-- " << "col1=" << iter-1 << endl;
	  statsfile << "# --LABEL-- " << iter-1 
		    << ",ymin " << "\"cool\"" << endl;
	}

        /* change the temperature every 10 steps during cooling run */
        if ((subiter > ndecor+1) && (!(subiter % 10))) {
	  moveT = highT - (iter - start - ndecor)*(highT - temperature)/ncool;
	  //	  cout << iter-start-ndecor << "  " << moveT << endl;
        }

        /* optimize the sequence if required */
        if (seqopt) {
          forcefield.OptimizeSeq(sys);
          altseq = forcefield.AltSeq();
        }

        /* move the system */
	//cout << "---------- SA Iteration " << subiter << " --------" << endl;
        //      mcmoves.Move(sys,forcefield,temperature,true);
        mcmoves.Move(sys,forcefield,moveT,false);

        /* save a snapshot of the system and the random number seed */
	//	sys.WritePDBSnapshot();

        if((frameinterval)&&(!(subiter%frameinterval))) {dumpstats = true;}
        if(subiter <= 100) {dumpstats = true;}

        /* derive information */
        potential = forcefield.CurrentPotential();
	rmsd = sys.RMSD(refsys,interface,"CA");
        landscape.Process(sys,rmsd,potential,40.0,"minnrg");

	/* make sure the molecules haven't separated */
        if (rmsd > maxrmsd) {
	  cout << "WARNING: RMSD > " << Num(maxrmsd) << " Ang, subiter = " 
	       << subiter << ", molecules have flown apart!, " 
	       << endl << "   restarting with initial ";
	  cout << "configuration and current random seed" << endl;
	  seed = sys.GetSeed();
	  sys = runstartconfig;  
	  sys.SetSeed(seed);
	  mcmoves.Reset();
        }

	/* do periodic dumps if necessary */
        if (dumpstats) {
	  // cout << "  rmsd= " << Num(rmsd) << "  iter= " << iter << endl;
	  mcmoves.Disp(movelist,stepsizes);
	  statsfile << iter << "  " << Num(rmsd) << "  " << Num(potential) 
		    << "  " << stepsizes; 
          if (seqopt) {
	    statsfile << "   " << altseq << endl;
          } else {
	    statsfile << endl;
          }

          if (potential < minpot) {
	    sprintf(cstring,"frame %d, RMSD= %s E= %s T= %s",iter,
		    (Num(rmsd)).c_str(),(Num(potential)).c_str(),
		    (Num(moveT)).c_str());
	    commentline = cstring;
	    sys.WriteXYZ("minnrg.xyz",commentline);
          }
        }

        if (dumpmovie) {
	  xtra = "T= " + Num(moveT);
          dumpmovieframe(movie,iter,sys,rmsd,potential,xtra);
        }
	lastiter = iter;
      } /* simulated annealing iteration loop */

      /* save intermediate configuration for further restarts */
      intermediate = sys;  

      /* put a comment in the statsfile at end of cooling steps */
      statsfile << "# finished " << ncool+ndecor << " steps, " 
		<< " beginning equilibration steps" << endl;
      statsfile << "# --VERTICAL_BAR-- " << "col1=" << lastiter << endl;
      statsfile << "# --LABEL-- " << lastiter << ",ymin " << "\"equil\"" << endl;

    } else {
      ncool = 0;
      ndecor = 0;
      intermediate = runstartconfig;  
    }

    /* perform the equilibrium moves on the system */
    cout << "Beginning Equilibrium Monte Carlo moves, simulation length = " 
         << nequil << endl;
    for(int iter=start+ncool+ndecor; iter<start+nequil+ncool+ndecor; iter++) {
      dumpstats = false;

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

      /* save a snapshot of the system and the random number seed */
      //      sys.WritePDBSnapshot();

      if((frameinterval)&&(!((iter-start)%frameinterval))) {dumpstats = true;}
      if(iter-start <= 100) {dumpstats = true;}

      /* derive information */
      potential = forcefield.CurrentPotential();
      rmsd = sys.RMSD(refsys,interface,"CA");
      landscape.Process(sys,rmsd,potential,20.0,"minnrg");

      /* make sure the molecules haven't separated */
      if (rmsd > maxrmsd) {
	cout << "WARNING: RMSD > " << Num(maxrmsd) << " Ang, subiter = " 
	     << iter-start << ", molecules have flown apart!, " << endl;
	cout << "  restarting ";
	if (ncool > 0) {
	  cout << "with initial configuration and current random seed"
	       << endl;	
	} else {
	  cout << "with post-annealing intermediate configuration " 
	       << "and current random seed" << endl;	
	}
	seed = sys.GetSeed();
	sys = intermediate;  
	sys.SetSeed(seed);
	mcmoves.Reset();
      }

      /* do periodic dumps if necessary */
      if (dumpstats) {
        //      cout << "  rmsd= " << Num(rmsd) << "  iter= " << iter << endl;
	mcmoves.Disp(movelist,stepsizes);
	statsfile << iter << "  " << Num(rmsd) << "  " << Num(potential) 
		  << "  " << stepsizes; 
        if (seqopt) {
	  statsfile << "   " << altseq << endl;
        } else {
	  statsfile << endl;
        }

        if (potential < minpot) {
	  sprintf(cstring,"frame %d, RMSD= %s E= %s",iter,
		  (Num(rmsd)).c_str(),(Num(potential)).c_str());
	  commentline = cstring;
	  sys.WriteXYZ("minnrg.xyz",commentline);
        }
      }

      if (dumpmovie) {
        xtra = "equil";
        dumpmovieframe(movie,iter,sys,rmsd,potential,xtra);
      }
      lastiter = iter;
    } /* equilibrium iteration loop */

    /* put a comment in the statsfile at end of equilibrium steps */
    if (runno < nruns-1) {
      statsfile << "# finished " << nequil << " equilibration steps" << endl;
    }

    /* dump cumulative landscape information to file */
    landscape.Write("landscape.txt");

  } /* run number loop */

  /* make final res-res potential and distance matrices */
  potential = forcefield.Potential(sys,potmtx,commentline);
  potmtx.Write("resresmtx_final",commentline);
  rmsd = sys.RMSD(refsys,interface,"CA");
  sys.DistMtx("CA","CA",interface,distmtx,commentline);
  distmtx.Write("distmtx_final",commentline);

  /* dump final structure */
  sprintf(cstring,"frame FINAL, RMSD= %s E= %s",
	  (Num(rmsd)).c_str(),(Num(potential)).c_str());
  commentline = cstring;
  sys.WriteXYZ("final.xyz",commentline);
  sys.WritePDB("final.pdb",commentline);

  movie.close();
  statsfile.close();

  cout << "Done" << endl;

  return 0;
}

void dumpmovieframe(ofstream& xyzfile, int num, CSystem& sys, 
		    CForcefield& ff, string& xtrastring) {
  double     pot;
  char       cstring[80];
  string     comment;

  /* get the system potential for use in header of movie frame */
  pot = ff.Potential(sys);

  /* write a frame to the movie file */
  sprintf(cstring,"frame %d, E= %s",num,(Num(pot)).c_str());
  comment = cstring;
  if (xtrastring.length() > 0) {
    comment += ", " + xtrastring;
  }
  sys.WriteXYZ(xyzfile,comment);
}

void dumpmovieframe(ofstream& xyzfile, int num, CSystem& sys, 
		    double rmsd, double pot, string& xtrastring) {
  char       cstring[80];
  string     comment;

  /* write a frame to the movie file */
  sprintf(cstring,"frame %d, RMSD= %s E= %s",num,
	  (Num(rmsd)).c_str(),(Num(pot)).c_str());
  comment = cstring;
  if (xtrastring.length() > 0) {
    comment += ", " + xtrastring;
  }
  sys.WriteXYZ(xyzfile,comment);
}

void parseflags(int argc, char* argv[], long& startseed, double&
		temperature, double& highT, double& maxrmsd, double&
		binsize, int& nruns, int& nequil, int& ncool, int&
		ndecor, string& inputfile, string& abdefn, string&
		agdefn, string& residuesfile, string& potentialsdir,
		string& centroids, string& referencefile, int& fftype,
		int& iframedump, int& resetfreq, bool& caonly, bool&
		seqopt, bool& movie, bool& usechi, bool& unkref) {

  bool     simann = false, ncool_set = false, ndecor_set = false;

  startseed = 0;
  temperature = 300.0;
  binsize = 0.2;
  highT = 0.0;  /* default, implies no simulated annealing */
  maxrmsd = 50.0;  /* maximum rmsd before reseting to previous config */
  nruns = 1;
  nequil = 1000;
  ncool = 0;
  ndecor = 0;
  iframedump = 20;
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
  movie = false;   /* if true, dump movie frames */
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
	case 'f' :
	  iframedump = ToInt(argv[i+1]);
	  i++;
	  break;
	case 'm' :
	  movie = true;
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
	case 'U' :
	  unkref = false;
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
  cout << "Number of runs: " << nruns << endl;
  if (highT > 0.0) {
    cout << "Simulated annealing run information:" << endl;
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
  cout << "  " << " maximum allowed RMSD from initial structure: " 
	 << maxrmsd << " Angstroms" << endl;
  cout << "  " << " bin size for landscape accumulation: " 
	 << binsize << " Angstroms" << endl;
  if (movie) {
    cout << "  " << " movie frames dumped every " << iframedump 
	 << " steps" << endl;
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
  cout << "use generic LJ potentials with radii from Huang, Love and Mayo \n";
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
  cout << "specify number of runs, default is 1 [-N X]\n";
  cout << "  f - ";
  cout << "frame statistics dump frequency, default is 20 [-f X]\n";
  cout << "  m - ";
  cout << "dump movie frames with same frequency as statistics\n";
  cout << "  R - ";
  cout << "Reset to initial configuration every X runs [-R X]\n";
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
  cout << "  U - ";
  cout << "use Unknown residue type as reference state\n";
  cout << endl;

  exit(1);
}
