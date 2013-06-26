    int     bin;
    double  rawbin,pot;
    std::cout << "Energy routine not yet ready" << std::endl;

    potential = 0.0;
    double r = sepvec.Length();
    if (r > finalsep) {return(true);}
    if (r >= firstsep) {
      rawbin = (r - firstsep)/sepinc;
      bin = int(rawbin);
      if (bin == nentries) {return(table[bin]));
      pot = table[bin] + ((rawbin - bin)/sepinc)*(table[bin+1] - table[bin]);
    }

