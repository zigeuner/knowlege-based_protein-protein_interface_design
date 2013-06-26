// These are general routines for doing string manipulation.

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "include/stringutils.h"
#include "include/config.h"

using std::ifstream;
using std::string;
using std::cout;
using std::endl;
using std::ios;

/* Trim leading and trailing spaces from a string */
void Trim(string& str) {
  /* trim the string to drop leading spaces */
  while (str[0] == ' ') {
    if (str.length() == 0) {return;}  // RH7 seems to need this
    str.erase(0,1);
  }
  /* trim the string to drop trailing spaces */
  while (str[str.length() - 1] == ' ') {
    if (str.length() == 0) {return;}  // RH7 seems to need this
    str.erase(str.length() - 1,1);
  }
}

/* Extract the first word from a string */
string FirstWord(string str) {
  int       firstblank, i;
  string    newstr = str;
  
  Trim(newstr);
  firstblank = newstr.length();
  for(i=0;i<firstblank;i++) {
    if (newstr.substr(i,1) == " ") {firstblank = i; break;}
  }
  newstr.erase(i,(newstr.length() - i));

  return(newstr);
}

/* Extract the file ending (ie .xyz, .pdb etc) from a string */
string FileEnding(string str) {
  int       lastdot, len;
  string    newstr = str;
  
  Trim(newstr);
  len = newstr.length();
  lastdot = -1;
  for(int i=len;i>=0;i--) {
    if (newstr.substr(i,1) == ".") {lastdot = i; break;}
  }

  if (lastdot == -1) {newstr = ""; return(newstr);}

  newstr.erase(0,lastdot+1);

  return(newstr);
}

/* Convert a string to a integer, uses ToDbl and a type cast */
int ToInt(string str) {  
  double num = ToDbl(str);
  return((int) num);
}

/* Convert a string to a double (doesn't recognize scientific notation yet */
double ToDbl(string str) {
  int            pos,pos1,pos2,pos3,power;
  unsigned long  i,digit;
  double         output = 0.0e0;
  bool           negative = false;
  string         trimmed, nums = "0123456789+-.";

  /* pos1 = position immediately before decimal point or end of string 
     pos2 = position immediately after decimal point
     pos3 = end of non-space portion of string */

  //  cout << "input string to ToDbl: '" << str << "'" << endl;

  /* trim the string to drop leading and trailing spaces */
  Trim(str);

  /* trim the string to drop anything after an intermediate space */
  i = str.find(" ",0);
  if (i == string::npos) {i = str.length();}
  //  if (i < 0) {i = str.length();}
  trimmed = str.substr(0,i);
  pos3 = trimmed.length() - 1;

  //  cout << "trimmed string: '" << trimmed << "'" << endl;
  //  cout << "position 3: " << pos3 << endl;

  /* use decimal place, if there is one, to determine other positions */
  pos1 = pos3;
  pos2 = pos3 + 1;
  i = str.find('.',0);
  if (i != string::npos) {
    //  if (i >= 0) {
    pos1 = i - 1;
    pos2 = i + 1;
  }

  //  cout << "position 1: " << pos1 << endl;
  //  cout << "position 2: " << pos2 << endl;

  /* process digits before the decimal point */
  power = 0;
  for(pos=pos1;pos>=0;pos--) {
    digit = nums.find(trimmed[pos],0);
    if ((digit == string::npos)||(digit > 12)) {
      //    if ((digit < 0)||(digit > 12)) {
      cout << " ERROR on line: " << __LINE__ << " of " << __FILE__ << endl;
      cout << " Encountered an unexpected non-numeric character '" << 
	trimmed[pos] << "' in string '" << str << "'" << endl;
      exit(-1);
    }
    if (digit == 11) {
      negative = true;
    }
    if (digit < 11) {  // normal digit
      output += digit*pow(10,power);
      power += 1;
    }
  }
  if (negative) {output *= -1;}

  //  cout << "number before decimal place: " << output << endl;

  /* process digits after the decimal point */
  power = -1;
  for(pos=pos2;pos<=pos3;pos++) {
    digit = nums.find(trimmed[pos],0);
    if ((digit == string::npos)||(digit > 12)) {
      //    if ((digit < 0)||(digit > 12)) {
      cout << " ERROR on line: " << __LINE__ << " of " << __FILE__ << endl;
      cout << " Encountered an unexpected non-numeric character '" << 
	trimmed[pos] << "' in string '" << str << "'" << endl;
      exit(-1);
    }
    if (negative) {
      output -= digit*pow(10,power);
    } else {
      output += digit*pow(10,power);
    }
    power -= 1;
  }

  //  cout << "output double: " << output << endl;

  return(output);
}

/* Convert a string to a integer, ignoring non-digit characters */
int FilterToInt(string str) {  
  string         shortstr = "";
  string         nums = "0123456789";
  int            len = str.length();
  unsigned long  digit;

  for(int i=0;i<len;i++) {
    digit = nums.find(str[i],0);
    if (digit == string::npos) {continue;}
    shortstr += Int2String(digit);
  }
  return(ToInt(shortstr));
}

/* Convert a 3-letter amino acid type to a 1-letter symbol */
/* Returns a blank string if it couldn't match the input */
string AAType2Symbol(string str) {
  int     i;
  string  symbol = "z";

  for(i=0;i<NAATypes;i++) {
    //    cout << i << "  " << AAType[i] << " =? " << str << endl;
    if (AAType[i] == str) {
      symbol = AASymbol[i];
      break;
    }
  }
  
  return(symbol);
}

/* Convert a 1-letter amino acid symbol to an integer index */
/* The index is the array index in the AAIntType array (in config.h) */
/* Returns -1 if it couldn't match the input */
int AAIntType2Index(string str) {
  int     i,index;

  index = -1;
  for(i=0;i<NIntResTypes;i++) {
    if (AAIntType[i] == str) {
      index = i;
      break;
    }
  }
  
  return(index);
}

/* Convert a 3-letter amino acid code to an integer index */
/* The index is the array index in the AAIntType array (in config.h) */
/* Returns -1 if it couldn't match the input */
int AAType2IntIndex(string str) {
  string symbol = AAType2Symbol(str); 
  return(AAIntType2Index(symbol));
}

/* Convert a double to a string using appropriate notation */
string Num(const double number) {
  char      cstring[80];
  string    line;

  if (number == 0.0) {
    sprintf(cstring,"%.2f",number);
  } else if (fabs(number) < 1.0e-2) {
    sprintf(cstring,"%.3e",number);
  } else if (fabs(number) > 1000) {
    sprintf(cstring,"%.3e",number);
  } else {
    sprintf(cstring,"%.2f",number);
  }
  line = cstring;
  return(line);
}

/* Convert a double to a string using appropriate notation */
string Num(const int number) {
  char      cstring[80];
  string    line;

  sprintf(cstring,"%d",number);
  line = cstring;
  return(line);
}

/* Convert an integer to a string using appropriate notation */
string Int2String(const long int number) {
  char      cstring[80];
  string    line;

  sprintf(cstring,"%ld",number);
  line = cstring;
  return(line);
}
string Int2String(const long unsigned int number) {
  char      cstring[80];
  string    line;

  sprintf(cstring,"%ld",number);
  line = cstring;
  return(line);
}
string Int2String(const int number) {
  char      cstring[80];
  string    line;

  sprintf(cstring,"%d",number);
  line = cstring;
  return(line);
}

/* Convert a double to a string using appropriate notation */

/* Convert a 1-letter element type to a mass */
double GetMass(string str) {
  int     i;

  for(i=0;i<NElements;i++) {
    if (Element[i] == str) {
      return(ElementMass[i]);
    }
  }
  
  return(-1.0);
}

/* Resize an array of strings */
void Resize(string* str, const int oldsize, const int newsize) {
  string*     temp;

  /* copy old array */
  temp = new string[oldsize];
  for(int i=0;i<oldsize;i++) {
    temp[i] = str[i];
  }
 
  /* deallocate and allocate new */
  delete[] str;

  /* copy old array into new array */
  str = new string[newsize];
  for(int i=0;i<oldsize;i++) {
    str[i] = temp[i];
  }
}

/* Read contents of file into string array */
/* Return -1 if it wasn't possible to open file */
/* NOTE: array must already be sized!! */
int ReadFile(const char* rfile, string* strline, int arraysize) {
  int      lineno = 0;
  char     buffer[1000];
  string   line;
  ifstream file;

  /* Open the file */
  file.open(rfile, ios::in);
  if(!file) {
    cout << " ERROR: ReadFile cannot open file '" << rfile << "'" << endl;
    return(-1);
  }

  /* Read file */
  while (!file.eof()) {
    /* read the line */
    file.getline(buffer,1000);
    line = buffer;

    //    cout << "'" << line << "'" << endl;
    StripComment(line);
    if (line.length() == 0) {continue;}
    //    cout << "'" << line << "'" << endl;

    strline[lineno] = line;
    lineno += 1;

    if (lineno >= arraysize+1) {
      cout << "ERROR: passed array was too small, returning incomplete file\n";
      return(lineno);
    }
  }  

  file.close();

  return(lineno);
}

/* Read contents of FASTA-formatted file into two string arrays */
/* Return -1 if it wasn't possible to open file */
/* NOTE: arrays must already be sized!! */
int ReadFastaFile(const char* rfile, string* header, string* seq, 
		  int arraysize) {
  int      nseqs;
  char     buffer[1000];
  string   line;
  ifstream file;

  /* Open the file */
  file.open(rfile, ios::in);
  if(!file) {
    cout << " ERROR: ReadFile cannot open file '" << rfile << "'" << endl;
    return(-1);
  }

  /* Read file */
  nseqs = -1;
  while (!file.eof()) {
    /* read the line */
    file.getline(buffer,1000);
    line = buffer;

    //    cout << "'" << line << "'" << endl;
    StripComment(line);
    if (line.length() == 0) {continue;}
    //    cout << "'" << line << "'" << endl;

    /* if the line begins with a '>', then it's a header */
    if (line.substr(0,1) == ">") {
      line.erase(0,1);
      Trim(line);
      nseqs += 1;
      if (nseqs >= arraysize) {
	cout << "ERROR: passed arrays were too small, incomplete sequence!\n";
	file.close();
	return(nseqs);
      }
      header[nseqs] = line;
      continue;
    }

    seq[nseqs] += line;
  }  

  file.close();

  return(nseqs+1);
}

/* Read contents of file and determines the number of non-blank */
/* lines after they have been stripped of comments */
/* Returns the number of non-blank (after comments stripped) lines */
int NLines(const char* rfile) {
  int      lineno = 0;
  char     buffer[1000];
  string   line;
  ifstream file;

  /* Open the file */
  file.open(rfile, ios::in);
  if(!file) {
    cout << " ERROR: NLines cannot open file '" << rfile << "'" << endl;
    return(-1);
  }

  /* Read file */
  while (!file.eof()) {
    file.getline(buffer,1000);
    line = buffer;
    StripComment(line);
    if (line.length() == 0) {continue;}
    lineno += 1;
  }  

  file.close();

  return(lineno);
}

/* Read contents of matrix file into a 2D string array */
/* Returns false if it wasn't possible to open file */
/* NOTE: 2D array must already be sized!! */
bool ReadMtxFile(const char* rfile, string**& mtx, int arraysize, 
		 string& lastcomment) {
  int      lineno = 0, reallineno = 1, nfields;
  char     buffer[1000];
  string   templine,field[1000];
  ifstream file;

  /* Size and initialize the matrix */
  mtx = new string*[arraysize];
  for(int i=0; i<arraysize; i++) {
    mtx[i] = new string[arraysize];
    for(int j=0; j<arraysize; j++) {
      mtx[i][j] = "";
    }
  }

  /* Open the file */
  file.open(rfile, ios::in);
  if(!file) {
    cout << " ERROR: ReadMtxFile cannot open file '" << rfile << "'" << endl;
    return(false);
  }

  /* Read file */
  while (!file.eof()) {
    /* read the line */
    file.getline(buffer,1000);
    templine = buffer;

    //    cout << "'" << templine << "'" << endl;
    if (templine[0] == '#') {
      lastcomment = templine;
      //      cout << "stored comment line: '" << lastcomment << "'" << endl;
    }
    StripComment(templine);
    if (templine.length() == 0) {continue;}
    //    cout << "'" << templine << "'" << endl;

    /* split the line and store in matrix */
    nfields = Split(templine,field,1000);
    if (nfields != arraysize) {
      cout << "ERROR: ReadMtxFile read the wrong number of entries " 
	   << "on line " << reallineno << " of file '" << rfile << "'" << endl;
      cout << "  expected " << arraysize << " but read " << nfields << endl;
      exit(-1);
      return(false);
    }
    for(int i=0; i<nfields; i++) {
      mtx[lineno][i] = field[i];
    }

    lineno += 1;
    reallineno += 1;

    if (lineno >= arraysize+1) {
      cout << "ERROR: somehow ReadMtxFile screwed-up file length\n";
      exit(-1);
      return(false);
    }
  }  

  file.close();

  return(true);
}

/* Split contents of a string into substrings using the given character */
/* as separators.  Returns the number of substrings */
int Split(string input, const char separator, string* str, int arraysize) {
  bool          lastsep = true;
  int           nsub = 0;
  unsigned int  length = input.length();
  string        letter, sep;

  if (length == 0) {return(0);}
  sep = separator;
  str[0].clear();
  for(unsigned int i=0; i<length; i++) {
    letter = input.substr(i,1);
    if (letter == sep) {
      if (!lastsep) {nsub += 1; str[nsub].clear();}
      lastsep = true;
    } else {
      str[nsub] += letter;
      lastsep = false;
    }
  }
  nsub += 1;

  return(nsub);
}

/* Split contents of a string into substrings using space separator */
/* Returns the number of substrings */
int Split(string input, string* str, int arraysize) {
  return(Split(input,' ',str,arraysize));
}

/* Join the contents of an array using the separator */
string Join(string* input, const char separator, int nitems) {
  string   output = "";
  for(int i=0; i<nitems; i++) {
    output += input[i];
    if (i != nitems-1) {
    output += separator;
    }
  }
  return(output);
}

/* Join the contents of an array using a default space */
string Join(string* input, int nitems) {
  return(Join(input,' ',nitems));
}

/* Return the index of a given string in an array */
/* Returns -1 if none of the entries matches */
int IndexInArray(string input, string* str, int arraysize) {
  for(int i=0;i<arraysize;i++) {
    if (input == str[i]) {return(i);}
  }  
  return(-1);
}


/*------------------------------------------------------------------------*/
/* Strips the comments from a string class variable.  It simply removes */
/* any characters from the string after a '#' character. */
/* Requires:  str -- character string to strip */
/*------------------------------------------------------------------------*/
void StripComment(string& str) {
  for(unsigned int i=0; i<str.length(); i++) {
    if (str[i] == '#') {
      str.erase(i,str.length() - i);
    }
  }
}


/*------------------------------------------------------------------------*/
/* Strips the comments from a character string.  For now, this means that */
/* it simply removes any characters from the string after a '#' character. */
/* Returns the number of non-space characters before a '#' was encountered, */
/* so if it returns 0, then the whole line is a comment */
/* Requires:  original -- character string to strip */
/*            new -- resultant character string */
/*------------------------------------------------------------------------*/
/*int stripcomment(char *original, char *newstr) {
  int  i,length,nchars = 0;

  length = (int) strlen(original);
  for(i=0; i<length; i++) {
    if (original[i] == '#') {
      newstr[i] = '\0';
      return(nchars);
    } else {
      newstr[i] = original[i];
      newstr[i+1] = '\0';
      if (original[i] != ' ') nchars += 1;
    }
  }

//  cout << "number of characters = " << nchars << endl;
  return(nchars);

}
*/
