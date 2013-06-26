// These are the prototype statements for the numerical utilities

#ifndef __STRINGUTILS_H
#define __STRINGUTILS_H

#include<string>

void Trim(std::string&);
std::string FirstWord(std::string);
std::string FileEnding(std::string);
int ToInt(std::string);
int FilterToInt(std::string);
double ToDbl(std::string);
double GetMass(std::string);
std::string AAType2Symbol(std::string);
int AAIntType2Index(std::string);
int AAType2IntIndex(std::string);
std::string Num(const double);
std::string Num(const int);
std::string Int2String(const long int);
std::string Int2String(const long unsigned int);
std::string Int2String(const int);
void Resize(std::string*, const int, const int);
int ReadFile(const char*, std::string*, int);
int ReadFastaFile(const char*, std::string*, std::string* seq, int);
int NLines(const char*);
bool ReadMtxFile(const char*, std::string**&, int, std::string&);
int IndexInArray(std::string, std::string*, int arraysize);
int Split(std::string, const char, std::string*, int);
int Split(std::string, std::string*, int);
std::string Join(std::string*, const char, int);
std::string Join(std::string*, int);
void StripComment(std::string&);

#endif
