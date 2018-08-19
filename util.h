#ifndef UTIL_H
#define UTIL_H

#include <fstream>

inline bool readLine( std::ifstream &in, float f[], int n)
{
  char c;
  for( int i=0; i<n-1; i++ ) in >> f[i] >> c;
  in >> f[n-1];
  return in.good();
}

inline bool readLine( std::ifstream &in, double f[], int n)
{
  char c;
  for( int i=0; i<n-1; i++ ) in >> f[i] >> c;
  in >> f[n-1];
  return in.good();
}

#endif
