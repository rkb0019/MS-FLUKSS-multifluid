#ifndef  MHDAM_UTILS_H
#define  MHDAM_UTILS_H

#ifndef  MHDAM_DEFS_H
#include "MHDAMDefs.H"
#endif   // MHDAM_DEFS_H

template <class T> inline
bool isEven( T i )   { return ((i % 2) == 0); }

template <class T> inline
bool isOdd ( T i )   { return ((i % 2) != 0); }

template <class T> inline
T isInRange( T a, T b, T c )
{
  return ((b < a) ? a : ((b < c) ? b : c ));
}

template <class T> inline
int compare( T a, T b )
{
  return ((b == a) ? 0 : ((a < b) ? -1 : 1 ));
}

template <class T> inline
void swap2( T & x, T & y )
{
  T tmp  = x;
  x      = y;
  y      = tmp;
}

template <class T> inline
void swap3( T & x, T & y, T & z )
{
  T tmp  = x;
  x      = y;
  y      = z;
  z      = tmp;
}

template <class T> inline
T absolute( T a )      { return ((a < (T)0) ? -a : a); }

template <class T> inline
bool isSameSign( T a, T b ) { return (((a > 0) && (b > 0)) || ((a < 0) && (b < 0))); }

template <class T> inline
T maximum( T a, T b ) { return ((a >= b) ? a : b); }

template <class T> inline
T maximum( T a, T b, T c )
{
  T maxValue = a;
  if( b > maxValue ) maxValue  = b;
  if( c > maxValue ) maxValue  = c;
  return (maxValue);
}

template <class T> inline
T maximum( T a, T b, T c, T d )
{
  T maxValue = a;
  if( b > maxValue ) maxValue  = b;
  if( c > maxValue ) maxValue  = c;
  if( d > maxValue ) maxValue  = d;
  return (maxValue);
}

template <class T> inline
T minimum( T a, T b ) { return ((a <= b) ? a : b); }

template <class T> inline
T minimum( T a, T b, T c )
{
  T minValue = a;
  if( b < minValue ) minValue  = b;
  if( c < minValue ) minValue  = c;
  return (minValue);
}

template <class T> inline
T minimum( T a, T b, T c, T d )
{
  T minValue = a;
  if( b < minValue ) minValue  = b;
  if( c < minValue ) minValue  = c;
  if( d < minValue ) minValue  = d;
  return (minValue);
}

template <class T> inline
T minmod( T a, T b )
{
  return ((a > (T)0) ? ((b <= (T)0) ? (T)0 : ((a > b) ? b : a )) :
                       ((b >= (T)0) ? (T)0 : ((a < b) ? b : a )));
}

template <class T> inline
T signum  ( T a )      { return ((a < (T)0) ? (T)(-1) : (T)1); }

template <class T> inline
T signum  ( T a, T b )
{
  if( b >= T(0) )
    return ((a < (T)0) ? -a : a);
  else
    return ((a < (T)0) ? a : -a);
}

#endif // MHDAM_UTILS_H
