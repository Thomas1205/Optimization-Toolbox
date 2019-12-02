/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#ifndef MAKROS_HH
#define MAKROS_HH

#include <cassert>
#include <iostream>
#include <limits> //includes numeric_limits
#include <sstream>
#include <string>
#include <iomanip>
#include <cstdlib> //includes the exit-function
#include <typeinfo>
#include <cmath>
#include <algorithm>

#include <string.h> //memcpy

#ifdef WIN32
namespace {
  inline bool isnan(double x)
  {
    return (x != x);
  }
}
#define M_PI 3.1415926535897931
#else
using std::isnan;
using std::isinf;
#endif

#ifdef GNU_COMPILER

#define attr_restrict __restrict
#define ref_attr_restrict __restrict
//pointers returned by new are guaranteed to have an address that is divisible by 16 if the type is a basic one
//it is convenient to give the compiler this hint so that he need not handle unaligned cases
#define ALIGNED16 __attribute__ ((aligned (16)))
#define assertAligned16(p) assert( ((size_t)p) % 16 == 0);

#include <execinfo.h>

inline void print_trace (void)
{
  void* array[15];
  size_t size;
  char** strings;
  size_t i;

  size = backtrace (array, 15);
  strings = backtrace_symbols (array, size);

  std::cerr << "Obtained " << size << " stack frames" << std::endl;

  for (i = 0; i < size; i++)
    std::cerr << strings[i] << std::endl;

  free (strings);
}

#else
#define attr_restrict
#define ref_attr_restrict
#define ALIGNED16
#define assertAligned16(p)
inline void print_trace (void) {}
#endif

//because c++ is missing those keywords:
#define abstract
#define overide
#define overrides

/******************** Data Macros *****************************/
typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned char uchar;
typedef double ALIGNED16 double_A16;
typedef float ALIGNED16 float_A16;

#define MIN_DOUBLE -1.0*std::numeric_limits<double>::max()
#define MAX_DOUBLE std::numeric_limits<double>::max()
#define MIN_LONGDOUBLE -1.0*std::numeric_limits<long double>::max()
#define MAX_LONGDOUBLE std::numeric_limits<long double>::max()
#define HIGH_DOUBLE (0.1*MAX_DOUBLE)
#define EPS_DOUBLE std::numeric_limits<double>::epsilon()
#define EPS_LONGDOUBLE std::numeric_limits<long double>::epsilon()
#define MIN_FLOAT  -1.0f*std::numeric_limits<float>::max()
#define MAX_FLOAT  std::numeric_limits<float>::max()
#define HIGH_FLOAT (0.1f*MAX_FLOAT)
#define EPS_FLOAT  std::numeric_limits<float>::epsilon()
#define MAX_UINT std::numeric_limits<uint>::max()
#define MAX_USHORT std::numeric_limits<ushort>::max()

enum NormType {L1,L2,L0_5};

enum DifferenceType {SquaredDiffs,AbsDiffs};

enum RegularityType {SquaredDiffReg,AbsDiffReg,TVReg};


/**** helpful routines ****/

namespace Makros {

  //making log, exp, pow and abs a template is convenient when you want to call the proper function inside your own template

  template<typename T>
  inline T log(T arg)
  {
    return T(::log(double(arg)));
  }

  //specializations:
  template<>
  inline float log(float arg)
  {
    return logf(arg);
  }

  template<>
  inline double log(double arg)
  {
    return ::log(arg);
  }

  template<>
  inline long double log(long double arg)
  {
    return logl(arg);
  }

  template<typename T>
  inline T exp(T arg)
  {
    return T(::exp(double(arg)));
  }

  //specializations:
  template<>
  inline float exp(float arg)
  {
    return expf(arg);
  }

  template<>
  inline double exp(double arg)
  {
    return ::exp(arg);
  }

  template<>
  inline long double exp(long double arg)
  {
    return expl(arg);
  }

  template<typename T>
  inline T pow(T base, T exponent)
  {
    return T(::pow(double(base),double(exponent)));
  }

  //specializations:
  template<>
  inline float pow(float base, float exponent)
  {
    return powf(base,exponent);
  }

  template<>
  inline double pow(double base, double exponent)
  {
    return ::pow(base,exponent);
  }

  template<>
  inline long double pow(long double base, long double exponent)
  {
    return powl(base,exponent);
  }

  template<typename T>
  inline T abs(T arg)
  {
    return std::abs(arg);
  }

  template<>
  inline uchar abs(uchar arg)
  {
    return arg;
  }

  template<>
  inline ushort abs(ushort arg)
  {
    return arg;
  }

  template<>
  inline uint abs(uint arg)
  {
    return arg;
  }

#ifndef _32BIT_OS
  template<>
  inline size_t abs(size_t arg)
  {
    return arg;
  }
#endif

  template<>
  inline float abs(float arg)
  {
    return fabsf(arg);
  }

  template<>
  inline double abs(double arg)
  {
    return fabs(arg);
  }

  template<>
  inline long double abs(long double arg)
  {
    return fabsl(arg);
  }

  template<typename T>
  inline void unified_assign(T* attr_restrict dest, const T* attr_restrict source, size_t size)
  {
    for (size_t i=0; i < size; i++)
      dest[i] = source[i];
  }
  
  template<>
  inline void unified_assign(char* attr_restrict dest, const char* attr_restrict source, size_t size) 
  {
    memcpy(dest, source, size * sizeof(char));
  }  
  
  template<>
  inline void unified_assign(uchar* attr_restrict dest, const uchar* attr_restrict source, size_t size) 
  {
    memcpy(dest, source, size * sizeof(uchar));
  }  
  
  template<>
  inline void unified_assign(short* attr_restrict dest, const short* attr_restrict source, size_t size) 
  {
    memcpy(dest, source, size * sizeof(short));
  }  

  template<>
  inline void unified_assign(ushort* attr_restrict dest, const ushort* attr_restrict source, size_t size) 
  {
    memcpy(dest, source, size * sizeof(ushort));
  }  
  
  template<>
  inline void unified_assign(int* attr_restrict dest, const int* attr_restrict source, size_t size) 
  {
    memcpy(dest, source, size * sizeof(int));
  }  

  template<>
  inline void unified_assign(uint* attr_restrict dest, const uint* attr_restrict source, size_t size) 
  {
    memcpy(dest, source, size * sizeof(uint));
  }  

  template<>
  inline void unified_assign(float* attr_restrict dest, const float* attr_restrict source, size_t size) 
  {
    memcpy(dest, source, size * sizeof(float));
  }  
  
  template<>
  inline void unified_assign(double* attr_restrict dest, const double* attr_restrict source, size_t size) 
  {
    memcpy(dest, source, size * sizeof(double));
  }  

  template<>
  inline void unified_assign(long double* attr_restrict dest, const long double* attr_restrict source, size_t size) 
  {
    memcpy(dest, source, size * sizeof(long double));
  }  
}

template<typename T>
std::string toString(T obj, uint width=1)
{

  std::ostringstream s;

  s << std::setw(width) << std::setfill('0') << obj;
  return s.str();
}

namespace Makros {

  void register_typename(const std::string& id, const std::string& fullname);

  std::string get_typename(const std::string& id);

  template<typename T>
  class Typename {
  public:

    std::string name() const;
  };

  template<typename T>
  std::string Typename<T>::name() const
  {
    return get_typename(typeid(T).name());
  }

  //specializations:

  template<typename T>
  class Typename<const T> {
  public:

    std::string name() const
    {
      return "const " +Typename<T>().name();
    }
  };


  template<typename T>
  class Typename<T*> {
  public:

    std::string name() const
    {
      return Typename<T>().name() + "*";
    }
  };
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const Makros::Typename<T>& t)
{
  out << t.name();
  return out;
}

template<typename T>
std::string operator+(std::string s, const Makros::Typename<T>& t)
{
  return s + t.name();
}

/***********************/

template <typename T>
inline T convert(const std::string s)
{

  std::istringstream is(s);
  T result;

  is >> result;
  if (is.bad() || is.fail()) {
    std::cerr << "ERROR: conversion of \"" << s << "\" to " << Makros::Typename<T>().name()
              << " failed. Exiting." << std::endl;
    exit(1);
  }
  if (!is.eof()) {

    //check if the string contains additional characters that are not whitespace
    char c;
    while (is >> c) {
      if (c != ' ' && c != '\n' && c != 13 && c != 10) {
        std::cerr << "WARNING AFTER CONVERSION: string contains additional characters" << std::endl;
        break;
      }
    }
  }

  return result;
}

template<>
inline uint convert<uint>(const std::string s)
{

  uint result = 0;
  char c;
  uint i=0;
  for (; i < s.size(); i++) {
    c = s[i];

    if (c < '0' || c > '9') {
      std::cerr << "ERROR: conversion of \"" << s << "\" to uint failed. Exiting." << std::endl;
      exit(1);
    }
    result = 10*result + (c - '0');
  }

  return result;
}

template<typename T1, typename T2>
void operator+=(std::pair<T1,T2>& x, const std::pair<T1,T2>& y)
{
  x.first += y.first;
  x.second += y.second;
}


/********************* Code Macros ****************************/
#define TODO(s) { std::cerr << "TODO ERROR[" << __FILE__ << ":" << __LINE__ << "]: feature \"" << (s) << "\" is currently not implemented. Exiting..." << std::endl; exit(1); }
#define EXIT(s) { std::cerr << __FILE__ << ":" << __LINE__ << ": " <<  s << std::endl; exit(1); }
#define MAKENAME(s) std::string(#s) + std::string("[") + std::string(__FILE__) + std::string(":") + toString(__LINE__) + std::string("]")

#ifdef SAFE_MODE
#define OPTINLINE
#else
#define OPTINLINE inline
#endif

#define INTERNAL_ERROR std::cerr << "INTERNAL ERROR[" << __FILE__ << ":" << __LINE__ << "]:" << std::endl
#define USER_ERROR std::cerr << "ERROR: "
#define IO_ERROR std::cerr << "I/O ERROR[" << __FILE__ << ":" << __LINE__ << "]:" << std::endl
#define WARNING std::cerr << "WARNING[" << __FILE__ << ":" << __LINE__ << "]:" << std::endl

template<typename T>
inline T sign(T arg)
{

  if (arg < ((T) 0.0) )
    return ((T) -1.0);
  else if (arg == ((T) 0.0))
    return ((T) 0.0);
  else
    return ((T) 1.0);
}

template<typename T>
inline T robust_sign(T arg, T tolerance)
{

  if (arg < ((T) -tolerance) )
    return ((T) -1.0);
  else if (arg > ((T) tolerance))
    return ((T) 1.0);
  else
    return ((T) 0.0);
}

//NOTE: prefetch will only work on an x86/x86_64 architecture with SSE1 (Pentium 3 or higher)
// if you are compiling on a different architecture simply remove the asm-statements
template<typename T>
inline void prefetcht0(const T* ptr)
{
  asm __volatile__ ("prefetcht0 %[ptr]" : : [ptr] "m" (ptr[0]));
}

template<typename T>
inline void prefetcht1(const T* ptr)
{
  asm ("prefetcht1 %[ptr]" : : [ptr] "m" (ptr[0]));
}

template<typename T>
inline void prefetcht2(const T* ptr)
{
  asm ("prefetcht2 %[ptr]" : : [ptr] "m" (ptr[0]));
}

namespace Makros {

  inline float max(const float_A16* data, size_t nData)
  {
    float max_val=MIN_FLOAT;
    float cur_datum;
    size_t i;

    //#if !defined(USE_SSE) || USE_SSE < 2
#if 1 // g++ 4.8.5 uses avx instructions automatically
    for (i=0; i < nData; i++) {
      cur_datum = data[i];
      max_val = std::max(max_val,cur_datum);
    }
#else
    //movups is part of SSE2

    float tmp[4] = {MIN_FLOAT,MIN_FLOAT,MIN_FLOAT,MIN_FLOAT};
    const float* fptr;

    asm __volatile__ ("movups %[tmp], %%xmm6" : : [tmp] "m" (tmp[0]) : "xmm6");
    for (i=0; (i+4) <= nData; i += 4) {
      fptr = data+i;
      asm __volatile__ ("movups %[fptr], %%xmm7\n\t"
                        "maxps %%xmm7, %%xmm6" : : [fptr] "m" (fptr[0]) : "xmm6", "xmm7");

    }
    asm __volatile__ ("movups %%xmm6, %[tmp]" : [tmp] "=m" (tmp[0]) : : );
    for (i=0; i < 4; i++)
      max_val = std::max(max_val,tmp[i]);

    for (i= nData - (nData % 4); i < nData; i++) {
      cur_datum = data[i];
      if (cur_datum > max_val)
        max_val = cur_datum;
    }
#endif

    return max_val;
  }

  inline float min(const float_A16* data, size_t nData)
  {
    float min_val=MAX_FLOAT;
    float cur_datum;
    size_t i;

    //#if !defined(USE_SSE) || USE_SSE < 2
#if 1 // g++ 4.8.5 uses avx instructions automatically
    for (i=0; i < nData; i++) {
      cur_datum = data[i];
      min_val = std::min(min_val,cur_datum);
    }
#else
    //movups is part of SSE2

    float tmp[4] = {MAX_FLOAT,MAX_FLOAT,MAX_FLOAT,MAX_FLOAT};
    const float* fptr;

    asm __volatile__ ("movups %[tmp], %%xmm6" : : [tmp] "m" (tmp[0]) : "xmm6");
    for (i=0; (i+4) <= nData; i += 4) {
      fptr = data+i;
      asm __volatile__ ("movups %[fptr], %%xmm7 \n\t"
                        "minps %%xmm7, %%xmm6 \n\t" : : [fptr] "m" (fptr[0]) : "xmm6","xmm7");
    }
    asm __volatile__ ("movups %%xmm6, %[tmp]" : [tmp] "=m" (tmp[0]) :  : "xmm6");
    for (i=0; i < 4; i++)
      min_val = std::min(min_val,tmp[i]);

    for (i= nData - (nData % 4); i < nData; i++) {
      cur_datum = data[i];
      if (cur_datum < min_val)
        min_val = cur_datum;
    }
#endif

    return min_val;
  }


  inline void find_max_and_argmax(const float_A16* data, const size_t nData, float& max_val, size_t& arg_max)
  {

    max_val = MIN_FLOAT;
    arg_max = MAX_UINT;

#if !defined(USE_SSE) || USE_SSE < 4

    if (nData > 0) {
      const float* ptr = std::max_element(data,data+nData);
      max_val = *ptr;
      arg_max = ptr - data;
    }

    // for (i=0; i < nData; i++) {
    //   cur_val = data[i];

    //   if (cur_val > max_val) {
    //     max_val = cur_val;
    //     arg_max = i;
    //   }
    // }
#else
    //blendvps is part of SSE4

    size_t i;
    float cur_val;

    assert(nData <= 17179869183);

    float tmp[4] = {MIN_FLOAT,MIN_FLOAT,MIN_FLOAT,MIN_FLOAT};
    const float* fptr;

    wchar_t itemp[4] = {1,1,1,1};

    asm __volatile__ ("movups %[tmp], %%xmm6 \n\t"
                      "xorps %%xmm5, %%xmm5 \n\t" //sets xmm5 (= argmax) to zero
                      "movups %[itemp], %%xmm4 \n\t"
                      "xorps %%xmm3, %%xmm3 \n\t"
                      : : [tmp] "m" (tmp[0]), [itemp] "m" (itemp[0]) : "xmm3", "xmm4", "xmm5", "xmm6");

    for (i=0; (i+4) <= nData; i += 4) {
      fptr = data+i;

      asm __volatile__ ("movups %[fptr], %%xmm7 \n\t"
                        "movaps %%xmm7, %%xmm0 \n\t"
                        "cmpnleps %%xmm6, %%xmm0 \n\t"
                        "blendvps %%xmm7, %%xmm6 \n\t"
                        "blendvps %%xmm3, %%xmm5 \n\t"
                        "paddd %%xmm4, %%xmm3 \n\t"
                        : : [fptr] "m" (fptr[0]) : "xmm0", "xmm3", "xmm5", "xmm6", "xmm7");
    }

    asm __volatile__ ("movups %%xmm6, %[tmp] \n\t"
                      "movups %%xmm5, %[itemp]"
                      : [tmp] "=m" (tmp[0]), [itemp] "=m" (itemp[0]) : : );

    for (i=0; i < 4; i++) {
      cur_val = tmp[i];
      if (cur_val > max_val) {
        max_val = cur_val;
        arg_max = 4*itemp[i] + i;
      }
    }

    for (i= nData - (nData % 4); i < nData; i++) {
      cur_val = data[i];
      if (cur_val > max_val) {
        max_val = cur_val;
        arg_max = i;
      }
    }
#endif
  }

  inline void find_max_and_argmax(const double_A16* data, const size_t nData, double& max_val, size_t& arg_max)
  {

    max_val = MIN_DOUBLE;
    arg_max = MAX_UINT;

#if !defined(USE_SSE) || USE_SSE < 4
    if (nData > 0) {
      const double* ptr = std::max_element(data,data+nData);
      max_val = *ptr;
      arg_max = ptr - data;
    }

    // for (i=0; i < nData; i++) {
    //   cur_val = data[i];

    //   if (cur_val > max_val) {
    //     max_val = cur_val;
    //     arg_max = i;
    //   }
    // }
#else

    size_t i;
    double cur_val;

    assert(nData < 8589934592);

    volatile double tmp[2] = {MIN_DOUBLE,MIN_DOUBLE};
    const double* dptr;

    volatile wchar_t itemp[4] = {0,1,0,1};


    asm __volatile__ ("movupd %[tmp], %%xmm6 \n\t"
                      "xorpd %%xmm5, %%xmm5 \n\t" //sets xmm5 (= argmin) to zero
                      "movupd %[itemp], %%xmm4 \n\t"
                      "xorpd %%xmm3, %%xmm3 \n\t" //contains candidate argmin
                      : : [tmp] "m" (tmp[0]), [itemp] "m" (itemp[0]) :  "xmm3", "xmm4", "xmm5", "xmm6");

    for (i=0; (i+4) <= nData; i += 4) {
      dptr = data+i;

      asm __volatile__ ("movupd %[dptr], %%xmm7 \n\t"
                        "movapd %%xmm7, %%xmm0 \n\t"
                        "cmpnlepd %%xmm6, %%xmm0 \n\t"
                        "blendvpd %%xmm7, %%xmm6 \n\t"
                        "blendvpd %%xmm3, %%xmm5 \n\t"
                        "paddd %%xmm4, %%xmm3 \n\t"
                        : : [dptr] "m" (dptr[0]) : "xmm0", "xmm3", "xmm5", "xmm6", "xmm7");
    }

    asm __volatile__ ("movupd %%xmm6, %[tmp] \n\t"
                      "movupd %%xmm5, %[itemp] \n\t"
                      : [tmp] "=m" (tmp[0]), [itemp] "=m" (itemp[0]) : : "xmm5", "xmm6");

    assert(itemp[0] == 0);
    assert(itemp[2] == 0);

    for (i=0; i < 2; i++) {
      cur_val = tmp[i];
      //std::cerr << "cur val: " << cur_val << std::endl;
      if (cur_val > max_val) {
        max_val = cur_val;
        arg_max = 2*itemp[2*i+1] + i;
      }
    }

    //std::cerr << "minval: " << min_val << std::endl;

    if ((nData % 2) == 1) {
      cur_val = data[nData-1];
      if (cur_val > max_val) {
        max_val = cur_val;
        arg_max = nData-1;
      }
    }


#endif
  }


  inline void find_min_and_argmin(const float_A16* data, const size_t nData, float& min_val, size_t& arg_min)
  {

    min_val = MAX_FLOAT;
    arg_min = MAX_UINT;

#if !defined(USE_SSE) || USE_SSE < 4

    if (nData > 0) {
      const float* ptr = std::min_element(data,data+nData);
      min_val = *ptr;
      arg_min = ptr - data;
    }

    // for (i=0; i < nData; i++) {
    //   cur_val = data[i];

    //   if (cur_val < min_val) {
    //     min_val = cur_val;
    //     arg_min = i;
    //   }
    // }
#else
    //blendvps is part of SSE4

    size_t i;
    float cur_val;

    assert(nData <= 17179869183);

    volatile float tmp[4] = {MAX_FLOAT,MAX_FLOAT,MAX_FLOAT,MAX_FLOAT};
    const float* fptr;

    volatile wchar_t itemp[4] = {1,1,1,1};

    asm __volatile__ ("movups %[tmp], %%xmm6 \n\t"
                      "xorps %%xmm5, %%xmm5 \n\t" //sets xmm5 (= argmin) to zero
                      "movups %[itemp], %%xmm4 \n\t"
                      "xorps %%xmm3, %%xmm3 \n\t" //contains candidate argmin
                      : : [tmp] "m" (tmp[0]), [itemp] "m" (itemp[0]) :  "xmm3", "xmm4", "xmm5", "xmm6");

    for (i=0; (i+4) <= nData; i += 4) {
      fptr = data+i;

      asm __volatile__ ("movups %[fptr], %%xmm7 \n\t"
                        "movaps %%xmm7, %%xmm0 \n\t"
                        "cmpltps %%xmm6, %%xmm0 \n\t"
                        "blendvps %%xmm7, %%xmm6 \n\t"
                        "blendvps %%xmm3, %%xmm5 \n\t"
                        "paddd %%xmm4, %%xmm3 \n\t"
                        : : [fptr] "m" (fptr[0]) : "xmm0", "xmm3", "xmm5", "xmm6", "xmm7");
    }

    asm __volatile__ ("movups %%xmm6, %[tmp] \n\t"
                      "movups %%xmm5, %[itemp] \n\t"
                      : [tmp] "=m" (tmp[0]), [itemp] "=m" (itemp[0]) : : "xmm5", "xmm6");

    //std::cerr << "intermediate minval: " << min_val << std::endl;

    for (i=0; i < 4; i++) {
      cur_val = tmp[i];
      //std::cerr << "cur val: " << cur_val << std::endl;
      if (cur_val < min_val) {
        min_val = cur_val;
        arg_min = 4*itemp[i] + i;
      }
    }

    //std::cerr << "minval: " << min_val << std::endl;

    for (i= nData - (nData % 4); i < nData; i++) {
      cur_val = data[i];
      if (cur_val < min_val) {
        min_val = cur_val;
        arg_min = i;
      }
    }
#endif
  }

  inline void find_min_and_argmin(const double_A16* data, const size_t nData, double& min_val, size_t& arg_min)
  {

    min_val = MAX_DOUBLE;
    arg_min = MAX_UINT;

#if !defined(USE_SSE) || USE_SSE < 4

    if (nData > 0) {
      const double* ptr = std::min_element(data,data+nData);
      min_val = *ptr;
      arg_min = ptr - data;
    }

    // for (i=0; i < nData; i++) {
    //   cur_val = data[i];

    //   if (cur_val < min_val) {
    //     min_val = cur_val;
    //     arg_min = i;
    //   }
    // }
#else

    size_t i;
    double cur_val;

    assert(nData < 8589934592);

    volatile double tmp[2] = {MAX_DOUBLE,MAX_DOUBLE};
    const double* dptr;

    volatile wchar_t itemp[4] = {0,1,0,1};


    asm __volatile__ ("movupd %[tmp], %%xmm6 \n\t"
                      "xorpd %%xmm5, %%xmm5 \n\t" //sets xmm5 (= argmin) to zero
                      "movupd %[itemp], %%xmm4 \n\t"
                      "xorpd %%xmm3, %%xmm3 \n\t" //contains candidate argmin
                      : : [tmp] "m" (tmp[0]), [itemp] "m" (itemp[0]) :  "xmm3", "xmm4", "xmm5", "xmm6");

    for (i=0; (i+4) <= nData; i += 4) {
      dptr = data+i;

      asm __volatile__ ("movupd %[dptr], %%xmm7 \n\t"
                        "movapd %%xmm7, %%xmm0 \n\t"
                        "cmpltpd %%xmm6, %%xmm0 \n\t"
                        "blendvpd %%xmm7, %%xmm6 \n\t"
                        "blendvpd %%xmm3, %%xmm5 \n\t"
                        "paddd %%xmm4, %%xmm3 \n\t"
                        : : [dptr] "m" (dptr[0]) : "xmm0", "xmm3", "xmm5", "xmm6", "xmm7");
    }

    asm __volatile__ ("movupd %%xmm6, %[tmp] \n\t"
                      "movupd %%xmm5, %[itemp] \n\t"
                      : [tmp] "=m" (tmp[0]), [itemp] "=m" (itemp[0]) : : "xmm5", "xmm6");

    assert(itemp[0] == 0);
    assert(itemp[2] == 0);

    for (i=0; i < 2; i++) {
      cur_val = tmp[i];
      //std::cerr << "cur val: " << cur_val << std::endl;
      if (cur_val < min_val) {
        min_val = cur_val;
        arg_min = 2*itemp[2*i+1] + i;
      }
    }

    //std::cerr << "minval: " << min_val << std::endl;

    if ((nData % 2) == 1) {
      cur_val = data[nData-1];
      if (cur_val < min_val) {
        min_val = cur_val;
        arg_min = nData-1;
      }
    }

#endif
  }


  inline void mul_array(float_A16* data, const size_t nData, const float constant)
  {

    size_t i;
#if !defined(USE_SSE) || USE_SSE < 2
    for (i=0; i < nData; i++) {
      data[i] *= constant;
    }
#else
    float temp[4];
    float* fptr;
    for (i=0; i < 4; i++)
      temp[i] = constant;
    asm volatile ("movups %[temp], %%xmm7" : : [temp] "m" (temp[0]) : "xmm7" );

    for (i=0; i+4 <= nData; i+=4) {
      fptr = data + i;
      asm volatile ("movups %[fptr], %%xmm6 \n\t"
                    "mulps %%xmm7, %%xmm6 \n\t"
                    "movups %%xmm6, %[fptr] \n\t"
                    : [fptr] "+m" (fptr[0]) : : "xmm6");
    }

    for (i= nData - (nData % 4); i < nData; i++) {
      data[i] *= constant;
    }
#endif
  }

  inline void mul_array(double_A16* data, const size_t nData, const double constant)
  {

    size_t i;
#if !defined(USE_SSE) || USE_SSE < 2
    for (i=0; i < nData; i++) {
      data[i] *= constant;
    }
#else
    double temp[2];
    double* dptr;
    for (i=0; i < 2; i++)
      temp[i] = constant;
    asm volatile ("movupd %[temp], %%xmm7" : : [temp] "m" (temp[0]) : "xmm7" );

    for (i=0; i+2 <= nData; i+=2) {
      dptr = data + i;

      asm volatile ("movupd %[dptr], %%xmm6 \n\t"
                    "mulpd %%xmm7, %%xmm6 \n\t"
                    "movupd %%xmm6, %[dptr] \n\t"
                    : [dptr] "+m" (dptr[0]) : : "xmm6");
    }

    for (i= nData - (nData % 2); i < nData; i++) {
      data[i] *= constant;
    }
#endif
  }

  //performs data[i] -= factor*data2[i] for each i
  //this is a frequent operation in the conjugate gradient algorithm
  inline void array_subtract_multiple(double_A16* attr_restrict data, const size_t nData, double factor,
                                      const double_A16* attr_restrict data2)
  {

    size_t i;
#if !defined(USE_SSE) || USE_SSE < 2
    for (i=0; i < nData; i++)
      data[i] -= factor*data2[i];
#else
    double temp[2];
    double* dptr;
    const double* cdptr;
    for (i=0; i < 2; i++)
      temp[i] = factor;
    asm volatile ("movupd %[temp], %%xmm7" : : [temp] "m" (temp[0]) : "xmm7" );
    for (i=0; i+2 <= nData; i+=2) {
      cdptr = data2+i;
      dptr = data+i;

      asm volatile ("movupd %[cdptr], %%xmm6 \n\t"
                    "mulpd %%xmm7, %%xmm6 \n\t"
                    "movupd %[dptr], %%xmm5 \n\t"
                    "subpd %%xmm6, %%xmm5 \n\t"
                    "movupd %%xmm5, %[dptr] \n\t"
                    : [dptr] "+m" (dptr[0]) : [cdptr] "m" (cdptr[0]) : "xmm5", "xmm6");
    }

    for (i= nData - (nData % 2); i < nData; i++)
      data[i] -= factor*data2[i];
#endif
  }

  template<typename T1, typename T2>
  class first_lower {
  public:
    bool operator()(const std::pair<T1,T2>& p1, const std::pair<T1,T2>& p2)
    {
      return (p1.first < p2.first);
    }
  };

  template<typename T1, typename T2>
  class first_higher {
  public:
    bool operator()(const std::pair<T1,T2>& p1, const std::pair<T1,T2>& p2)
    {
      return (p1.first > p2.first);
    }
  };


} //end of namespace Makros



#endif
