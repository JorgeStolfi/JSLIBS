/* jswsize.h -- general definitions that depend on the word size. */
/* Last edited on 2011-12-22 00:49:39 by stolfi */

#ifndef jswsize_H
#define jswsize_H

#define _GNU_SOURCE
#include <stdint.h>

#if __WORDSIZE == 64
typedef uint64_t uaddress_t; 
  /* An unsigned integer type with the same size as an address. */
#define uaddress_u_fmt "lu"
  /* Format to print an {uaddress_t} in decimal. */
#define uint64_u_fmt "lu"
  /* Format to print an {uint64_t} in decimal. */
#define int64_d_fmt "ld"
  /* Format to print an {int64_t} in decimal. */
#endif

#if __WORDSIZE == 32
typedef uint32_t uaddress_t;
  /* An unsigned integer type with the same size as an address. */
#define uaddress_u_fmt "llu"
  /* Format to print an {uaddress_t} in decimal. */
#define uint64_u_fmt "llu"
  /* Format to print an {uint64_t} in decimal. */
#define int64_d_fmt "lld"
  /* Format to print an {int64_t} in decimal. */
#endif

#endif
