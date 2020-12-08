/* In-core byte-oriented toy cryptography tools. */
/* Last edited on 2013-10-31 01:43:07 by stolfilocal */

#ifndef cryptoy_mem_H
#define cryptoy_mem_H

#include <stdio.h>
#include <stdint.h>

#define cryptoy_mem_H_copyright \
  "Copyright © 2013 Jorge Stolfi, State University of Campinas"

typedef uint8_t byte_t;
  /* An 8-bit byte, also an integer in {0..255}. */
  
/* In all operations of this interface, each output area must either 
  coincide exactly with an input area, or must be completely disjoint 
  from it; otherwise the result is undefined. */

void cryptoy_mem_xor(size_t n, byte_t a[], byte_t b[], byte_t r[]);
  /* Stores into bytes {r[0..n-1]} the exclusive OR of bytes {a[0..n-1]}
    and {b[0..n-1]}. The input areas had better be distinct.
    The input data can be recovered by {cryptoy_mem_xor(n,a,r,b)}
    or {cryptoy_mem_xor(n,r,b,a)}. */
    
void cryptoy_mem_mix(size_t n, byte_t a[], byte_t b[], byte_t x[], byte_t s[], byte_t r[]);
  /* Stores into bytes {r[0..n-1],s[0..n-1]} the bytes
    {a[0..n-1],b[0..n-1]}, respectively, after swapping any bits that
    correspond to '1' bits of {x[0..n-1]}. The three input areas must
    be distinct, and ditto for the two output areas.  The input data
    can be recovered by {cryptoy_mem_mix(n,r,s,x,a,b)}. */

/* RIGHTS AND DISCLAIMERS */

#define cryptoy_mem_H_rights \
  "This file can be used and modified for any purpose as long as the source" \
  " is made available with any compiled code and the copyright and this" \
  " 'rights' note are preserved in all copies and derived versions."

#define cryptoy_mem_H_warranty \
  "This file is provided 'as is'. The author and his employer are" \
  " not liable for any damages that may result from its use."

#endif
