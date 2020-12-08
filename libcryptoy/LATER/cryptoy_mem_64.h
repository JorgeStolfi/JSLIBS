/* In-core 64-bit-oriented toy cryptography tools. */
/* Last edited on 2013-10-30 15:20:18 by stolfilocal */

#ifndef cryptoy_mem_64_H
#define cryptoy_mem_64_H

#include <stdio.h>
#include <stdint.h>

#define cryptoy_H_copyright \
  "Copyright © 2013 Jorge Stolfi, State University of Campinas"

void cryptoy_mem_xor_64(size_t n, uint64_t a[], uint64_t b[], uint64_t r[]);
  /* Stores into the {n} 64-bit words {c[0..n-1]} the exclusive OR of words {a[0..n-1]}
    and {b[0..n-1]}.  Pairwise, these three memory areas must be either 
    coincident or disjoint. */
    
void 


/* RIGHTS AND DISCLAIMERS */

#define cryptoy_mem_H_rights \
  "This file can be used and modified for any purpose as long as the source" \
  " is made available with any compiled code and the copyright and this" \
  " 'rights' note are preserved in all copies and derived versions."

#define cryptoy_mem_H_warranty \
  "This file is provided 'as is'. The author and his employer are" \
  " not liable for any damages that may result from its use."

#endif
