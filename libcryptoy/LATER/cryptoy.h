/* Toy cryptography tools. */
/* Last edited on 2013-09-04 18:15:07 by stolfi */

#ifndef cryptoy_H
#define cryptoy_H

#include <stdio.h>

#define cryptoy_H_copyright \
  "Copyright © 2013 Jorge Stolfi, State University of Campinas"
  
/* FILE PROCESSING TOOLS */
  
#define cryptoy_emsg(m) \
  fprintf(stderr, "%s:%s: ** error - %s\n", __FILE__, __FUNCTION__, (m))
  
#define cryptoy_wmsg(m) \
  fprintf(stderr, "%s:%s: !! warning - %s\n", __FILE__, __FUNCTION__, (m))

#define cryptoy_ERROR (-2)

/* The procedures in this section require that all file arguments be
  opened by the caller in the correct direction (reading of writing).
  In case of failure (including read or write errors) they print a
  warning to {stderr} and return {cryptoy_ERROR}; otherwise they
  return the number of bytes processed. Output files are flushed with
  {fflush}, but all files are left open. The bits in each byte are
  read or written in big-to-small order (from highest-valued to
  lowest-valued). */

int cryptoy_xor_files(FILE *f0, FILE *f1, FILE *f);
  /* Writes for file {f} the bitwise XOR of the contents of files {f0}
    and {f1}. The files must be properly opened by the caller. Fails
    if file {f0} is longer than file {f1}. */

int cryptoy_bitsplit_file(FILE *f, FILE *k, FILE *f0, FILE *f1);
  /* Reads file {f} and splits it into two files {f0,f1}, bit by bit,
    according to the bits read from the key file {k}. Specifically,
    each bit of {f} goes into {f0} or to {f1} depending on whether the
    corresponding bit of {k} is 0 or 1, respectively. Fails (returning
    {cryptoy_ERROR}) if file {f} is longer than file {k}.
    
    When file {f} is exhausted, the last bytes of file {f0} and {f1}
    may be both partly filled. In that case, the last incomplete byte
    of {f1} is not written, and its bits are appended to the last
    incomplete byte of {f0}, completing it. Thus, the last 1 to 7 bits
    of file {f0} may be inconsistent with the contents of the key file {k}.
    
    This special case does not occur if, for example, every {2*m}
    bytes of {k} contain equal numbers of 0s and 1s, and the
    byte-length of file {f} is a multiple of {2*m}. Then files {f0}
    and {f1} will have the same byte-length, which will be a multiple
    of {m}.
    
    The effect of this procedure can be exactly reversed with
    {cryptoy_bitmerge_files}. */

int cryptoy_bitmerge_files(FILE *f0, FILE *f1, FILE *k, FILE *f);
  /* Reads files {f0} and {f1} and merges them into a single file {f},
    bit by bit, according to the bits read from the key file {k}.
    Specifically, each bit of {f} is copied from {f0} or from {f1}
    depending on whether the corresponding bit of {k} is 0 or 1,
    respectively.
    
    Once either file is exhausted, then the remainder of the other
    file is copied to {f}, ignoring the remainder of the key file {k}.
    May fail (returning {cryptoy_ERROR}) if files {f0} and {f1}
    together are longer than file {k}.
    
    If {f0} and {f1} were created by
    {cryptoy_bitsplit_file(g,k,f0,f1)}, then
    {cryptoy_bitmerge_files(f0,f1,k,f)}, with the same key file {k},
    will succeed and the file {f} will be identical to {g}.
    
    More generally, the operation {cryptoy_bitmerge_files(f0,f1,k,f)}
    can be undone with {cryptoy_bitsplit_file(f,k,g0,g1)} only if the
    lengths of {f0} and {f1} are ``nearly consistent'' with the number
    of zeros and ones in the {k} file. Specifically, if the procedure
    reaches the end of file {f1} when the next bit of file {k} is 1,
    then file {f0} must be exhausted too, or must have at most 7 bits
    left (which, as per above, are copied to {f} ignoring the key
    {k}). In that case, {g0,g1} will be identical to the merged files
    {f0,f1}. In any other case, {cryptoy_bitsplit_file(f,k,g0,g1)}
    will split the trailing bits of {f} incorrectly.
    
    This ``nearly consistent'' condition is always satisfied for files
    {f0,f1} created with {cryptoy_bitsplit_file(f,k,f0,f1)}. It is
    also satisfied if every {2*m} bytes of {k} contain equal numbers
    of 0s and 1s, and {f0,f1} have the same length (in bytes) that is
    divisible by {m}. */

#define cryptoy_H_rights \
  cryptoy_H_copyright ".\n" \
  "This file is provided 'as is'; the author and his employer are" \
  " not liable for any damages that may result from its use.  This" \
  " file can be used and modified for any purpose as long as" \
  " the copyright and this 'rights' note are preserved in" \
  " all copies and derived versions.";

#endif
