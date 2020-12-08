#ifndef rdag_io_H
#define rdag_io_H

/* Last edited on 2009-10-28 17:53:42 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <rdag.h>

/* DAG IO */

void rdag_write(FILE *wr, rdag_t *D);
  /* Writes the dag {D} to the given writer, in a format that can be
    read back with {rdag_read} below. The procedure only writes the
    existing nodes, from {0} to {rdag_max_node(D)}. The output is in
    ASCII but not meant to be human-readable. */

rdag_t *rdag_read(FILE *rd);
  /* Reads from {rd} a dag, in the format used by {rdag_write}. */

/* STRING IO AND CLEANUP */

int rdag_line_read (FILE *rd, uint32_t *nbuf, char buf[], uint32_t nmax);
  /* Reads the next line from {rd}, sets {*n} to its length in bytes 
    (excluding the final end-of-line, if any), and stores those bytes into 
    {buf[0..(*n)-1]}.
    
    Assumes that the end-of-line marker is either LF ('\012'), 
    CR ('\015'), or CR-LF. In any case, the entire end-of-line marker is 
    discarded.
    
    In particular, sets {*buf_byte_count} to 0 iff the line was empty 
    (i.e. iff, upon entry, {rd} was positioned just before an end-of-line marker).

    Returns 1 if the read succeeded (even if the line was empty).

    Returns -1 if the there is no next line (i.e. iff, upon entry, the
    file was positioned just after its last character).
    
    The procedure silently provides (and discards) an implied
    end-of-line marker for the last line, if that marker is missing
    (that is, if the file ends with any character other than LF or CR).
    
    Returns -2 the line has more than {max_bytes} bytes (not counting
    the end-of-line marker). In that case {buf_byte_count} is set to
    {max_bytes} and only the first {max_bytes} are stored in {buf}.
    However the entire line is skipped (including the final
    end-of-line marker if present). */

void rdag_string_iso_latin_1_cleanup(uint32_t *nbuf, char buf[]);
  /* Basic cleanup of the string {buf[0..nbuf-1]}. Namely, deletes a
    trailing ASCII NUL ('\000') is present, and converts ASCII whitespace
    characters ('\011'..'\015' = HT, LF, VT, NP, CR) and the ISO
    Latin-1 non-breaking space ('\240') to plain ASCII space (SP =
    '\040'). Then removes leading and trailing whitespace, and
    squeezes any embedded whitespace to a single space. */

bool_t rdag_string_iso_latin_1_encode(uint32_t nbuf, char buf[], uint32_t *nstr, rdag_symbol_t str[], uint32_t nmax);
  /* Converts a string {buf[0..nbuf-1]} of iso-latin-1 character codes
    to a string of {rdag_symbol_t} values. The number of symbols is
    stored in {*nstr} and the symbols in {str[0..(*nstr)-1]}.
    
    Fails if the number of symbols exceeds {nmax}, or the string
    contains any whitespace other than SP, or any non-printable
    iso-latin-1 characters, or the soft hyphen '\255'. In those cases,
    prints a warning message to {stderr}, sets {*nstr} to the number
    of symbols before the error, and returns FALSE. Otherwise returns
    TRUE. */

void rdag_string_iso_latin_1_write(FILE *wr, uint32_t nstr, rdag_symbol_t str[]);
  /* Writes to {wr} a string {str[0..nstr-1]} of input symbols,
    each interpreted as the numeric code of an iso-latin-1 character. */

#endif
