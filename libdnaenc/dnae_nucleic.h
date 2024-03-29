#ifndef dnae_nucleic_H
#define dnae_nucleic_H

/* Raw DNA/RNA sequences */
/* Last edited on 2022-10-31 11:23:02 by stolfi */

#define dnae_nucleic_H_COPYRIGHT \
  "Copyright � 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdint.h>

#include <vec.h>

#include <msm_rung.h>

bool_t is_dna_basis(char c);
  /* TRUE iff {c} is a valid nucleotide letter, namely one of [ATCGUatcgu]. */

/* CONVERSION BETWEEN LETTERS AND SAMPLES */

void dnae_nucleic_value(char b, int *A, int *T, int *C, int *G);
  /* If {b} is a nucleotide letter code in [ATCGUatcgu], sets
    the appropriate variable {A,T,C,G} to 1, and the rest to 0.
    Otherwise sets all four variables to 0.  Ignores the 
    case of {b}. Treats 'U' and 'u' as equivalent to 'T'. */

/* RANDOM NUCLEOTIDE SEQUENCES */

char dnae_nucleic_throw(bool_t upper, bool_t dna);
  /*  Returns a random nucleotide code.  
    The code will be uppercase iff {upper} is true.  The result is chosen unformly 
    from [ATCG] or [atcg] if {dna} is TRUE, or from [AUCG] or [aucg] if {dna} is FALSE. */

char *dnae_nucleic_string_throw(int32_t nb, bool_t upper, bool_t dna);
  /* Generates a string containing a random sequence of {nb} nucleotide
    codes.  Each letter is generated with {dnae_nucleic_throw(upper, dna)}. */
    
void dnae_nucleic_string_mutate
  ( char *s, 
    bool_t dna, 
    bool_t upper, 
    double mutProb, 
    double delProb,
    char **rp,
    msm_rung_vec_t *gv
  );
  /* Generates a mutated copy {r} of nucleotide string {s},
    and puts its address in {*rp}.
  
    The new string {r} is built one nucleotide at a time. With probability
    {mutProb}, the next nucleotide is replaced by a random nucleotide.
    With probability {delProb/2}, the next nucleotide is deleted. With
    probabilty {delProb/2}, a random nucleotide is inserted before the
    next one. The random nucleotides are generated by
    {dnae_nucleic_throw(upper,dna)}.
    
    The procedure also returns in {*gv} a list of the index pairs
    {(ix,iy)} such that nucleotide {s[ix]} was copied without mutation
    to {r[iy]}. The string {r} and the vector {*gv} are allocated by
    the procedure. */
  
/* NUCLEOTIDE SEQUENCE I/O */

void dnae_nucleic_string_write(FILE *wr, char *s, char *cmt);
  /* Writes a DNA/RNA sequence {s} to file {wr}. The file format is
    described by {dnae_nucleic_file_format_INFO} below. The comment
    {cmt} is written out at the top of the file. */

void dnae_nucleic_string_write_named(char *fname, char *tag, char *ext, char *s, char *cmt);
  /* Same as {dnae_nucleic_string_write} except that writes to a disk file
    called "{fname}{tag}{ext}". */

void dnae_nucleic_string_read(FILE *rd, char **sp, char **cmtp);
  /* Reads a DNA/RNA sequence from file {rd}, and returns its bases
    as the NUL-terminated {*basp}. The file format is described by
    {dnae_nucleic_file_format_INFO} below. The comment lines in the
    file, if any are returned in {*cmtp}. Storage for both strings is
    allocated by the procedure. */

void dnae_nucleic_string_read_named(char *fname, char *tag, char *ext, char **sp, char **cmtp);
  /* Same as {dnae_nucleic_string_read} except that reads from a disk file
    called "{fname}{tag}{ext}". */

#define dnae_nucleic_file_format_INFO \
  "  A DNA/RNA file starts with zero or more comment" \
  " lines beginning with \"> \".  The remaining lines" \
  " contain a sequence of nucleotide codes" \
  " \"A\", \"T\", \"C\", \"G\", or \"U\", in upper or" \
  " lower case.  Blanks, tabs, newlines, and" \
  " carriage-returns are not significant."

#endif
