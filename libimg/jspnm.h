/* jspnm.h - basic definitions for reading/writing PBM/PGM/PBM files. */
/* Last edited on 2024-12-04 23:34:28 by stolfi */ 

#ifndef jspnm_H
#define jspnm_H

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

#include <bool.h>

#define PNM_MAX_SAMPLE 65535u
  /* Maximum value that can be stored in a {uint16_t}. */

typedef uint16_t pnm_format_t;
  /* The ``magic number'' (first two bytes) of a PBM/PGM/PPM file,
    packed as a big-endian short. */

/* ********************************************************************** */
/* PBM (BILEVEL) FILE FORMAT */

#define PBM_FORMAT ('P' * 256 + '1')
#define RPBM_FORMAT ('P' * 256 + '4')
  /* Magic numbers of plain and raw PBM file formats. */

/* ********************************************************************** */
/* PGM (GRAYSCALE) FILE FORMAT */

#define PGM_FORMAT ('P' * 256 + '2')
#define RPGM_FORMAT ('P' * 256 + '5')
  /* Magic numbers of plain and raw PGM file formats. */

/* ********************************************************************** */
/* PPM (COLOR) FILE FORMAT */

#define PPM_FORMAT ('P' * 256 + '3')
#define RPPM_FORMAT ('P' * 256 + '6')
  /* Magic numbers of plain and raw PPM file formats. */

#define PNM_FILE_MAX_MAXVAL 65535u
  /* The maximum {maxval} that can be used in PGM/PPM file, raw or plain. */

/* ********************************************************************** */
/* FILE HEADER I/O */

uint16_t pnm_file_max_maxval(pnm_format_t format);
  /* The maximum {maxval} allowed in a file with the given {format}. */

void pnm_choose_output_format
  ( uint16_t maxval, 
    int32_t chns, 
    bool_t forceplain,
    pnm_format_t *formatP,
    bool_t *rawP,
    bool_t *bitsP
  );
  /* Chooses a suitable output format {*formatP} for an image with {chns}
    channels and the given {maxval}. The result will be 

      {PBM_FORMAT} or {RPBM_FORMAT} if {chns == 1} and {maxval == 1}
      {PGM_FORMAT} or {RPGM_FORMAT} if {chns == 1} and {maxval > 1}
      {PPM_FORMAT} or {RPPM_FORMAT} if {chns == 3}

    The procedure also sets {*rawP = TRUE} iff `raw' version was
    chosen, and {*bitsP = TRUE} iff one of the PBM versions was
    chosen. It chooses the `raw' version if {forceplain = FALSE} and
    {maxval} small enough; otherwise it chooses the plain version. */

void pnm_read_header
  ( FILE *rd, 
    int32_t *colsP, 
    int32_t *rowsP, 
    int32_t *chnsP, 
    uint16_t *maxvalP, 
    bool_t *rawP, 
    bool_t *bitsP,
    pnm_format_t *formatP
  );
  /* Reads the header of a PBM/PGM/PPM file from {rd}.  Returns its parameters in 
    {*colsP,*rowsP,*maxvalP,*formatP}.  Also returns the number of channels
    (1 or 3) in {*chnsP}, and sets {*rawP} to TRUE iff the format is a `raw'
    variant (with binary samples).  Also sets {*bitsP} iff the format is 
    a PBM format (raw or plain). */

void pnm_write_header(FILE *wr, int32_t cols, int32_t rows, uint16_t maxval, pnm_format_t format);
  /* Writes the header of a PBM/PGM/PPM file to {wr}, with the given dimensions,
    {maxval}, and {format}. */

/* ********************************************************************** */
/* READING AND WRITING MAGIC NUMBERS (FILE FORMAT TAGS) */

pnm_format_t pnm_read_format(FILE *rd);
  /* Reads the ``magic number'' of a PBM/PGM/PPM file {rd}. */

void pnm_write_format(FILE *wr, pnm_format_t format);
  /* Writes {format} to {wr} as the ``magic number'' of a PBM/PGM/PPM file. */

/* ********************************************************************** */
/* SKIPPING COMMENTS */

void pnm_skip_whitespace_and_comments(FILE *rd, bool_t verbose);
  /* If the next character in {rd} is whitespace (blank, TAB, CR, or
    LF), consumes that character; else, if the next character is '#',
    consumes it and any additional characters until LF or EOF. Repeats
    these steps until EOF or until the next character is neither
    whitespace nor '#'. Leaves that next character unread. If
    {verbose} is TRUE, prints the '#' comments to {stderr}. */

/* ********************************************************************** */
/* READING AND WRITING IMAGE SIZES */

int32_t pnm_read_plain_int32_t(FILE *rd);
  /* Reads a decimal integer in ASCII from {rd}.  Skips whitespace
    and stops reading at the first non-numeric char.  Bombs out
    if the number is missing (e.g. at EOF) or malformed. */

void pnm_write_plain_int32_t(FILE *wr, int32_t ival);
  /* Writes the integer {ival} in decimal ASCII to {wr}. */

/* ********************************************************************** */
/* READING AND WRITING INDIVIDUAL SAMPLES

  All these procedures complain if the sample value {ival} that
  was read or is to be written is greater than the given {maxval}. */

void pnm_check_sample_range(uint16_t *ival, uint16_t maxval);
  /* If {ival > maxval}, prints a warning and sets {*ival = maxval}.
    Otherwise does nothing. Bombs out after a certain number of warnings. */

uint16_t pnm_read_plain_sample(FILE *rd, uint16_t maxval);
  /* Reads a pixel sample as a signed decimal integer in ASCII from
    {rd}. Skips whitespace and stops reading at the first non-numeric
    char. Bombs out if the number is negative, missing (e.g. at EOF) or
    malformed. */

int32_t pnm_write_plain_sample(FILE *wr, uint16_t ival, uint16_t maxval);
  /* Writes the pixel sample {ival} as an signed decimal integer in ASCII to
    {rd}. Returns the number of digits written. */

uint16_t pnm_read_plain_bit(FILE *rd);
  /* Reads a pixel sample as a single ASCII character '0' or '1' from
    {rd}. Skips whitespace, and stops reading at the first non-numeric
    char. Bombs out if the digit is missing (e.g. at EOF). */

void pnm_write_plain_bit(FILE *wr, uint16_t ival);
  /* Writes the pixel sample {ival} (which must be 0 or 1) 
    as a single ASCII digit '0' or '1' to {rd}. */

uint16_t pnm_read_raw_byte(FILE *rd, uint16_t maxval);
  /* Reads a single byte from {rd}, interprets it as an 8-bit unsigned int32_t.  
    Bombs out at EOF. */

void pnm_write_raw_byte(FILE *wr, uint16_t ival, uint16_t maxval);
  /* Writes the 8-bit unsigned int32_t {ival} as single byte to {wr}. */

uint16_t pnm_read_raw_short(FILE *rd, uint16_t maxval);
  /* Reads two bytes from {rd}, interprets them as a 16-bit unsigned int32_t.  
    Bombs out at EOF. */

void pnm_write_raw_short(FILE *wr, uint16_t ival, uint16_t maxval);
  /* Writes the 16-bit unsigned int32_t {ival} to {wr}, as two bytes,
    big-endian way. */

bool_t pnm_uint_leq(uint32_t x, uint32_t xmax);
  /* Returns TRUE iff {x <= xmax}.  Handy hack to avoid warnings
     about "comparison always true due to limited data range". */

/* ********************************************************************** */
/* READING AND WRITING PIXEL ROWS */

void pnm_read_pixels
  ( FILE *rd, 
    uint16_t *smp, 
    int32_t cols, 
    int32_t chns,
    uint16_t maxval, 
    bool_t raw,
    bool_t bits
  );
  /* Reads a row of pixels from the PNM file {rd} and
    stores them in {smp[0..cols*chns-1]}. Assumes that each sample
    ranges in {0..maxval}, and that the file is positioned right
    before the first sample of the row.
    
    If {raw} is true, assumes the `raw' file format variants
    ({RPBM_FORMAT}, {RPGM_FORMAT}, or {RPPM_FORMAT}), where samples
    are encoded in binary. In that case, assums that each sample is
    encoded as 

      * a single bit (packed 8 to a byte, big-endian) if {bits==TRUE}
      and {maxval==1};

      * as an unsigned byte if {bits==FALSE} and {1<=maxval<=255};
      and

      * as an unsigned big-endian short if {bits==FALSE} and
      {256<=maxval<=65535}.

    If {raw} is false, assumes the `plain' {format} variants
    ({PBM_FORMAT}, {PGM_FORMAT}, or {PPM_FORMAT}), where the samples
    are encoded as ASCII decimal numbers. In this case, {bits} is 
    false, the numbers mus be separated by whitespace and/or newlines; if {bits}
    then {maxval} must be 1, and the whitespaces are optional. 
    
    In any case, if {bits} is true, the bits are complemented after 
    being read from the file. */

void pnm_write_pixels
  ( FILE* wr, 
    uint16_t *smp, 
    int32_t cols, 
    int32_t chns,
    uint16_t maxval, 
    bool_t raw,
    bool_t bits
  );
  /* Writes grayscale pixels {smp[0..cols*chns-1]} to the PGM or PPM file {wr},
    assuming it has the given {maxval} and {format}. Assumes that each
    sample ranges in {0..maxval} and that the previous line
    has just been written.  
    
    If {raw} is true, assumes the `raw' file format variants
    ({RPBM_FORMAT}, {RPGM_FORMAT}, or {RPPM_FORMAT}). In that case,
    writes each sample in binary as 

      * a single bit (packed 8 per byte, big-endian) if {bits==TRUE}
      and {maxval==1};

      * as an unsigned byte if {bits==FALSE} and {1<=maxval<=255},
      and

      * as an unsigned big-endian short if {bits==FALSE} and
      {256<=maxval<=65535}.

    If {raw} is false, assumes the `plain' {format} variants
    ({PGM_FORMAT} or {PPM_FORMAT}), and writes the samples as ASCII
    decimal numbers. In this case, if {bits} is false, each sample is preceded by a
    blank or a newline; if {bits} isthen then {maxval} must be 1, and
    the whitespace may be omitted.
    
    In any case, if {bits} is true, the samples are complemented before
    beng written out. */

/* ********************************************************************** */
/* CONVERSION TO/FROM FLOAT VALUES 

  These procedures convert between integer samples in {0..maxval} 
  and floating-point samples in [0-1].  
  
  Conceptually, the range [0_1] is divided into {N=maxval+1} equal
  intervals. A float sample value {fval} in a given interval gets
  mapped to the interval's index {ival}, counting from 0. Conversely,
  an integer sample value {ival} gets mapped to the center {fval} of
  of the interval with that index.

  The argument {badval} is an optional `invalid' sample value. If
  {badval} is in {0..maxval}, then the integer sample value {badval}
  is mapped to the {double} value {NAN}, which does not
  correspond to any other integer sample value. In that case, the
  value {badval} is implicitly excluded from the integer range, as
  explained below. */

#define PNM_NO_BADVAL ((uint32_t)65536)
  /* Pass this value to the {badval} parameter if all integer sample
    values in {0..maxval} are valid. */

uint16_t pnm_quantize(double fval, uint16_t maxval, bool_t isMask, uint32_t badval);
  /* Converts the floating-point sample value {fval} to a {uint16_t}
    in the range {0..maxval}. 
    
    If {badval} is not in {0..maxval}, the procedure maps any
    {fval} in [0_1) to an integer {ival} in {0..maxval}.
    It also maps any {fval<0} to 0, and {fval>=1.0} to {maxval},
    and {NAN} to {maxval/2}.
    
    If {badval} is in {0..maxval}, the procedure returns {badval}
    whenever {fval} is {NAN}, and effectively excludes the value
    {badval} from the output range. That is, it converts {fval} to an integer 
    in the range {0..maxval-1} instead of {0..maxval} , and adds 1 to
    the computed output sample value {ival} if {ival >= badval}.
    
    In any case, the conversion from {fval} in {[0_1)} to {0..M} is a linear
    scaling followed by rounding. If {iMask} is TRUE, {fval} is mapped
    to {fval*M}, rounded to the nearest integer. Halfway values are
    rounded down if {fval<0.5}, and up if {fval>=0.5}. This choice is
    appropriate e.g. when writing opacity masks.
    
    If {iMask} is FALSE, {fval} is instead mapped to {fval*(M+1)},
    rounded down.  This choice is preferrable when the sample values 0 and {maxval}
    have no special significance, e.g. for properly exposed digital photos. */

double pnm_floatize(uint16_t ival, uint16_t maxval, bool_t isMask, uint32_t badval);
  /* Converts the integer sample value {ival}, assumed to be in {0..maxval},
    to a floating-point number in the range [0_1].  
    
    If {badval} is outside the range {0..maxval}, the integers {0..maxval}
    are converted to float values in the range {[0_1]}. 
    
    If {badval} is in {0..maxval}, the procedure maps the
    sample value {ival==badval} to {NAN}; subtracts 1 from any {ival}
    in {badval..maxval}, and then maps the integers {0..maxval-1}
    to float values in the range {[0_1]}.
    
    In any case, the conversion from {ival} in {0..M} to a float in {[0_1)} is 
    essentially a linear scaling. If {iMask} is TRUE, {ival} is mapped
    to {ival/M}. This choice is appropriate e.g. when writing opacity masks.
    
    If {iMask} is FALSE, {ival} is assumed to represent all real
    numbers in the interval {[ival_ival+1]}, out of a total range of
    {[0_M+1]}; so it is mapped to the float value {(ival+0.5)/(M+1)}.
    This choice is preferrable when the sample values 0 and {maxval}
    have no special significance, e.g. for properly exposed digital
    photos.*/

double *pnm_make_floatize_table(uint16_t maxval, bool_t isMask, uint32_t badval);
  /* Creates a pixel-to-float conversion table {cvt} such that
    {cvt[iv] = floatize(iv, maxval, badval)} for all {iv} in {0..maxval}. */

/* ********************************************************************** */
/* DIAGNOSTIC MESSAGES */

void pnm_message(char* msg, ...);
void pnm_error(char* msg, ...);
  /* Prints the error message {msg} to {stderr}.  
    If additional arguments are given, prints them too, interpreting {msg} as 
    a format spec in the same way as {printf}. The {pnm_error}
    function also aborts the perogram with {exit(1)}. */

/* ********************************************************************** */
/* CHECKED MEMORY ALLOCATION */

#define pnm_malloc(nbytes) \
  notnull(malloc(nbytes), "not enough memory for malloc")
  /* Like {malloc}, but bombs out with a message
    if there is not enough memory. */

#endif
