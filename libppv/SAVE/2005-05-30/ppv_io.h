/* Input and output of portable 6-dimensional arrays. */
/* Last edited on 2005-05-29 10:46:25 by stolfi */
/* Copyright © 2005 by Jorge Stolfi, from University of Campinas, Brazil. */
/* See the rights and conditions notice at the end of this file. */

#ifndef ppv_io_H
#define ppv_io_H

#include <ppv_types.h>
#include <ppv_array.h>
#include <stdio.h>

/* 
  FILE FORMAT
  
    This module provides functions that read and write PPV arrays from
    open file descriptors, in a standard format. The external
    representation consists of a fixed header line; the array sizes
    and the number of bits per sample; the values of all samples; and
    a fixed footer (trailer) line.  The details of the format are specified at the
    end of this file.
    
    Two variant formats are supported: /plain/, where the samples are
    stored in free-format ASCII decimal, and /raw/, where their binary
    values are packed or chopped into 8-bit bytes. The plain format is
    provided for the benefit of programmers who need to read sample
    arrays but lack a PPV library or facilities to process binary
    files. The raw format provides faster I/O and much smaller files
    (from 1/3 to 1/16 of the plain equivalent, depending on the number
    of bits per sample). */
  
/* 
   READ/WRITE OF WHOLE ARRAYS */

void ppv_write_array(FILE *wr, ppv_array_t *A, bool_t plain);
  /* Outputs the array {A} to file {wr}, with samples in the 
    format specified by {plain}. */

ppv_array_t ppv_read_array(FILE *rd, ppv_nbits_t bpw);
  /* Reads an array {A} from file {rd}, assumed to be in the standard
    format (plain or raw). 
    
    Automatically allocates the sample array {A.el} and computes the
    steps {A.step}. Aborts the program, with an error message to
    {stderr}, if the file is malformed, or there is not enough
    memory. */ 

/*
   FILE FORMAT

    The external representation of an array {A} consists of:

      * a header line "begin ppv_array_t (format of 2005-05-28)".

      * a line "axes = {N}" where {N} does not exceed {ppv_NAX};

      * a line "size = {S[0]} ... {S[N-1]}" where {S[i] = A.size[i]};
      
      * a line "repeat = {R[0]} ... {R[N-1]}" where {R[i]} is either 0 or 1;

      * a line "bps = {M}" where {M = A.bps};

      * a line "plain = {B}" where {B} is 0 or 1;

      * the sample values, in a format determined by {B};

      * a footer line "end ppv_array_t".

    If {N} is less than {ppv_NAX}, the input functions automatically
    extend the index list with {ppv_NAX-N} trivial indices (with 
    {size[k] = 1}. Conversely, the output functions will omit
    any trailing indices which have {size[k] = 1}. 
    
    When {R[k]} is 1, the corresponding index is assumed to be 
    virtually replicated.
    
    If the flag {B} is 1, sample values are represented in the /plain/
    format: namely, as decimal numbers with ASCII digits, separated by
    one or more ASCII formatting chars (ASCII 'SPACE', 'TAB', 'NUL',
    'CR', 'LF', 'FF', or 'VT'). In particular, {write_XXX} procedures of
    this module will insert end-of-line terminators (ASCII 'NL')
    regularly in order to avoid long lines. Extra formatting characters
    may be present before and/or after the sample block.  

    If the flag {B} is 0, sample values are represented in the /raw/
    (packed binary) format, detailed at the end of this file. Extra
    formatting characters may be present after the sample block,
    but not before it.

    In either case, the sample block is terminated by an end-of-line
    marker.
  
  DETAILS OF THE PACKED BINARY FORMAT

    If the number {bps} of bits per sample is 8 or less, each byte of
    the file will contain {K=8/bps} samples from the high end to the
    low end. For example, the sample block of an array with nine samples
    {A,B,...,I}, with {bps=3}, will have the following layout

      | byte#   <--0---> <--1---> <--2---> <--3---> <--4--->
      | sample#   <-0-->   <-1-->   <-2-->   <-3-->   <-4-->      
      |         ==AAABBB ==CCCDDD ==EEEFFF ==GGGHHH ==III===
      | bit#    7      0 7      0 7      0 7      0 7      0   

    where '=' denotes a padding zero bit. 

    On the other hand, if {bps} is greater than 8, then each sample
    will span {K=(bps+7)/8} bytes, in big-endian order.  If {bps}
    is not a multiple of 8, the first byte of each sample will be
    padded with zeros at the high end.  If
    {bps=18}, for example, the layout will be

      | byte#   <--0---> <--1---> <--2---> <--3---> <--4---> <--5---> ...
      | sample#       <-------0---------->       <--------1---------> ...
      |         ======AA AAAAAAAA AAAAAAAA ======BB BBBBBBBB BBBBBBBB ...
      | bit#    7    2 0 7      0 7      0 7    2 0 7      0 7      0 ... 

*/
  
#endif
