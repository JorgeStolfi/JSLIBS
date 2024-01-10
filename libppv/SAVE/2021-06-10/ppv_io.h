/* Input and output of portable 6-dimensional arrays. */
/* Last edited on 2016-03-10 18:46:10 by stolfilocal */
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

#define ppv_FILE_FORMAT_INFO \
  "The external representation of an array {A} consists of: \n" \
  "\n" \
  "    * a header line \"begin ppv_array_t (format of 2005-05-28)\".\n" \
  "\n" \
  "    * a line \"dim = {N}\" where {N} is in {0..ppv_NAXES};\n" \
  "\n" \
  "    * a line \"size = {sz[0]} ... {sz[N-1]}\";\n" \
  "\n" \
  "    * a line \"asize = {asz[0]} ... {asz[N-1]}\";\n" \
  "\n" \
  "    * a line \"bps = {M}\" where {M = A.bps};\n" \
  "\n" \
  "    * a line \"plain = {B}\" where {B} is 0 or 1;\n" \
  "\n" \
  "    * the sample values, in a format determined by {B};\n" \
  "\n" \
  "    * a footer line \"end ppv_array_t\".\n" \
  "\n" \
  "  The parameters {sz[0..N-1]} define the nominal shape of the array, i.e. " \
  " {A.size[k] = sz[k]}. If {N} is less than {ppv_NAXES}, the input functions " \
  " will automatically extend the index list with {ppv_NAXES-N} trivial indices " \
  " (with {A.size[k] = 1}. Conversely, the output functions will omit any " \
  " trailing indices which have {A.size[k] = 1}.\n" \
  "\n" \
  "  The parameters {asz[0..N-1]} are the actual shape of the array, apart from " \
  " virtually replicated indices. For each {k}, either {asz[k] == A.size[k]}, or " \
  " {asz[k] == 1} and {A.size[k] >= 2}. In the second case, {A.step[k]} is zero, " \
  " that is, {A} is virtually replicated along the axis {k}.\n" \
  "\n" \
  "  Note that {N} may be zero, but the tags \"size = \" and \"asize = \" " \
  " are always present.\n" \
  "\n" \
  "  The samples are written in C/Pascal order (with the last index in " \
  " the innermost loop), except that virtually replicated indices are " \
  " fixed at 0. Therefore, the number of samples actually written is " \
  " {PROD { asz[i] : i \\in 0..N-1 } }.\n" \
  "\n" \
  "  If the flag {B} is 1, sample values are represented in the /plain/ " \
  " format: namely, as decimal numbers with ASCII digits, separated by " \
  " one or more ASCII formatting chars (ASCII 'SPACE', 'TAB', 'NUL', " \
  " 'CR', 'LF', 'FF', or 'VT'). In particular, {write_XXX} procedures of " \
  " this module will insert end-of-line terminators (ASCII 'NL') " \
  " regularly in order to avoid long lines. Extra formatting characters " \
  " may be present before and/or after the sample block.\n" \
  "\n" \
  "  If the flag {B} is 0, sample values are represented in the /raw/ " \
  " (packed binary) format, detailed at the end of this file. Extra " \
  " formatting characters may be present after the sample block, " \
  " but not before it.\n" \
  "\n" \
  "  In either case, the sample block is terminated by an end-of-line " \
  " marker.\n" \
  "\n" \
  "DETAILS OF THE PACKED BINARY FORMAT\n" \
  "\n" \
  "  If the number {bps} of bits per sample is 8 or less, each byte of " \
  " the file will contain {K=8/bps} samples from the high end to the " \
  " low end. For example, the sample block of an array with nine samples " \
  " {A,B,...,I}, with {bps=3}, will have the following layout:\n" \
  "\n" \
  "    | byte#   <--0---> <--1---> <--2---> <--3---> <--4--->   \n" \
  "    | sample#   <0><1>   <2><3>   <4><5>   <6><7>   <8>      \n" \
  "    |         ==AAABBB ==CCCDDD ==EEEFFF ==GGGHHH ==III===   \n" \
  "    | bit#    7      0 7      0 7      0 7      0 7      0   \n" \
  "\n" \
  "  where '=' denotes a padding zero bit.\n " \
  "\n" \
  "  On the other hand, if {bps} is greater than 8, then each sample" \
  " will span {K=(bps+7)/8} bytes, in big-endian order.  If {bps}" \
  " is not a multiple of 8, the first byte of each sample will be" \
  " padded with zeros at the high end.  If" \
  " {bps=18}, for example, the layout will be:\n" \
  "\n" \
  "    | byte#   <--0---> <--1---> <--2---> <--3---> <--4---> <--5---> ... \n" \
  "    | sample#       <-------0---------->       <--------1---------> ... \n" \
  "    |         ======AA AAAAAAAA AAAAAAAA ======BB BBBBBBBB BBBBBBBB ... \n" \
  "    | bit#    7    2 0 7      0 7      0 7    2 0 7      0 7      0 ... " \

#endif
