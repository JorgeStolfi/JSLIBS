/* Input and output of portable 6-dimensional arrays. */
/* Last edited on 2021-06-12 00:58:51 by jstolfi */
/* Copyright © 2005 by Jorge Stolfi, from University of Campinas, Brazil. */

#ifndef ppv_array_read_H
#define ppv_array_read_H

#include <ppv_types.h>
#include <ppv_array.h>
#include <stdio.h>

ppv_array_desc_t *ppv_array_read_file ( FILE *rd, ppv_nbits_t bpw );
  /* Reads an array {A} from file {rd}, assumed to be in the standard
    format (plain or raw) as described in {ppv_array_read_FORMAT_INFO}. 
    
    The procedure automatically allocates the descriptor for {A} and the
    sample storage area {A.el}. The base displacement {A.base} and the
    addressing increments {A.step}, which are not in the file, are
    defined by the procedure based on the bits per sample and the array
    sizes specified in the file. The bits-per-word sample packing
    parameter {A.bpw}, which is also omitted in the file. is set to the
    value {bpw} given (8, 16, or 32).
    
    The procedure aborts the program, with an error message to {stderr},
    if the file is malformed, or there is not enough memory. */ 

#define ppv_array_read_FORMAT_INFO \
  "The external representation of an array {A} consists of a fixed header" \
  " line; the array sizes and the number of bits per sample; the" \
  " values of all samples; and a fixed footer (trailer) line. More" \
  " specifically: \n" \
  "\n" \
  "    * a header line \"begin ppv_array_t (format of 2005-05-28)\".\n" \
  "\n" \
  "    * a line \"dim = {d}\" where {d} is in {0..ppv_MAX_DIM};\n" \
  "\n" \
  "    * a line \"size = {sz[0]} ... {sz[d-1]}\";\n" \
  "\n" \
  "    * a line \"asize = {asz[0]} ... {asz[d-1]}\";\n" \
  "\n" \
  "    * a line \"bps = {M}\" where {M = A.bps};\n" \
  "\n" \
  "    * a line \"plain = {B}\" where {B} is 0 or 1;\n" \
  "\n" \
  "    * the sample values, in a format determined by {B};\n" \
  "\n" \
  "    * a footer line \"end ppv_array_t\".\n" \
  "\n" \
  "  The parameters {sz[0..d-1]} define the nominal shape of the array, i.e. " \
  " {A.size[k] = sz[k]}.\n" \
  "\n" \
  "  The parameters {asz[0..d-1]} are the actual shape of the array, apart from " \
  " virtually replicated indices. For each {k}, either {asz[k] == A.size[k]}, or " \
  " {asz[k] == 1} and {A.size[k] >= 2}. In the second case, {A.step[k]} is zero, " \
  " that is, {A} is virtually replicated along the axis {k}.\n" \
  "\n" \
  "  Note that {d} may be zero, but the tags \"size = \" and \"asize = \" " \
  " are always present.\n" \
  "\n" \
  "  The samples are in C/Pascal order (with the last index in " \
  " the innermost loop), except that virtually replicated indices are " \
  " fixed at 0. Therefore, the number of samples actually written is " \
  " {PROD { asz[i] : i \\in 0..d-1 } }.\n" \
  "\n" \
  "  Two variant formats are supported: /plain/, where the samples are" \
  " stored in free-format ASCII decimal, and /raw/, where their binary" \
  " values are packed or chopped into 8-bit bytes. The plain format" \
  " is provided for the benefit of programmers who need to read" \
  " sample arrays but lack a PPV library or facilities to" \
  " process binary files. The raw format, on the other hand, provides faster" \
  " I/O and much smaller files (from 1/3 to 1/16 of the" \
  " plain equivalent, depending on the number of bits per sample).\n" \
  "\n" \
  "  Specifically, if the plain flag {B} in the file is 1, sample" \
  " values are represented in the /plain/ " \
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
