#ifndef wt_table_args_parse_H
#define wt_table_args_parse_H

/* Weight tables for filtering digital signals */
/* Last edited on 2024-11-15 19:17:09 by stolfi */

#define wt_table_args_parse_H_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <argparser.h>
#include <wt_table_generic.h>

wt_table_kind_t wt_table_args_parse_kind(argparser_t *pp);
   /* Parses the next string from the command line, and converts it to 
     a {wt_table_kind_t} using {wt_table_kind_from_string}. */

double_vec_t wt_table_args_parse_weights(argparser_t *pp, bool_t unitNorm);
   /* Parses a sequence of numbers from the command line, in the
     format described by {wt_table_args_parse_weights_HELP} below. \
     
     If {unitNorm} is {TRUE}, scales all weights so that their sum is 1;
     the optional "/ {DENOM}" is parsed but ignored. 
     
     If {unitNorm} is {FALSE}, does no normalization; but, if the
     optional "/ {DENOM}" is present, divides all weights by
     {DENOM}. */
     
#define wt_table_args_parse_weights_HELP \
  "{WEIGHT}.. [ \"/\" {DENOM} ]"

#define wt_table_args_parse_weights_INFO \
  "The weights may be integer or fractional," \
  " and may include \"e\" exponent parts. If the" \
  " optional slash and denominator are present," \
  " all weights are divided by {DENOM}."

#define wt_table_args_parse_weights_norm_sum_INFO \
  "The weights may be integer or fractional," \
  " and may include \"e\" exponent parts. They are" \
  " implicitly normalized to unit sum. (The slash" \
  " and {DENOM}, if present, are parsed but ignored.)"

#endif
