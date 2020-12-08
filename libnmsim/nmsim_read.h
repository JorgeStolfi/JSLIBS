#ifndef nmsim_read_H
#define nmsim_read_H
 
/* Value reading procedures. */
/* Last edited on 2019-03-21 12:46:40 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <nmsim_basic.h>

int64_t nmsim_read_int64_value(FILE *rd, char *name, int64_t vmin, int64_t vmax);
  /* Reads an integer from {rd}, which must be on the current line. 
    Aborts if the value read is in the range {vmin..vmax}.  */

double nmsim_read_double_value(FILE *rd, char *name, double vmin, double vmax);
  /* Reads a {double} value from {rd}, which must be on the current line. 
    Aborts if the value read is not not {NAN} and not in the range {[vmin _ vmax]}.  */

int64_t nmsim_read_int64_param(FILE *rd, char *name, int64_t vmin, int64_t vmax);
  /* Reads a line \"{name} = {value}\", including the end-of-line.
    Aborts if the {value} is not an integer in the range {vmin..vmax}..  */

double nmsim_read_double_param(FILE *rd, char *name, double vmin, double vmax);
  /* Reads a line \"{name} = {value}\", including the end-of-line.
    Aborts if the {value} is not {NAN} and not in the range {[vmin _ vmax]}. */

#endif
