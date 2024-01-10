#ifndef msm_basic_H
#define msm_basic_H

/* Basic types for multiscale matching of univariate signals. */
/* Last edited on 2008-01-11 12:15:40 by stolfi */

#define msm_basic_H_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <bool.h>
#include <sign.h>
#include <values.h>
#include <math.h>

#define INF INFINITY

void msm_check_malloc(char *tag);
  /* A procedure that calls {malloc} and {free} several times, in order
    to test the integrity of the {malloc} heap. */

int msm_choose(int lo, int hi, double dProb);
  /* Generates a random integer in the range {lo..hi}.
    If {dProb} is zero, the result is the center {(lo+hi)/2}
    of the range, rounded randomly. As {dProb} increases,
    the distribution widens.  If {dprob} is 1.0, the 
    result is either {lo} or {hi}, each with probebility 1/2. */

FILE *msm_open_read(char *name, char *tag, char *ext, bool_t verbose);
FILE *msm_open_write(char *name, char *tag, char *ext, bool_t verbose);
  /* These procedure open a disk file called "{name}{tag}{ext}",
    for reading or writing, respectively.  A message is printed
    to {stderr} if {verbose} is true.  The procedures abort on
    any errors. */

#endif
