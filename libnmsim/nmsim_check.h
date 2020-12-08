#ifndef nmsim_check_H
#define nmsim_check_H
 
/* Checking range of values. */
/* Last edited on 2019-03-21 12:45:36 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <nmsim_basic.h>

void nmsim_check_int64_value(char *name, int64_t v, int64_t vmin, int64_t vmax);
  /* Aborts if the value {v} not in the range {vmin..vmax}. */

void nmsim_check_double_value(char *name, double v, double vmin, double vmax);
  /* Aborts if the value {v} is {NAN} is not in the range {[vmin _ vmax]}. */
    
#endif
