/* See {nmsim_read.h} */
/* Last edited on 2019-03-21 12:48:14 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <fget.h>
#include <nget.h>

#include <nmsim_basic.h>
#include <nmsim_check.h>

#include <nmsim_read.h>
 
int64_t nmsim_read_int64_value(FILE *rd, char *name, int64_t vmin, int64_t vmax)
  { int64_t v = fget_int64(rd);
    nmsim_check_int64_value(name, v, vmin, vmax);
    return v;
  }

double nmsim_read_double_value(FILE *rd, char *name, double vmin, double vmax)
  { double v = fget_double(rd);
    if (! isnan(v)) { nmsim_check_double_value(name, v, vmin, vmax); }
    return v;
  }

int64_t nmsim_read_int64_param(FILE *rd, char *name, int64_t vmin, int64_t vmax)
  { int64_t v = nget_int64(rd, name);
    nmsim_check_int64_value(name, v, vmin, vmax);
    fget_eol(rd);
    return v;
  }

double nmsim_read_double_param(FILE *rd, char *name, double vmin, double vmax)
  { double v = nget_double(rd, name);
    if (! isnan(v)) { nmsim_check_double_value(name, v, vmin, vmax); }
    fget_eol(rd);
    return v;
  }
