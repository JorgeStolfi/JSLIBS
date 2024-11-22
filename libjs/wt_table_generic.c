/* See wt_table_generic.h */
/* Last edited on 2024-11-16 10:27:27 by stolfi */

#define wt_table_generic_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <values.h>

#include <affirm.h>
#include <jsstring.h>
#include <argparser.h>
#include <gauss_distr.h>

#include <wt_table.h>
#include <wt_table_gaussian.h>
#include <wt_table_triangular.h>
#include <wt_table_uniform.h>
#include <wt_table_binomial.h>
#include <wt_table_hann.h>

#include <wt_table_generic.h>

wt_table_kind_t wt_table_kind_from_string(char *name)
  { demand(strlen(name) > 0, "invalid empty {name}");
    if (isprefix(name, "gaussian"))
      { return wt_table_kind_GAUSSIAN; }
    else if (isprefix(name, "binomial"))
      { return wt_table_kind_BINOMIAL; }
    else if (isprefix(name, "triangular"))
      { return wt_table_kind_TRIANGULAR; }
    else if (isprefix(name, "hann"))
      { return wt_table_kind_HANN; }
    else if (isprefix(name, "uniform"))
      { return wt_table_kind_UNIFORM; }
    else 
      { return wt_table_kind_INVALID; }
  }

char *wt_table_kind_to_string(wt_table_kind_t kind)
  {
    switch(kind)
      { case wt_table_kind_GAUSSIAN:   return "gaussian";
        case wt_table_kind_HANN:       return "hann";       
        case wt_table_kind_BINOMIAL:   return "binomial";   
        case wt_table_kind_TRIANGULAR: return "triangular"; 
        case wt_table_kind_UNIFORM:    return "uniform";
        case wt_table_kind_INVALID:    return "invalid";
        default: demand(FALSE, "invalid {kind}");
      }
  }

void wt_table_fill(wt_table_kind_t kind, uint32_t n, double parm, double wt[], bool_t norm, uint32_t *stride_P)
  {
    switch(kind)
      { case wt_table_kind_GAUSSIAN:   wt_table_gaussian_fill(n, parm, wt, stride_P); break;
        case wt_table_kind_HANN:       wt_table_hann_fill(n, parm, wt, stride_P); break;       
        case wt_table_kind_BINOMIAL:   wt_table_binomial_fill(n, wt, stride_P); break;   
        case wt_table_kind_TRIANGULAR: wt_table_triangular_fill(n, wt, stride_P); break; 
        case wt_table_kind_UNIFORM:    wt_table_uniform_fill(n, parm, wt, stride_P); break;
        default: demand(FALSE, "invalid {kind}");
      }
    if (norm) { wt_table_normalize_sum(n, wt); }
  }
  
double_vec_t wt_table_make(wt_table_kind_t kind, uint32_t n, double parm, bool_t norm, uint32_t *stride_P)
  { double_vec_t wt = double_vec_new(n);
    wt_table_fill(kind, n, parm, wt.e, norm, stride_P);
    return wt;
  }
