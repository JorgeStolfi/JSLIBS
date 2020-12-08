/* See pst_basic.h */
/* Last edited on 2016-04-01 01:48:56 by stolfilocal */

#define _GNU_SOURCE
#include <math.h>
#include <assert.h>
#include <values.h>
#include <stdlib.h>

#include <float_image.h>
#include <vec.h>
#include <r2.h>
#include <r3.h>
#include <argparser.h>

#include <pst_basic.h>

vec_typeimpl(name_vec_t,name_vec,char *);
  
vec_typeimpl(image_vec_t,image_vec,float_image_t *);

void pst_double_vec_regularize(double_vec_t *v, int NC, double defval)
  { int KC = v->ne;
    if (KC != NC) 
      { demand(KC <= 1, "range info specified for the wrong number of channels");
        /* Provide default value: */
        double val = (KC == 1 ? v->e[0] : defval);
        double_vec_trim(v, NC);
        int c;
        for (c = 0; c < NC; c++) { v->e[c] = val; }
      }
  }

void pst_double_vec_uniformize(double_vec_t *v, double defval)
  { int KC = v->ne;
    if (KC == 0) 
      { double_vec_trim(v, 1);
        v->e[0] = defval;
      }
    else if (KC > 1)
      { defval =  v->e[0];
        int c;
        for (c = 1; c < KC; c++) 
          { demand(v->e[c] == defval, "inconsistent values"); }
        double_vec_trim(v, 1);
      }
  }

double_vec_t pst_double_vec_parse(argparser_t *pp, int *NC)
  { double_vec_t v = double_vec_new(0);
    int NP = 0; /* Number of values actually parsed. */
    int NPMAX = ((NC == NULL) || ((*NC) < 0) ? MAXINT : (*NC)); /* Max to parse. */
    while ((NP < NPMAX) && argparser_next_is_number(pp))
      { double_vec_expand(&v, NP);
        v.e[NP] = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX);
        NP++;
      }
    if ((NC != NULL) && (NP != 1))
      { /* Set/check {*NC} against {NP}: */
        if ((*NC) < 0)
          { (*NC) = NP; }
        else
          { if (NP != (*NC)) { argparser_error(pp, "wrong number of elements"); } }
      }
    double_vec_trim(&v, NP);
    if (argparser_keyword_present_next(pp, "/"))
      { double den = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX);
        if (den == 0.0) { argparser_error(pp, "bad denominator"); } 
        int c;
        for (c = 0; c < NP; c++) { v.e[c] /= den; }
      }
    return v;
  }

void pst_int_vec_regularize(int_vec_t *v, int NC, int defval)
  { int KC = v->ne;
    if ((KC != NC) && (KC <= 1))
      { int val = (KC == 1 ? v->e[0] : defval);
        /* Provide default value: */
        int_vec_expand(v, NC);
        int c;
        for (c = 0; c < NC; c++) { v->e[c] = val; }
        int_vec_trim(v, NC);
      }
  }

int_vec_t pst_int_vec_parse(argparser_t *pp, int *NC)
  { int_vec_t v = int_vec_new(0);
    int NP = 0; /* Number of values actually parsed. */
    int NPMAX = ((NC == NULL) || ((*NC) < 0) ? MAXINT : (*NC)); /* Max to parse. */
    while ((NP < NPMAX) && argparser_next_is_number(pp))
      { int_vec_expand(&v, NP);
        v.e[NP] = (int)argparser_get_next_int(pp, -MAXINT, +MAXINT);
        NP++;
      }
    if ((NC != NULL) && (NP != 1))
      { /* Set/check {*NC} against {NP}: */
        if ((*NC) < 0)
          { (*NC) = NP; }
        else
          { if (NP != (*NC)) { argparser_error(pp, "wrong number of elements"); } }
      }
    int_vec_trim(&v, NP);
    return v;
  }
          
