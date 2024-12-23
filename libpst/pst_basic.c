/* See pst_basic.h */
/* Last edited on 2024-12-22 12:03:16 by stolfi */

#include <math.h>
#include <assert.h>
#include <values.h>
#include <stdlib.h>
#include <stdint.h>

#include <float_image.h>
#include <vec.h>
#include <r2.h>
#include <r3.h>
#include <argparser.h>

#include <pst_basic.h>

vec_typeimpl(name_vec_t,name_vec,char *);
  
vec_typeimpl(image_vec_t,image_vec,float_image_t *);

void pst_double_vec_regularize(double_vec_t *v, uint32_t NC, double defval)
  { uint32_t KC = v->ne;
    if (KC != NC) 
      { demand(KC <= 1, "range info specified for the wrong number of channels");
        /* Provide default value: */
        double val = (KC == 1 ? v->e[0] : defval);
        double_vec_trim(v, NC);
        for (uint32_t c = 0; c < NC; c++) { v->e[c] = val; }
      }
  }

void pst_double_vec_uniformize(double_vec_t *v, double defval)
  { uint32_t KC = v->ne;
    if (KC == 0) 
      { double_vec_trim(v, 1);
        v->e[0] = defval;
      }
    else if (KC > 1)
      { defval =  v->e[0];
        for (uint32_t c = 1; c < KC; c++) 
          { demand(v->e[c] == defval, "inconsistent values"); }
        double_vec_trim(v, 1);
      }
  }

double_vec_t pst_double_vec_parse(argparser_t *pp, uint32_t *NC)
  { double_vec_t v = double_vec_new(0);
    uint32_t NP = 0; /* Number of values actually parsed. */
    uint32_t NPMAX = ((NC == NULL) || ((*NC) < 0) ? MAXINT : (*NC)); /* Max to parse. */
    while ((NP < NPMAX) && argparser_next_is_number(pp))
      { double_vec_expand(&v, (vec_index_t)NP);
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
        for (uint32_t c = 0; c < NP; c++) { v.e[c] /= den; }
      }
    return v;
  }

void pst_int_vec_regularize(int32_vec_t *v, uint32_t NC, int32_t defval)
  { uint32_t KC = v->ne;
    if ((KC != NC) && (KC <= 1))
      { int32_t val = (KC == 1 ? v->e[0] : defval);
        /* Provide default value: */
        int32_vec_expand(v, (vec_index_t)NC);
        for (uint32_t c = 0; c < NC; c++) { v->e[c] = val; }
        int32_vec_trim(v, NC);
      }
  }

int32_vec_t pst_int_vec_parse(argparser_t *pp, uint32_t *NC)
  { int32_vec_t v = int32_vec_new(0);
    uint32_t NP = 0; /* Number of values actually parsed. */
    uint32_t NPMAX = ((NC == NULL) || ((*NC) < 0) ? MAXINT : (*NC)); /* Max to parse. */
    while ((NP < NPMAX) && argparser_next_is_number(pp))
      { int32_vec_expand(&v, (vec_index_t)NP);
        v.e[NP] = (int32_t)argparser_get_next_int(pp, -MAXINT, +MAXINT);
        NP++;
      }
    if ((NC != NULL) && (NP != 1))
      { /* Set/check {*NC} against {NP}: */
        if ((*NC) < 0)
          { (*NC) = NP; }
        else
          { if (NP != (*NC)) { argparser_error(pp, "wrong number of elements"); } }
      }
    int32_vec_trim(&v, NP);
    return v;
  }
          
