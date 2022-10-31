/* See {dnae_sample.h}. */
/* Last edited on 2022-10-31 09:41:34 by stolfi */

#define dnae_sample_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>
#include <assert.h>

#include <vec.h>
#include <bool.h>
#include <affirm.h>

#include <dnae_nucleic.h>
#include <dnae_seq.h>

#define dnae_sample_enc_VERY_MIN (-32768)
#define dnae_sample_enc_VERY_MAX (+32767)
  /* True range of any {dnae_sample_t} variable, including 
    invalid values. */

#define dnae_sample_BIAS (-dnae_sample_enc_VERY_MIN)
  /* An integer such that {v + dnae_sample_BIAS >= 0} for any 
    {dnae_sample_t} {v}, valid or not. */

/* INTERNAL PROTOTYPES: */

static double *dnae_decode_table = NULL;
  /* Sample decoding table.  The decoded value of signed byte {v}
    is {dnae_decode_table[dnae_sample_BIAS+v]}. */

void dnae_decode_table_setup(void); 
  /* Allocates and initializes the table {dnae_decode_table}. */

/* IMPLEMENTATIONS: */

void dnae_decode_table_setup(void)
  { 
    auto double cumpr(double x);
      /* Cumulative probability distribution for one-sided
        (0,1)-normal variable. */
      
    double cumpr(double x) { return erf(x/M_SQRT2); }
    
    /* Paranoia: */
    assert(dnae_sample_enc_VERY_MIN <= dnae_sample_enc_VALID_MIN);
    assert(dnae_sample_enc_VERY_MAX >= dnae_sample_enc_VALID_MAX);
    assert(dnae_sample_enc_VALID_MIN + dnae_sample_enc_VALID_MAX == 0);
    
    int32_t N = dnae_sample_BIAS + dnae_sample_enc_VERY_MAX + 1;
    double *tb = (double *)notnull(malloc(N*sizeof(double)), "no mem");
    int32_t v; /* Must be {int32_t} and not {dnae_sample_t} to avoid overflow. */
    /* Fill the table with {NaN}s: */
    for (v = dnae_sample_enc_VERY_MIN; v <= dnae_sample_enc_VERY_MAX; v++)
      { tb[dnae_sample_BIAS + v] = NAN; }
    /* Now set valid entries: */
    tb[dnae_sample_BIAS + 0] = 0.0;
    double xprev = 0.0;
    for (v = 1; v <= dnae_sample_enc_VALID_MAX; v++)
      { /* Map {v} linearly to a probability between 0 and 1: */
        double fv = ((double)v)/(dnae_sample_enc_VALID_MAX + 0.5);
        /* Find {x} such that {erf(x/sqrt(2)) = fv}: */
        double xinf = xprev;
        double xsup = 10.0;
        /* fprintf(stderr, "\n"); */
        /* fprintf(stderr, "looking for fv = %24.16e\n", fv); */
        /* fprintf(stderr, "xinf = %24.16e cumpr(xinf) = %24.16e\n", xinf, cumpr(xinf)); */
        /* fprintf(stderr, "xsup = %24.16e cumpr(xsup) = %24.16e\n", xsup, cumpr(xsup)); */
        double x;
        assert(cumpr(xinf) <= fv);
        assert(cumpr(xsup) >= fv);
        while(TRUE)
          { double xrad = (xsup - xinf)/2;
            x = xinf + xrad;
            x = fabs(x); /* Hack to ensure that {x} is truncated to double instead of extended. */
            /* fprintf(stderr, "%24.16e %24.16e %24.16e\n", xinf, x, xsup); */
            if ((x <= xinf) || (x >= xsup)) { break; }
            double ex = cumpr(x);
            if (fv < ex) 
              { xsup = x; }
            else if (fv > ex)
              { xinf = x; }
            else
              { break; }
          }
        /* fprintf(stderr, "%03d = %24.16e\n", v, x); */
        /* Store in table: */
        tb[dnae_sample_BIAS + v] = +x;
        tb[dnae_sample_BIAS - v] = -x;
        xprev = x;
      }
    dnae_decode_table = tb;
  }       

double dnae_sample_decode(dnae_sample_enc_t ev, double scale)
  { int32_t iev = ev;
    demand((iev >= dnae_sample_enc_VALID_MIN) && (iev <= dnae_sample_enc_VALID_MAX), "invalid encoded value"); 
    if (dnae_decode_table == NULL) { dnae_decode_table_setup(); }
    int32_t k = ev + dnae_sample_BIAS; /* Should be in {1..2^16-1}. */
    return scale*dnae_decode_table[k];
  }

dnae_sample_enc_t dnae_sample_encode(double dv, double scale)
  { /* Must check whether this is the best rounding: */
    double fv = erf((dv/scale)/M_SQRT2); /* Result in {[-1 _ +1]}. */
    int32_t v = (int32_t)round(fv*(dnae_sample_enc_VALID_MAX + 0.5));
    if (v > dnae_sample_enc_VALID_MAX) { v = dnae_sample_enc_VALID_MAX; }
    if (v < dnae_sample_enc_VALID_MIN) { v = dnae_sample_enc_VALID_MIN; }
    return (dnae_sample_enc_t)v;
  }

double dnae_sample_diffsq(dnae_sample_enc_t xs, double xscale, dnae_sample_enc_t ys, double yscale)
  { double D = dnae_sample_decode(xs, xscale) - dnae_sample_decode(ys, yscale);
    return D*D;
  }

vec_typeimpl(dnae_sample_enc_vec_t,dnae_sample_enc_vec,dnae_sample_enc_t);

