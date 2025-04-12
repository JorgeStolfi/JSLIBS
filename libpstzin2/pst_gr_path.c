/* See {pst_gr_path.h} */
/* Last edited on 2025-03-14 16:12:13 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <r2.h>
#include <fget.h>

#include <pst_gr_path.h>
 
void pst_gr_path_free(pst_gr_path_t P)
  { free(P.v); }

r2_t pst_gr_path_start_dir(r2_t *org, pst_gr_path_t P, r2_t *dst)
  { r2_t v0;
    if (P.n == 0)
      { v0 = (*dst); }
    else
      { v0 = P.v[0]; }
    r2_t d; r2_sub(&v0, org, &d);
    (void)r2_dir(&d, &d);
    return d;
  }

pst_gr_path_t pst_gr_path_reverse(pst_gr_path_t P)
  { P.reverse = !P.reverse;
    return P;
  }

pst_gr_path_t pst_gr_path_concatenate(pst_gr_path_t P0, r2_t *mid, pst_gr_path_t P1)
  { uint32_t n0 = P0.n;
    uint32_t n1 = P1.n;
    uint32_t n = n0 + n1 + 1;

    r2_t *v = (r2_t*)malloc(sizeof(r2_t)*n);
    for (uint32_t i = 0; i < n; i++)
      { r2_t c;
        if (i < n0)
          { c = P0.v[i]; }
        else if (i == n0)
          { c = (*mid); }
        else
          { uint32_t j = (uint32_t)(i - n0 - 1);
            c = P1.v[j];
          }
        v[i] = c;
      }
    return (pst_gr_path_t) { .n=n, .v=v, .reverse=FALSE };
  }

void pst_gr_path_write(FILE *wr, pst_gr_path_t P)
  {
    fprintf(wr,"%d %d", P.n, P.reverse);
    for (uint32_t i = 0; i < P.n; i++)
      { fprintf(wr," (%9.6f,%9.6f)",P.v[i].c[0],P.v[i].c[1]); }
  }

pst_gr_path_t pst_gr_path_read(FILE *rd)
  {
    pst_gr_path_t P = pst_gr_path_NULL;
    P.n = fget_uint32(rd, 10);
    demand(P.n < pst_gr_path_MAX_VERTICES, "too many vertices in path");
    
    uint32_t rev = fget_uint32(rd, 10);
    demand((rev == 0) || (rev == 1), "invalid path reversal bit");
    P.reverse = (rev == 1);
    if (P.n == 0) 
      { P.v = NULL; }
    else
      { P.v = talloc(P.n, r2_t);
        for (uint32_t i = 0; i < P.n; i++)
          { fget_skip_spaces_and_match(rd, "(");
            double X = fget_double(rd);
            fget_skip_spaces_and_match(rd, ",");
            double Y = fget_double(rd);
            fget_skip_spaces_and_match(rd, ")");
            demand(isfinite(X) && isfinite(Y), "invalid path coordinate");
            P.v[i] = (r2_t){{ X, Y }};
          }
      }
    return P;
  }

void pst_gr_path_ctr_dir(r2_t *org, pst_gr_path_t P, r2_t *dst, r2_t *ctr_P, r2_t *dir_P)
  { r2_t ctr, dir;
    uint32_t h = P.n/2;
    if (P.n == 0)
      { r2_mix(0.5, org, 0.5, dst, &ctr); 
        r2_sub(dst, org, &dir);
      }
    else if (P.n == 1)
      { ctr = P.v[h];
        r2_sub(dst, org, &dir);
      }
    else if ((P.n % 2) == 1)
      { ctr = P.v[h];
        r2_t *pm = &(P.v[h-1]); 
        r2_t *pp = &(P.v[h + 1]); 
        r2_sub(pp, pm, &dir); 
      }
    else
      { r2_t *p0 = &(P.v[h]); 
        r2_t *p1 = &(P.v[h + 1]); 
        r2_mix(0.5, p0, 0.5, p1, &ctr);
        r2_sub(p1, p0, &dir);
      }
    (*ctr_P) = ctr;
    (*dir_P) = dir;
  }
