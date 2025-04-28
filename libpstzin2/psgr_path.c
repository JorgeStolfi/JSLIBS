/* See {psgr_path.h} */
/* Last edited on 2025-04-27 11:30:16 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <jsrandom.h>
#include <i2.h>
#include <r2.h>
#include <fget.h>
#include <bool.h>

#include <psgr_path.h>
 
psgr_path_t psgr_path_copy(psgr_path_t P)
  { psgr_path_t P_copy = (psgr_path_t){ .n = P.n, .v = NULL, .reverse = P.reverse }; ;
    if (P.n == 0)
      { assert(P.v == NULL); }
    else
      { P_copy.v = talloc(P.n, r2_t);
        for (int32_t i = 0; i < P.n; i++) { P_copy.v[i] = P.v[i]; }
      }
    return P_copy;
  }

bool_t psgr_path_equal(psgr_path_t P0, psgr_path_t P1)
  {    
    if (P1.n != P0.n) { return FALSE; }
    if (P1.reverse != P0.reverse) { return FALSE; }
    if (P0.n == 0)
      { assert(P0.v == NULL); assert(P1.v == NULL); }
    else
      { for (uint32_t k = 0; k < P0.n; k++)
          { if (P0.v[k].c[0] != P1.v[k].c[0])  { return FALSE; }
            if (P0.v[k].c[1] != P1.v[k].c[1])  { return FALSE; }
          }
      }
    return TRUE;
  }


void psgr_path_free(psgr_path_t P)
  { free(P.v); }

r2_t psgr_path_starting_dir(i2_t *ipo, psgr_path_t P, i2_t *ipd)
  { r2_t p0 = (r2_t){{ ipo->c[0], ipo->c[1] }};
    r2_t p1;
    if (P.n == 0)
      { p1 = (r2_t){{ ipd->c[0], ipd->c[1] }}; }
    else
      { p1 = P.v[0]; }
    r2_t d; r2_sub(&p1, &p0, &d);
    (void)r2_dir(&d, &d);
    return d;
  }

psgr_path_t psgr_path_reverse(psgr_path_t P)
  { P.reverse = !P.reverse;
    return P;
  }

psgr_path_t psgr_path_concatenate(psgr_path_t P0, r2_t *pm, psgr_path_t P1)
  { uint32_t n0 = P0.n;
    uint32_t n1 = P1.n;
    uint32_t n = n0 + n1 + 1;

    r2_t *v = talloc(n, r2_t);
    for (uint32_t i = 0; i < n; i++)
      { r2_t c;
        if (i < n0)
          { uint32_t i0 = i; 
            c = (P0.reverse ? P0.v[n0-1-i0] : P0.v[i]);
          }
        else if (i == n0)
          { c = *pm; }
        else
          { uint32_t i1 = (uint32_t)(i - n0 - 1);
            c = (P1.reverse ? P1.v[n1-1-i1] : P1.v[i1]);
          }
        v[i] = c;
      }
    return (psgr_path_t) { .n=n, .v=v, .reverse=FALSE };
  }

void psgr_path_write(FILE *wr, psgr_path_t P)
  {
    fprintf(wr,"%d %d", P.n, P.reverse);
    for (uint32_t i = 0; i < P.n; i++)
      { fprintf(wr," (%5.3f,%5.3f)", P.v[i].c[0], P.v[i].c[1]); }
  }

psgr_path_t psgr_path_read(FILE *rd, uint32_t NX, uint32_t NY)
  {
    psgr_path_t P = psgr_path_NULL;
    P.n = fget_uint32(rd, 10);
    demand(P.n < psgr_path_MAX_NODES, "too many nodes in path");
    
    uint32_t rev = fget_uint32(rd, 10);
    demand((rev == 0) || (rev == 1), "invalid path reversal bit");
    P.reverse = (rev == 1);
    if (P.n == 0) 
      { P.v = NULL; }
    else
      { P.v = talloc(P.n, r2_t);
        for (uint32_t i = 0; i < P.n; i++)
          { fget_skip_spaces_and_match(rd, "(");
            double x = fget_double(rd);
            fget_skip_spaces_and_match(rd, ",");
            double y = fget_double(rd);
            fget_skip_spaces_and_match(rd, ")");
            demand((x >= 0) && (x <= NX-1) && (y >= 0) && (y <= NY-1), "invalid path node coords");
            P.v[i] = (r2_t){{ x, y }};
          }
      }
    return P;
  }

r2_t psgr_path_center(i2_t *ipo, psgr_path_t P, i2_t *ipd)
  { r2_t ctr;
  
    uint32_t h = P.n/2;
    if (P.n == 0)
      { r2_t po = (r2_t){{ ipo->c[0], ipo->c[1] }};
        r2_t pd = (r2_t){{ ipd->c[0], ipd->c[1] }};
        r2_mix(0.5, &po, 0.5, &pd, &ctr); 
      }
    else if ((P.n % 2) == 1)
      { ctr = P.v[h]; }
    else
      { uint32_t i = h+1;
        r2_t ph = P.v[h]; 
        r2_t pi = P.v[i];
        r2_mix(0.5, &ph, 0.5, &pi, &ctr);
      }
    return ctr;
  }

r2_t psgr_path_central_dir(i2_t *ipo, psgr_path_t P, i2_t *ipd)
  { r2_t dir;
    uint32_t h = P.n/2;
    if ((P.n == 0) || (P.n == 1))
      { r2_t po = (r2_t){{ ipo->c[0], ipo->c[1] }};
        r2_t pd = (r2_t){{ ipd->c[0], ipd->c[1] }};
        r2_sub(&pd, &po, &dir);
        /* Not affected by {P.reverse}. */
      }
    else if ((P.n % 2) == 1)
      { uint32_t i = h-1;
        uint32_t j = h+1;
        r2_sub(&(P.v[j]), &(P.v[i]), &dir); 
        if (P.reverse) { r2_neg(&dir, &dir); }
      }
    else
      { uint32_t j = h+1;
        r2_sub(&(P.v[j]), &(P.v[h]), &dir);
        if (P.reverse) { r2_neg(&dir, &dir); }
      }
    (void)r2_dir(&dir, &dir);
    return dir;
  }

psgr_path_t psgr_path_in_star_wedge(i2_t *ip1, i2_t *ipm, i2_t *ip2)
  { r2_t pm = (r2_t){{ ipm->c[0], ipm->c[1] }};
    r2_t p1 = (r2_t){{ ip1->c[0], ip1->c[1] }};
    r2_t p2 = (r2_t){{ ip2->c[0], ip2->c[1] }};
    r2_t bar =
      (r2_t)
        {{  (pm.c[0] + p1.c[0] + p2.c[0] )/3.0, 
            (pm.c[1] + p1.c[1] + p2.c[1] )/3.0
        }};
        
    psgr_path_t P = (psgr_path_t){ .n = 1, .v = talloc(1, r2_t), .reverse = FALSE };
    P.v[0] = bar;
    
    return P;
  }

psgr_path_t psgr_path_throw(i2_t *ip0, uint32_t n, i2_t *ip1)
  { r2_t p0 = (r2_t){{ ip0->c[0], ip0->c[1] }};
    r2_t p1 = (r2_t){{ ip1->c[0], ip1->c[1] }};
    psgr_path_t P = psgr_path_NULL;

    if (n > 0)
      { P.n = n;
        P.v = talloc(n, r2_t);
        /* Get a transversal vector {u}: */
        r2_t u; r2_sub(&p1, &p0, &u);
        r2_cross(&u, &u); 
        for (int32_t k = 0; k < n; k++)
          { double r = ((double)k+1)/(n+1); /* Position along segment {ip0--ip1}. */
            r2_t pk; r2_mix(1-r, &p0, r, &p1, &pk);
            double t = (n == 1 ? 0.5 : ((k > 0) && (k < n-1) ? ((double)k)/(n-1) : 0.0)); /* Deviation. */
            double s = t*(1-t)*dabrandom(-1.0,+1.0);
            r2_mix(1.0, &pk, s, &u, &pk);
            P.v[k] = pk;
          }
      }
    return P;
  }
