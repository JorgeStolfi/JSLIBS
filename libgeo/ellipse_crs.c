/* See ellipse.h */
/* Last edited on 2010-03-20 09:42:08 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <assert.h>

#include <r2.h>
#include <r2x2.h>
#include <interval.h>
#include <affirm.h>

#include <ellipse_crs.h>
#include <ellipse_ouv.h>

#define Pr fprintf
#define Er stderr
  
void ellipse_crs_to_ouv(ellipse_crs_t *E, ellipse_ouv_t *F)
  {
    /* Grab the parameters: */
    double rad = fabs(E->rad); /* Just in case... */
    double xstr = E->str.c[0];
    double ystr = E->str.c[1];
    
    /* Compute the main ellipse axes {u,v} and the stretch length {strLen}: */
    double strLen = hypot(xstr, ystr);
    if (strLen == 0.0) 
      { /* Any orthonormal pair will do: */
        F->u = (r2_t) {{ 1.0, 0.0 }};
        F->v = (r2_t) {{ 0.0, 1.0 }};
      }
    else
      { /* Pick {v} at 90 degrees from {u}: */
        F->u = (r2_t) {{ xstr/strLen, ystr/strLen }};
        F->v = (r2_t) {{ - F->u.c[1], +F->u.c[0] }};
      }
    
    /* Compute the minor and major semidiameters {uRad,vRad}: */
    F->a = rad + strLen;
    F->b = rad;
  }

void ellipse_ouv_to_crs(ellipse_ouv_t *F, r2_t *ctr, ellipse_crs_t *E)
  {
    assert(F->a >= F->b);
    E->ctr = *ctr;
    E->rad = F->b;
    r2_scale(F->a - F->b, &(F->u), &(E->str));
  }

void ellipse_crs_bbox(ellipse_crs_t *E, interval_t bbox[])
  {
    ellipse_ouv_t F;
    ellipse_crs_to_ouv(E, &F);
    ellipse_ouv_bbox(&F, bbox);
    LO(bbox[0]) += E->ctr.c[0];
    HI(bbox[0]) += E->ctr.c[0];
    LO(bbox[1]) += E->ctr.c[1];
    HI(bbox[1]) += E->ctr.c[1];
  }

void ellipse_crs_int_bbox
  ( ellipse_crs_t *E, 
    double mrg, /* Extra margin. */
    int *xLoP,  /* (OUT) Min X of clip area. */
    int *xHiP,  /* (OUT) Max X of clip area. */
    int *yLoP,  /* (OUT) Min Y of clip area. */
    int *yHiP   /* (OUT) Max Y of clip area. */
  )
  {
    demand(mrg >= 0, "the extra margin {mrg} must be non-negative");
    interval_t bbox[2];
    ellipse_crs_bbox(E, bbox);
    (*xLoP) = (int)floor(LO(bbox[0]) - mrg);
    (*xHiP) = (int)ceil (HI(bbox[0]) + mrg);
    (*yLoP) = (int)floor(LO(bbox[1]) - mrg);
    (*yHiP) = (int)ceil (HI(bbox[1]) + mrg);
  }

bool_t ellipse_crs_inside(ellipse_crs_t *E, r2_t *p)
  { ellipse_ouv_t F;
    ellipse_crs_to_ouv(E, &F);
    r2_t hp; r2_sub(p, &(E->ctr), &hp);
    return ellipse_ouv_inside(&F, &hp);
  }

r2_t ellipse_crs_relative_coords(ellipse_crs_t *E, r2_t *p)
  { 
    bool_t debug = FALSE;
    ellipse_ouv_t F;
    ellipse_crs_to_ouv(E, &F);
    r2_t op;
    r2_sub(p, &(E->ctr), &op);
    r2_t uvp = (r2_t){{ r2_dot(&op, &(F.u)), r2_dot(&op, &(F.v)) }};
    if (debug) { r2_gen_print(stderr, &uvp, "%12.7f", "uvp abs = (", " ", ")\n"); }
    uvp.c[0] /= F.a;
    uvp.c[1] /= F.b;
    if (debug) { r2_gen_print(stderr, &uvp, "%12.9f", "uvp rel = (", " ", ")\n"); }
    r2_t res;
    r2_mix(uvp.c[0], &(F.u), uvp.c[1], &(F.v), &res);
    if (debug) { r2_gen_print(stderr, &uvp, "%12.9f", "res     = (", " ", ")\n"); }
    return res;
  }

double ellipse_crs_position(ellipse_crs_t *E, r2_t *p, r2_t *csp)
  { ellipse_ouv_t F;
    ellipse_crs_to_ouv(E, &F);
    r2_t hp; r2_sub(p, &(E->ctr), &hp);
    return ellipse_ouv_position(&F, &hp, csp);
  }
  
double ellipse_crs_nearest_point(ellipse_crs_t *E, r2_t *p, r2_t *q)
  { ellipse_ouv_t F;
    ellipse_crs_to_ouv(E, &F);
    r2_t hp; r2_sub(p, &(E->ctr), &hp);
    double dpq = ellipse_ouv_nearest_point(&F, &hp, q);
    if (q != NULL) { r2_add(&(E->ctr), q, q); }
    return dpq;
  }
  
double ellipse_crs_border_position(ellipse_crs_t *E, double hwd, r2_t *p)
  { ellipse_ouv_t F;
    ellipse_crs_to_ouv(E, &F);
    r2_t hp; r2_sub(p, &(E->ctr), &hp);
    return ellipse_ouv_border_position(&F, hwd, &hp);
  }  
  
void ellipse_crs_print(FILE *wr, ellipse_crs_t *E, char *fmt)
  { 
    fprintf(wr, "{ ctr: ");
    fprintf(wr, fmt, E->ctr.c[0]);
    fprintf(wr, " ");
    fprintf(wr, fmt, E->ctr.c[1]);
    fprintf(wr, "  rad: ");
    fprintf(wr, fmt, E->rad);
    fprintf(wr, "  str: ");
    fprintf(wr, fmt, E->str.c[0]);
    fprintf(wr, " ");
    fprintf(wr, fmt, E->str.c[1]);
    fprintf(wr, " }");
    
    fflush(wr);
  }
