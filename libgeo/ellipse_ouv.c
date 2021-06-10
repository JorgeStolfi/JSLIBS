/* See ellipse_ouv.h */
/* Last edited on 2021-06-09 19:38:42 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <assert.h>

#include <r2.h>
#include <r2x2.h>
#include <interval.h>
#include <affirm.h>

#include <ellipse_ouv.h>
#include <ellipse_aligned.h>

#define Pr fprintf
#define Er stderr
  
void ellipse_ouv_bbox(ellipse_ouv_t *F, interval_t bbox[])
  {
    /* Roundoff safety epsilon: */
    double eps = 1.0e-8*(fabs(F->a) + fabs(F->b));
    
    double dx = hypot(F->a*F->u.c[0], F->b*F->v.c[0]);
    double dy = hypot(F->a*F->u.c[1], F->b*F->v.c[1]);

    LO(bbox[0]) = - dx - eps;
    HI(bbox[0]) = + dx + eps;
    LO(bbox[1]) = - dy - eps;
    HI(bbox[1]) = + dy + eps;
  }

void ellipse_ouv_int_bbox
  ( r2_t *ctr,  /* Center of ellipse. */
    ellipse_ouv_t *F, 
    double mrg, /* Extra margin. */
    int32_t *xLoP,  /* (OUT) Min X of clip area. */
    int32_t *xHiP,  /* (OUT) Max X of clip area. */
    int32_t *yLoP,  /* (OUT) Min Y of clip area. */
    int32_t *yHiP   /* (OUT) Max Y of clip area. */
  )
  {
    demand(mrg >= 0, "the extra margin {mrg} must be non-negative");
    interval_t bbox[2];
    ellipse_ouv_bbox(F, bbox);
    double cx = ctr->c[0];
    double cy = ctr->c[1];
    (*xLoP) = (int32_t)floor(LO(bbox[0]) + cx - mrg);
    (*xHiP) = (int32_t)ceil (HI(bbox[0]) + cx + mrg);
    (*yLoP) = (int32_t)floor(LO(bbox[1]) + cy - mrg);
    (*yHiP) = (int32_t)ceil (HI(bbox[1]) + cy + mrg);
  }

double ellipse_ouv_position(ellipse_ouv_t *F, r2_t *p, r2_t *csp)
  { r2_t uvp = (r2_t){{ r2_dot(p, &(F->u)), r2_dot(p, &(F->v)) }};
    double *cosP = (csp == NULL ? NULL : &(csp->c[0]));
    double *sinP = (csp == NULL ? NULL : &(csp->c[1]));
    return ellipse_aligned_position(F->a, F->b, uvp.c[0], uvp.c[1], cosP, sinP);
  }
  
double ellipse_ouv_nearest_point(ellipse_ouv_t *F, r2_t *p, r2_t *q)
  { r2_t uvp = (r2_t){{ r2_dot(p, &(F->u)), r2_dot(p, &(F->v)) }};
    r2_t uvq;
    double dpq = ellipse_aligned_nearest_point(F->a, F->b, uvp.c[0], uvp.c[1], &(uvq.c[0]), &(uvq.c[1]));
    if (q != NULL) { r2_mix(uvq.c[0], &(F->u), uvq.c[1], &(F->v), q); }
    return dpq;
  }
  
bool_t ellipse_ouv_inside(ellipse_ouv_t *F, r2_t *p)
  { /* Quick test against square with side {2*a}: */
    if ((fabs(p->c[0]) > F->a) || (fabs(p->c[0]) > F->a)) { return FALSE; }
    /* Convert to UV coordinates and test as aligned: */
    r2_t uvp = (r2_t){{ r2_dot(p, &(F->u)), r2_dot(p, &(F->v)) }};
    return ellipse_aligned_inside(F->a, F->b, uvp.c[0], uvp.c[1]);
  }

double ellipse_ouv_border_position(ellipse_ouv_t *F, double hwd, r2_t *p)
  { /* Quick test against square with side {2*(a+hwd)}: */
    if ((fabs(p->c[0]) > F->a + hwd) || (fabs(p->c[0]) > F->a + hwd)) { return +1.0; }
    /* Convert to UV coordinates and test as aligned: */
    r2_t uvp = (r2_t){{ r2_dot(p, &(F->u)), r2_dot(p, &(F->v)) }};
    return ellipse_aligned_border_position(F->a, F->b, hwd, uvp.c[0], uvp.c[1]);
  }

void ellipse_ouv_print(FILE *wr, ellipse_ouv_t *F, char *fmt)
  { fprintf(wr, "{ u: ( ");
    fprintf(wr, fmt, F->u.c[0]);
    fprintf(wr, " ");
    fprintf(wr, fmt, F->u.c[1]);
    fprintf(wr, " )");
    fprintf(wr, "  a: ");
    fprintf(wr, fmt, F->a);
    fprintf(wr, "   v: ( ");
    fprintf(wr, fmt, F->v.c[0]);
    fprintf(wr, " ");
    fprintf(wr, fmt, F->v.c[1]);
    fprintf(wr, " )  b: ");
    fprintf(wr, fmt, F->b);
    fprintf(wr, "  }");
    fflush(wr);
  }
    
