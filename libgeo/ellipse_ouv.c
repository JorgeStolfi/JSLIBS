/* See ellipse_ouv.h */
/* Last edited on 2025-01-03 13:21:46 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <assert.h>

#include <r2.h>
#include <box.h>
#include <jsrandom.h>
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
    
bool_t ellipse_ouv_box_inside(ellipse_ouv_t *F, interval_t B[])
  { for (uint32_t kx = 0; kx <= 1; kx++)
      { for (uint32_t ky = 0; ky <= 1; ky++)
         { r2_t p = (r2_t){{ B[0].end[kx], B[1].end[ky] }};
           if (! ellipse_ouv_inside(F, &p)) { return FALSE; }
         }
      }
    return TRUE;
  }

double ellipse_ouv_border_position(ellipse_ouv_t *F, double hwd, r2_t *p)
  { /* Quick test against square with side {2*(a+hwd)}: */
    if ((fabs(p->c[0]) > F->a + hwd) || (fabs(p->c[0]) > F->a + hwd)) { return +1.0; }
    /* Convert to UV coordinates and test as aligned: */
    r2_t uvp = (r2_t){{ r2_dot(p, &(F->u)), r2_dot(p, &(F->v)) }};
    return ellipse_aligned_border_position(F->a, F->b, hwd, uvp.c[0], uvp.c[1]);
  }

double ellipse_ouv_box_coverage(ellipse_ouv_t *F, interval_t B[], uint32_t N)
  { demand(! box_is_empty(2, B), "box must not be empty");
    demand(box_dimension(2, B) == 2, "box must have dimension 2 (non-empty interior)");
    demand(N >= 1, "invalid num probes {N}");
    
    /* Get the ellipse's bounding box: */
    interval_t BF[2];
    ellipse_ouv_bbox(F, BF);
    
    if (box_disjoint(2, B, BF)) { return 0.0; }
    
    if (ellipse_ouv_box_inside(F, B)) { return 1.0; }

    interval_t BI[2]; /* Intersection of bounding boxes. */
    box_meet(2, BF, B, BI);

    double vol = box_measure(2, B);  /* Area of given box. */
    assert(isfinite(vol) && (vol > 0));
    
    double volI = box_measure(2, BI); /* Area of {BI}. */
    assert(isfinite(volI) && (volI > 0));
    assert(volI <= vol);
    
    /* Choose the counts {NXI,NYI} of sample cols and rows: */
    uint32_t NI = (uint32_t)ceil(N*volI/vol); /* Num of probes in intersection. */
    assert(NI >= 1);
    double dx = HI(BI[0]) - LO(BI[0]);
    double dy = HI(BI[1]) - LO(BI[1]);
    uint32_t NXI = (uint32_t)ceil(sqrt(NI*dx/dy) + 1.0e-12);
    uint32_t NYI = (uint32_t)ceil(sqrt(NI*dy/dx) + 1.0e-12); 
    assert(NXI*NYI >= NI);
    NI = NXI*NYI;
    
    /* Do {NXI × NYI} probes in {BI}: */
    uint32_t N_hit = 0;
    for (uint32_t kx = 0; kx < NXI; kx++)
      { for (uint32_t ky = 0; ky < NYI; ky++)
          { /* Divide {BI} into {NXI} by {NYI} cells and jitter a probe in each cell: */ 
            r2_t z = (r2_t){{ (kx + drandom())/NXI, (ky + drandom())/NYI }};
            r2_t p;
            box_point_map(2, z.c, BI, p.c);
            if (ellipse_ouv_inside(F, &p)) { N_hit++; }
          }
      }
    double cov = (N_hit + 0.5)/(NI + 1.0);
    assert(isfinite(cov) && (cov > 0.0) && (cov < 1.0));
    if (volI != vol)
      { /* Account for the area of original box outside intersection: */
        cov = cov*volI/vol;
      }
    return cov;
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
    
