/* See voxm_path.h */
/* Last edited on 2016-04-04 00:11:42 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <r3.h>
#include <r3_extra.h>
#include <r3_path.h>
#include <affirm.h>

#include <voxm_path.h>
#include <voxm_bezier.h>

voxm_path_state_t voxm_path_state_from_r3_path_state(double t, r3_path_state_t *S, double v)
  { voxm_path_state_t T;
    T.t = t;
    T.p = S->p;
    T.v = (r3_t){{ v*S->M.c[0][0], v*S->M.c[0][1], v*S->M.c[0][2] }}; 
    return T;
  }

r3_path_state_t voxm_path_state_to_r3_path_state(voxm_path_state_t *S)
  { 
    r3_path_state_t T;
    
    /* Copy the position: */
    (T.p) = (S->p);
    
    /* Compute the unit direction vector {u} of the path at {p}: */
    r3_t u;
    double vm = r3_dir(&(S->v), &u);
    if (vm == 0) { u = (r3_t){{ 1.0, 0.0, 0.0 }}; }

    /* Transverse direction {v} - any vector orthogonal to {w}: */
    r3_t v;
    (void)r3_pick_ortho(&u, &v);
    (void)r3_dir(&v, &v);

    /* Third direction {w}: */
    r3_t w;
    r3_cross(&u, &v, &w);
    (void)r3_dir(&w, &w);

    /* Assemble the matrix: */
    r3x3_from_rows(&u, &v, &w,  &(T.M));
    
    return T;
  }

void voxm_path_interpolate_some(voxm_path_state_t S[], int N)
  { 
    if (N >= 3)
      { double t0 = S[0].t;
        r3_t *p0 = &(S[0].p);
        r3_t p1, p2;
        double t3 = S[N-1].t;
        r3_t *p3 = &(S[N-1].p);
        voxm_bezier_from_path_states(&(S[0]), &(S[N-1]), &p1, &p2);
        int im;
        for (im = 1; im <= N-2; im++)
          { voxm_path_state_t *Sm = &(S[im]);
            double fm = ((double)im)/((double)N-1); /* Fractional position of {S[im]}. */
            double tm = (1-fm)*t0 + fm*t3; /* Time for {im}. */
            Sm->t = tm;
            r3_t p01, p12, p23, p012, p123; /* Intermediate points. */
            voxm_bezier_split(t0, t3, p0, &p1, &p2, p3, tm, &p01, &p12, &p23, &p012, &p123, &(Sm->p), &(Sm->v));      
          }
      }
  }

double voxm_path_length_estimate(voxm_path_state_t *S, voxm_path_state_t *T, int order)
  { r3_t *p0 = &(S->p);
    r3_t p1, p2;
    r3_t *p3 = &(T->p);
    voxm_bezier_from_path_states(S, T, &p1, &p2);
    return voxm_bezier_length_estimate(p0, &p1, &p2, p3, order);
  }

void voxm_path_state_print (FILE *f, voxm_path_state_t *S)
  { voxm_path_state_gen_print(f, S, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); }

void voxm_path_state_gen_print
 ( FILE *f,
   voxm_path_state_t *S,
   char *tfmt,
   char *pfmt,
   char *vfmt,
   char *olp, char *osep, char *orp, /* Outer delimiters for matrix. */
   char *ilp, char *isep, char *irp  /* Delimiters for matrix rows and position. */
  )
  {
    if (olp == NULL) { olp = "( "; }
    if (osep == NULL) { osep = "\n  "; }
    if (orp == NULL) { orp = "\n)"; }
    if (ilp == NULL) { ilp = "= ( "; }
    if (isep == NULL) { isep = " "; }
    if (irp == NULL) { irp = " )"; }
    if (tfmt == NULL) { pfmt = "%16.8e"; }
    if (pfmt == NULL) { pfmt = "%16.8e"; }
    if (vfmt == NULL) { vfmt = "%16.8e"; }
    fputs(olp, f);
    
    fputs("t", f);
    fputs(ilp, f);
    fprintf(f, tfmt, &(S->t));
    fputs(irp, f);
    
    fputs(osep, f);
    
    fputs("p", f);
    r3_gen_print(f, &(S->p), pfmt, ilp, isep, irp);
    
    fputs(osep, f);
    
    fputs("v", f);
    r3_gen_print(f, &(S->v), vfmt, ilp, isep, irp);
    
    fputs(orp, f);
    fflush(f);
  }

void voxm_path_state_debug(FILE *wr, voxm_path_state_t *S, char *indent, char *title)
  { char *olp; asprintf(&olp, "%s%s state:\n%s  ", indent, title, indent); 
    char *osep; asprintf(&osep, "\n%s  ", indent); 
    voxm_path_state_gen_print
      ( wr, S, 
        "%.4f", "%.2f", "%+.3f", 
        olp, osep, "\n",
        " = ( ", " ", " )"
      );
    free(olp); free(osep);
  }
