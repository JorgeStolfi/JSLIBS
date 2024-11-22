/* See r3_path.h */
/* Last edited on 2024-11-20 15:50:56 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <r3.h>
#include <r3_motion.h>
#include <affirm.h>
#include <jsprintf.h>

#include <r3_path.h>
#include <r3_bezier.h>

r3_path_state_t r3_path_state_from_r3_motion_state(double t, r3_motion_state_t *S, double v)
  { r3_path_state_t T;
    T.t = t;
    T.p = S->p;
    T.v = (r3_t){{ v*S->M.c[0][0], v*S->M.c[0][1], v*S->M.c[0][2] }}; 
    return T;
  }

r3_motion_state_t r3_path_state_to_r3_motion_state(r3_path_state_t *S)
  { 
    r3_motion_state_t T;
    
    /* Copy the position: */
    (T.p) = (S->p);
    
    /* Compute the unit direction vector {u} of the path at {p}: */
    r3_t u;
    double vm = r3_dir(&(S->v), &u);
    if (vm == 0) { u = (r3_t){{ 1.0, 0.0, 0.0 }}; }

    /* Transverse directions {v,w} orthonormal to {u}: */
    r3_t v, w;
    (void)r3_throw_ortho_pair(&u, &v, &w);

    /* Assemble the matrix: */
    r3x3_from_rows(&u, &v, &w,  &(T.M));
    
    return T;
  }

void r3_path_interpolate_some(r3_path_state_t S[], uint32_t N)
  { 
    if (N >= 3)
      { double t0 = S[0].t;
        r3_t *p0 = &(S[0].p);
        r3_t p1, p2;
        double t3 = S[N-1].t;
        r3_t *p3 = &(S[N-1].p);
        r3_path_bezier_from_states(&(S[0]), &(S[N-1]), &p1, &p2);
        for (int32_t im = 1; im <= N-2; im++)
          { r3_path_state_t *Sm = &(S[im]);
            double fm = ((double)im)/((double)N-1); /* Fractional position of {S[im]}. */
            double tm = (1-fm)*t0 + fm*t3; /* Time for {im}. */
            Sm->t = tm;
            r3_t p01, p12, p23, p012, p123; /* Intermediate points. */
            r3_bezier_split(t0, t3, p0, &p1, &p2, p3, tm, &p01, &p12, &p23, &p012, &p123, &(Sm->p), &(Sm->v));      
          }
      }
  }

void r3_path_bezier_from_states(r3_path_state_t *S, r3_path_state_t *T, r3_t *p1, r3_t *p2)
  { 
    bool_t debug = TRUE;
    if (debug) { fprintf(stderr, "enter %s\n", __FUNCTION__); }

    if (debug) 
      { r3_path_state_debug(stderr, S, "  ", "initial"); 
        r3_path_state_debug(stderr, T, "  ", "final"); 
      }
    
    r3_t *p0 = &(S->p);
    r3_t *p3 = &(T->p);
    double dt = T->t - S->t;  /* Duration of path. */
    r3_mix(1.0, p0, +dt/3.0, &(S->v), p1);
    r3_mix(1.0, p3, -dt/3.0, &(T->v), p2);
    
    if (debug) { fprintf(stderr, "exit %s\n", __FUNCTION__); }
  }

double r3_path_length_estimate(r3_path_state_t *S, r3_path_state_t *T, int32_t order)
  { r3_t *p0 = &(S->p);
    r3_t p1, p2;
    r3_t *p3 = &(T->p);
    r3_path_bezier_from_states(S, T, &p1, &p2);
    return r3_bezier_length_estimate(p0, &p1, &p2, p3, order);
  }

void r3_path_state_print (FILE *f, r3_path_state_t *S)
  { r3_path_state_gen_print(f, S, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); }

void r3_path_state_gen_print
 ( FILE *f,
   r3_path_state_t *S,
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

void r3_path_state_debug(FILE *wr, r3_path_state_t *S, char *indent, char *title)
  { char *olp = jsprintf("%s%s state:\n%s  ", indent, title, indent); 
    char *osep = jsprintf("\n%s  ", indent); 
    r3_path_state_gen_print
      ( wr, S, 
        "%.4f", "%.2f", "%+.3f", 
        olp, osep, "\n",
        " = ( ", " ", " )"
      );
    free(olp); free(osep);
  }
