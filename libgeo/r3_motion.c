/* See r3_motion.h */
/* Last edited on 2021-06-09 19:54:46 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <r3.h>
#include <bool.h>
#include <r3x3.h>
#include <affirm.h>

#include <r3_motion.h>

void r3_motion_state_canonical(r3_motion_state_t *C)
  { (C->p) = (r3_t){{ 0.0, 0.0, 0.0 }};
    r3x3_ident(&(C->M));
  }

void r3_motion_state_compose(r3_motion_state_t *S, r3_motion_state_t *T, r3_motion_state_t *C)
  { /* Beware that {C} may be the same as {S} or {T}, so use {S,T} before setting {C}: */
    r3_t q;
    r3x3_map_row(&(S->p), &(T->M), &q);
    r3_add(&(T->p), &q, &(C->p));
    r3x3_mul(&(S->M), &(T->M), &(C->M));
  }

void r3_motion_sample_uniform
  ( r3_motion_proc_t path,
    double t0,
    double t1,
    int32_t n,
    bool_t mids,
    double t[],
    r3_motion_state_t S[]
  )
  {
    if (n == 0) { return; }
    if ((t == NULL) && (S == NULL)) { return; }

    demand(t1 >= t0, "invalid {t} interval");
    demand(n >= 1, "{n} must be positive");

    /* Compute the time step {dt} and set {t[0]}: */
    double T = t1 - t0;
    double dt; /* Sampling step. */
    double ts0, ts1; /* First and last sampling times. */
    if (mids)
      { dt = T/n;
        ts0 = t0 + dt/2;
        ts1 = t1 - dt/2;
      }
    else
      { if (T == 0)
          { dt = 0; }
        else
          { demand(n >= 2, "{n} must be at least 2");
            dt = T/(n - 1);
          }
        ts0 = t0;
        ts1 = t1;
      }
    if (t != NULL) { t[0] = ts0; }
    if (S != NULL) { S[0] = path(ts0); }
    if (n >= 2) {
      if (t != NULL) { t[n-1] = ts1; }
      if (S != NULL) { S[n-1] = path(ts1); }
    }
    if (n <= 2) { return; }

    /* Get the other samples: */
    int32_t k;
    for (k = 1; k <= n-2; k++)
      { double rk = ((double) k)/((double) n - 1);
        double tk = (1 - rk)*ts0 + rk*ts1;
        if (t != NULL) { t[k] = tk; }
        if (S != NULL) { S[k] = path(tk); }
      }
  }

void r3_motion_circle(double t, double L, double A, r3_motion_state_t *S)
  {  /* Compute the position and direction of arc at {t}: */
     double X, Y; /* Position of {S(t)} on plane {XY}. */
     double dX, dY; /* Direction of path at {t}. */
     if ((t == 0) || (A == 0))
       { X = t*L;
         Y = 0;
         dX = 1;
         dY = 0.0;
       }
     else
       { double R = L/A; /* Radius of curvature. */
         double at = t*A; /* Angular argument at {t}. */
         double cat = cos(at);
         double sat = sin(at);
         X = R*sat;
         Y = R*(1-cat);
         dX = cat;
         dY = sat;
       }
     /* Build the state: */
     (S->p) = (r3_t){{ X, Y, 0.0 }}; /* Position. */
     r3x3_ident(&(S->M));
     S->M.c[0][0] = +dX;  S->M.c[0][1] = +dY; /* Tangent direction. */
     S->M.c[1][0] = -dY;  S->M.c[1][1] = +dX; /* Transverse horizontal direction. */
  }

void r3_motion_helix(double t, double L, double A, double H, r3_motion_state_t *S)
  {
    /* Generate a circle on the {XY} plane: */
    r3_motion_circle(t, L, A, S);
    if (H != 0)
      { /* Raise the point in {Z}: */
        S->p.c[2] += t*H;
        /* Adjust rows 0 and 2 of the matrix: */
        double a = atan2(H,L); /* Tilt angle of spiral. */
        double ca = cos(a);
        double sa = sin(a);
        int32_t j;
        for (j = 0; j < 3; j++)
          { double oldU = S->M.c[0][j];
            double oldW = S->M.c[2][j];
            double newU = + ca*oldU + sa*oldW;
            double newW = - sa*oldU + ca*oldW;
            S->M.c[0][j] = newU;
            S->M.c[2][j] = newW;
          }
      }
  }

void r3_motion_state_print (FILE *f, r3_motion_state_t *S)
  { r3_motion_state_gen_print(f, S, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); }

void r3_motion_state_gen_print
 ( FILE *f,
   r3_motion_state_t *S,
   char *pfmt,
   char *Mfmt,
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
    if (pfmt == NULL) { pfmt = "%16.8e"; }
    if (Mfmt == NULL) { Mfmt = "%16.8e"; }
    fputs(olp, f);
    fputs("p", f);
    r3_gen_print(f, &(S->p), pfmt, ilp, isep, irp);
    int32_t i;
    for (i = 0; i < 3; i++)
      { fputs(osep, f);
        r3_t vi; r3x3_get_row(&(S->M), i, &vi);
        fputc("uvw"[i], f);
        r3_gen_print(f, &vi, Mfmt, ilp, isep, irp);
      }
    fputs(orp, f);
    fflush(f);
  }

void r3_motion_state_debug(FILE *wr, r3_motion_state_t *S, char *indent, char *title)
  { char *olp; asprintf(&olp, "%s%s state:\n%s  ", indent, title, indent);
    char *osep; asprintf(&osep, "\n%s  ", indent);
    r3_motion_state_gen_print
      ( wr, S,
        "%.2f", "%+9.6f",
        olp, osep, "\n",
        " = ( ", " ", " )"
      );
    free(olp); free(osep);
  }
