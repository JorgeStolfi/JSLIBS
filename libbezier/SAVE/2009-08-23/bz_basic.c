/* See bz_basic.h */
/* Last edited on 2007-01-04 03:39:52 by stolfi */

#include <bz_basic.h>

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

double bz_bernstein(bz_degree_t g, bz_index_t i, double z)
  {
    double v; int e, k;
    double u = 1.0 - z;
    if (i+i > g) { u = z; z = 1.0 - z; i = g - i; }
    /* Compute {v = (1-z)^{g-i}}: */
    v = 1.0; e = g-i;
    while (e > 0) { if ((e & 1) == 1) { v *= u; } u = u*u; e >>= 1; }
    /* Multiply {v} by {choose(g,i)*z^i}: */
    for (k = 0; k < i; k++)
      { v = v*(g-k)/(k+1)*z; }
    return v;
  }

double bz_bernstein_max(bz_degree_t g, bz_index_t i)
  { if ((i == 0) || (i == g))
      { return 1.0; }
    else
      { double z = ((double)i)/((double)g);
        return bz_bernstein(g, i, z);
      }
  }
 
void bz_split
  ( bz_degree_t g,   /* Degree of curve. */
    double *c,       /* Bezier coeffs of a polynomial piece. */
    double wa,       /* Width of first half. */
    double *a,       /* OUT: Bézier coeffs of first half. */
    double wb,       /* Width of second half. */
    double *b        /* OUT: Bézier coeffs of second half. */
  )
  { /* If client wants nothing, do nothing: */
    if (a == b) { assert(a == NULL); return; }

    /* Compute relative widths: */
    double w = wa+wb;
    register double ra = wa/w, rb = wb/w;
    
    int i, j;
    if (b == NULL)
      { /* Clients wants only the first piece. */
        /* Copy coeffs to {a}, unless they are already there: */
        if (c != a) { for (i = 0; i <= g; i++) { a[i] = c[i]; } }
        for (i = 0; i <= g; i++)
          { for (j = g; j > i; j--) { a[j] = rb*a[j-1] + ra*a[j]; } }
      }
    else if (a == NULL)
      { /* Clients wants only the second piece. */
        /* Copy coeffs to {b}, unless it is already there: */
        if (c != b) { for (i = 0; i <= g; i++) { b[i] = c[i]; } }
        for (i = g; i >= 0; i--)
          { for (j = 0; j < i; j++) { b[j] = rb*b[j] + ra*b[j+1]; } }
      }
    else
      { /* Clients wants both pieces. */
        /* Copy coeffs to {a}, unless they are already there: */
        if (c != a) { for (i = 0; i <= g; i++) { a[i] = c[i]; } }
        for (i = 0; i <= g; i++)
          { b[g-i] = a[g];
            for (j = g; j > i; j--) { a[j] = rb*a[j-1] + ra*a[j]; }
          }
      }
  }

void bz_diff
  ( bz_degree_t g,   /* Degree of curve. */
    double w,        /* Actual width of interval. */
    double *c,       /* Bezier coeffs of a polynomial. */
    double *d        /* OUT: Bezier coeffs of derivative. */
  )
  { int i;
    double f = g;
    for (i = 0; i <= g; i++) { d[i] = f*(c[i+1] - c[i])/w; }
  }

void bz_integ
  ( bz_degree_t g,   /* Degree of curve. */
    double w,        /* Actual width of interval. */
    double *c,       /* Bezier coeffs of a polynomial. */
    double *d        /* OUT: Bezier coeffs of integral. */
  )
  { int i;
    double f = g + 1;
    double sum = 0;
    d[0] = 0;
    for (i = 0; i <= g; i++) { sum += c[i]; d[i+1] = w*sum/f; }
  }

void bz_eval
  ( bz_degree_t g,   /* Degree of curve. */
    double *c,       /* Bezier coeffs. */
    double u,        /* Start of interval. */
    double v,        /* End of interval */
    double x,        /* Argument value. */
    bz_degree_t ord, /* Max derivative order desired. */
    double f[]       /* OUT: Image vector (size {n}). */
  )
  { /* If {ord} is negative, there is nothing to compute: */
    if (ord < 0) { return; }
    /* Get relative coord {wa} and complement {wb}: */
    double wa = x - u, wb = v - x;
    double p[g+1];
    int k;
    if (wa < wb)
      { /* Split at {wa:wb}, differentiate *second* half: */
        bz_split(g, c, wa, NULL, wb, p);
        /* Now differentiate {p} up to order {ord}, and keep the *first* coeffs: */
        for (k = 0; ; k++)
          { f[k] = p[0];
            if (k == ord) { return; }
            bz_diff(g, wb, p, p);
            g--;
          }
      }
    else
      { /* Split the polynomial at {wa:wb}, keep the first half: */
        bz_split(g, c, wa, p, wb, NULL);
        /* Now differentiate {p} up to order {ord}, and keep the *last* coeffs: */
        for (k = 0; ; k++)
          { f[k] = p[g];
            if (k == ord) { return; }
            bz_diff(g, wa, p, p);
            g--;
          }
      }
  }
