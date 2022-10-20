/* See bz_basic.h */
/* Last edited on 2022-10-20 06:27:35 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>

#include <affirm.h>
#include <bool.h>

#include <bz_basic.h>

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

double bz_bernstein(bz_degree_t g, bz_index_t i, double z)
  {
    if ((i < 0) || (i > g)) { return NAN; }
    double v; int32_t e, k;
    double u = 1.0 - z;
    if (i+i > g) { u = z; z = 1.0 - z; i = (bz_index_t)(g - i); }
    /* Compute {v = (1-z)^{g-i}}: */
    v = 1.0; e = g-i;
    while (e > 0) { if ((e & 1) == 1) { v *= u; } u = u*u; e >>= 1; }
    /* Multiply {v} by {choose(g,i)*z^i}: */
    for (k = 0; k < i; k++) { v = v*(g-k)/(k+1)*z; }
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
    double c[],      /* Bezier coeffs of a polynomial piece. */
    double wa,       /* Width of first half. */
    double a[],      /* OUT: Bézier coeffs of first half. */
    double wb,       /* Width of second half. */
    double b[]       /* OUT: Bézier coeffs of second half. */
  )
  { /* If client wants nothing, do nothing: */
    if (a == b) { demand(a == NULL, "outputs are the same vector"); }

    /* Compute relative widths: */
    double w = wa + wb;
    register double ra, rb;
    if (w == 0)
      { ra = rb = 0.5; }
    else if (wa < wb)
      { ra = wa/w; rb = 1 - ra; }
    else
      { rb = wb/w; ra = 1 - rb; }
    
    int32_t i, j;
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
    double c[],      /* Bezier coeffs of a polynomial. */
    double w,        /* Actual width of interval. */
    double d[]       /* OUT: Bezier coeffs of derivative. */
  )
  { int32_t i;
    double f = g;
    for (i = 0; i < g; i++) { d[i] = f*(c[i+1] - c[i])/w; }
  }

void bz_integ
  ( bz_degree_t g,   /* Degree of curve. */
    double c[],       /* Bezier coeffs of a polynomial. */
    double w,        /* Actual width of interval. */
    double d[]        /* OUT: Bezier coeffs of integral. */
  )
  { int32_t i;
    double f = g + 1;
    double sum = 0;
    d[0] = 0;
    for (i = 0; i <= g; i++) { sum += c[i]; d[i+1] = w*sum/f; }
  }

void bz_eval
  ( bz_degree_t g,   /* Degree of curve. */
    double c[],      /* Bezier coeffs. */
    double u,        /* Start of interval. */
    double v,        /* End of interval */
    double x,        /* Argument value. */
    bz_degree_t ord, /* Max derivative order desired. */
    double f[]       /* OUT: Image vector (size {ord+1}). */
  )
  { /* If {ord} is negative, there is nothing to compute: */
    if (ord < 0) { return; }
    /* Get relative coord {wa} and complement {wb}: */
    double wa = x - u, wb = v - x;
    double p[g+1]; /* Bezier coeffs of either first or second piece. */
    bool_t second; /* Coeffs {p[0..g]} belong to the second piece. */
    double w;      /* Width of nominal domainof {p}. */
    if (wa < wb)
      { /* Split at {x}, keep and differentiate *second* half: */
        bz_split(g, c, wa, NULL, wb, p);
        w = wb;
        second = TRUE;
      }
    else
      { /* Split at {x}, keep and differentiate the first half: */
        bz_split(g, c, wa, p, wb, NULL);
        w = wa;
        second = FALSE;
      }
    /* Now differentiate {p} up to order {ord} or until {g} is constant: */
    int32_t k = 0;
    while (k <= ord)
      { f[k] = p[second ? 0: g]; k++;
        if (g <= 0) { break; }
        bz_diff(g, p, w, p); g--;
      }
    /* Fill the resmaining derivatives with zeros: */
    while (k <= ord) { f[k] = 0; k++; }
  }
