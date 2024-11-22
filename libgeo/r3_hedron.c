/* See r3_hedron.h */
/* Last edited on 2024-11-22 03:58:26 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <r3.h>
#include <r3_hedron.h>

#include <affirm.h>

#define N 3

/* INTERNAL PROTOTYPES */

void r3_hedron_icosa_vertices_generic(double B, double C, uint32_t n, r3_t r[]);
  /* Stores in {r[0..n-1]} the vertices of an icosahedron of the form
    {(0,±B,±C)} and its cyclic permutations.   Requires {n==12}. */

/* IMPLEMENTATIONS */

void r3_hedron_cylinder(double H, double R, uint32_t m, uint32_t n, double skew, r3_t r[])
  { demand((n >= 0) && (m >= 0), "invalid parameters");
    uint32_t k = 0;
    for (int32_t i = 0; i < m; i++)
      { double hi = (m == 1 ? 0.0 : H*(2*((double)i)/(m-1) - 1));
        double fi = (m == 1 ? 0.0 : skew*(((double)i)/(m-1) - 0.5));
        for (int32_t j = 0; j < n; j++) 
          { double t = 2*M_PI*(j + fi)/n;
            r[k] = (r3_t){{ R*cos(t), R*sin(t), hi }};
            k++;
          }
      }
    assert(k == m*n);
  }

void r3_hedron_tetra_vertices(double R, uint32_t n, r3_t r[])
  { demand(n == 4, "bad num vertices"); 
    double S = R/sqrt(3.0);
    r[0] = (r3_t){{ -S, -S, -S }};
    r[1] = (r3_t){{ +S, +S, -S }};
    r[2] = (r3_t){{ +S, -S, +S }};
    r[3] = (r3_t){{ -S, +S, +S }};
  }
  
void r3_hedron_octa_vertices(double R, uint32_t n, r3_t r[])
  { demand(n == 6, "bad num vertices"); 
    uint32_t k = 0;
    for (int32_t i = 0; i < N; i++) 
      { r[k] = (r3_t){{ 0,0,0 }}; r[k].c[i] = +R; k++;
        r[k] = (r3_t){{ 0,0,0 }}; r[k].c[i] = -R; k++;
      }
    assert(k == n);
  }
          
void r3_hedron_hexa_vertices(double R, uint32_t n, r3_t r[])
  { demand(n == 8, "bad num vertices"); 
    double S = R/sqrt(3.0);
    uint32_t k = 0;
    for (int32_t s0 = -1; s0 <= +1; s0 += 2) 
      { for (int32_t s1 = -1; s1 <= +1; s1 += 2) 
          { for (int32_t s2 = -1; s2 <= +1; s2 += 2) 
              { r[k] = (r3_t){{ s0*S, s1*S, s2*S }}; k++; }
          }
      }
    assert(k == n);
  }
          
void r3_hedron_icosa_vertices(double R, uint32_t n, r3_t r[])
  { demand(n == 12, "bad num vertices"); 
    double s = 1.0/sqrt(5.0);
    double c = sqrt((1 - s)/2);
    double b = sqrt((1 + s)/2);
    assert(fabs(2*c - sqrt((b-c)*(b-c) + b*b + c*c)) < 0.000001);
    assert(fabs(1 - hypot(b, c)) < 0.000001);
    double B = b*R;
    double C = c*R;
    r3_hedron_icosa_vertices_generic(B, C, n, r);
  }
           
void r3_hedron_dodeca_vertices(double R, uint32_t n, r3_t r[])
  { demand(n == 20, "bad num vertices");
    double f = (sqrt(5) + 1)/2;
    double s = sqrt(3);
    double B = R/s*f;
    double C = R/s/f;
    r3_hedron_icosa_vertices_generic(B, C, 12, r);
    r3_hedron_hexa_vertices(R, 8, &(r[12]));
  }

void r3_hedron_icosa_vertices_generic(double B, double C, uint32_t n, r3_t r[])
  {
    int32_t k = 0;
    for (int32_t i = 0; i < N; i++) 
      { for (int32_t sb = -1; sb <= +1; sb += 2) 
          { for (int32_t sc = -1; sc <= +1; sc += 2) 
              { int32_t ja = i;
                int32_t jb = (i + 1) % (int32_t)N;
                int32_t jc = (i + 2) % (int32_t)N;
                r[k].c[ja] = 0.0;
                r[k].c[jb] = sb*B;
                r[k].c[jc] = sc*C;
                k++;
              }
          }
      }
    assert(k == n);
  }
