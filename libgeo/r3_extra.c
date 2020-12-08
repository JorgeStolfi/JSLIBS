/* See r3_extra.h */
/* Last edited on 2014-01-12 15:38:49 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <r3.h>
#include <r3_extra.h>

#include <affirm.h>

#define N 3

/* INTERNAL PROTOTYPES */

void r3_icosahedron_gen_vertices(double B, double C, int n, r3_t r[]);
  /* Stores in {r[0..n-1]} the vertices of an icosahedron of the form
    {(0,±B,±C)} and its cyclic permutations.   Requires {n==12}. */

/* IMPLEMENTATIONS */

double r3_pick_ortho (r3_t *u, r3_t *r)
  {
    if (u !=  r) { (*r) = (*u); }
    double m = r3_L_inf_norm(u);
    if (m != 0)
      { int i;
        if (m == fabs(u->c[0]))
          { i = 0; }
        else if (m == fabs(u->c[1]))
          { i = 1; }
        else if (m == fabs(u->c[2]))
          { i = 2; }
        else
          { assert(FALSE); }
        int j = (i + 1) % N;
        int k = (i + 2) % N;
        double t = r->c[i];  
        r->c[i] = -(r->c[j]);
        r->c[j] = t; 
        r->c[k] = 0;
      }
    return m;
  }

void r3_cylindrical_grid(double H, double R, int m, int n, double skew, r3_t r[])
  { demand((n >= 0) && (m >= 0), "invalid parameters");
    int k = 0;
    int i, j;
    for (i = 0; i < m; i++)
      { double hi = (m == 1 ? 0.0 : H*(2*((double)i)/(m-1) - 1));
        double fi = (m == 1 ? 0.0 : skew*(((double)i)/(m-1) - 0.5));
        for (j = 0; j < n; j++) 
          { double t = 2*M_PI*(j + fi)/n;
            r[k] = (r3_t){{ R*cos(t), R*sin(t), hi }};
            k++;
          }
      }
    assert(k == m*n);
  }

void r3_tetrahedron_vertices(double R, int n, r3_t r[])
  { demand(n == 4, "bad num vertices"); 
    double S = R/sqrt(3.0);
    r[0] = (r3_t){{ -S, -S, -S }};
    r[1] = (r3_t){{ +S, +S, -S }};
    r[2] = (r3_t){{ +S, -S, +S }};
    r[3] = (r3_t){{ -S, +S, +S }};
  }
  
void r3_octahedron_vertices(double R, int n, r3_t r[])
  { demand(n == 6, "bad num vertices"); 
    int k = 0;
    int i;
    for (i = 0; i < N; i++) 
      { r[k] = (r3_t){{ 0,0,0 }}; r[k].c[i] = +R; k++;
        r[k] = (r3_t){{ 0,0,0 }}; r[k].c[i] = -R; k++;
      }
    assert(k == n);
  }
          
void r3_hexahedron_vertices(double R, int n, r3_t r[])
  { demand(n == 8, "bad num vertices"); 
    double S = R/sqrt(3.0);
    int k = 0;
    int s0, s1, s2;
    for (s0 = -1; s0 <= +1; s0 += 2) 
      { for (s1 = -1; s1 <= +1; s1 += 2) 
          { for (s2 = -1; s2 <= +1; s2 += 2) 
              { r[k] = (r3_t){{ s0*S, s1*S, s2*S }}; k++; }
          }
      }
    assert(k == n);
  }
          
void r3_icosahedron_vertices(double R, int n, r3_t r[])
  { demand(n == 12, "bad num vertices"); 
    double s = 1.0/sqrt(5.0);
    double c = sqrt((1 - s)/2);
    double b = sqrt((1 + s)/2);
    assert(fabs(2*c - sqrt((b-c)*(b-c) + b*b + c*c)) < 0.000001);
    assert(fabs(1 - hypot(b, c)) < 0.000001);
    double B = b*R;
    double C = c*R;
    r3_icosahedron_gen_vertices(B, C, n, r);
  }
           
void r3_dodecahedron_vertices(double R, int n, r3_t r[])
  { demand(n == 20, "bad num vertices");
    double f = (sqrt(5) + 1)/2;
    double s = sqrt(3);
    double B = R/s*f;
    double C = R/s/f;
    r3_icosahedron_gen_vertices(B, C, 12, r);
    r3_hexahedron_vertices(R, 8, &(r[12]));
  }

void r3_icosahedron_gen_vertices(double B, double C, int n, r3_t r[])
  {
    int k = 0;
    int i, sb, sc;
    for (i = 0; i < N; i++) 
      { for (sb = -1; sb <= +1; sb += 2) 
          { for (sc = -1; sc <= +1; sc += 2) 
              { int ja = i;
                int jb = (i + 1) % N;
                int jc = (i + 2) % N;
                r[k].c[ja] = 0.0;
                r[k].c[jb] = sb*B;
                r[k].c[jc] = sc*C;
                k++;
              }
          }
      }
    assert(k == n);
  }
