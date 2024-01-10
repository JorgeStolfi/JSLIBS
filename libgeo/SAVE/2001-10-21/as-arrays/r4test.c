/* r4test --- test program for r4.h, r4x4.h  */
/* Last edited on 2001-10-21 22:01:05 by stolfi */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ioprotos.h>
#include <r4.h>
#include <r4x4.h>
#include <js.h>

#define N 4

/* Internal prototypes */

int main (int argc, char **argv);
void test_r4(void);
void test_r4x4(void);
double frac (double x);

int main (int argc, char **argv)
{
  int i;
  srand(1993);
  
  for (i=0; i<100; i++) test_r4();
  for (i=0; i<100; i++) test_r4x4();
  fclose(stderr);
  fclose(stdout);
  return (0);
}

double frac (double x)
{ int i = (int)x;
  double f = x - (double)i;
  return (f);
}

void test_r4(void)
{
  r4_t a, b, c, d;
  double r, s, t;
  int i, j, k, l;
  
  /* r4_zero: */
  r4_zero(a);
  for (i=0; i<N; i++)
    { assert(a[i] == 0.0, "r4_zero error"); }
    
  /* r4_axis: */
  for (k=0; k<N; k++)
    { r4_axis(k, a);
      for (i=0; i<N; i++)
        { assert(a[i] == (i==k ? 1.0 : 0.0), "r4_axis error"); }
    }

  /* r4_throw_cube */
  r4_throw_cube(a);
  for (i=0; i<N; i++)
    { assert(a[i] != a[(i+1)%N], "r4_throw error(1)"); 
      assert(frac(a[i]*256.0) != 0.0, "r4_throw error(1)"); 
      assert((a[i] > -1.0) && (a[i] < 1.0), "r4_throw error(2)"); 
    }

  /* r4_add: */
  r4_throw_cube(a);
  r4_throw_cube(b);
  r4_add(a, b, c);
  for (i=0; i<N; i++)
    { assert(c[i] == a[i] + b[i], "r4_add error"); }

  /* r4_sub: */
  r4_throw_cube(a);
  r4_throw_cube(b);
  r4_sub(a, b, c);
  for (i=0; i<N; i++)
    { assert(c[i] == a[i] - b[i], "r4_sub error"); }

  /* r4_scale: */
  s = frandom();
  r4_throw_cube(a);
  r4_scale(s, a, c);
  for (i=0; i<N; i++)
    { assert(c[i] == s * a[i], "r4_scale error"); }

  /* r4_dist: */
  r4_throw_cube(a);
  r4_throw_cube(b);
  r = r4_dist(a, b);
  s = sqrt(r4_dot(a,a)) + sqrt(r4_dot(b,b));
  t = 0.0;
  for (i=0; i<N; i++)
    { t += (a[i] - b[i])*(a[i] - b[i]); }
  t = sqrt(t);
  assert(fabs(t - r) <= 0.000000001 * s, "r4_dist error");
  
  /* r4_dot: */
  r4_throw_cube(a);
  r4_throw_cube(b);
  r = r4_dot(a, b);
  s = sqrt(r4_dot(a,a)*r4_dot(b,b));
  t = 0.0;
  for (i=0; i<N; i++)
    { t += a[i]*b[i]; }
  assert(fabs(t-r) <= 0.000000001 * s, "r4_dot error(1)");
  for (i=0; i<N; i++)
    { r4_axis(i, a);
      for (j=0; j<N; j++)
        { r4_axis(j, b);
          r = r4_dot(a, b);
          t = (i==j ? 1.0 : 0.0);
          assert(r == t, "r4_dot error(2)");
        }
    }
  
  /* r4_cross: */
  r4_throw_cube(a);
  r4_throw_cube(b);
  r4_throw_cube(c);
  r4_cross(a, b, c, d);
  s = r4_dot(a, a);  
  t = r4_dot(d, a);
  assert (fabs(t) <= 0.000000001 * s, "r4_cross error(1)");
  s = r4_dot(b, b);  
  t = r4_dot(d, b);
  assert (fabs(t) <= 0.000000001 * s, "r4_cross error(2)");
  s = r4_dot(c, c);  
  t = r4_dot(d, c);
  assert (fabs(t) <= 0.000000001 * s, "r4_cross error(3)");
  /* Should check length and sign */
  /* r4_cross on basis vectors: */
  for (i=0; i<N; i++)
    { r4_axis(i, a);
      for (j=0; j<N; j++)
        { r4_axis(j, b);
          for (k=0; k<N; k++)
            { int ijk = 0; double val = 0.0;
              r4_axis(k, c);
              if ((i != j) && (i != k) && (j != k))
                { while ((ijk == i) || (ijk == j) || (ijk == k)) ijk++;
                  val = 1.0;
                  if (i > j) { val = -val; }
                  if (i > k) { val = -val; }
                  if (j > k) { val = -val; }
                  if (i > ijk) { val = -val; }
                  if (j > ijk) { val = -val; }
                  if (k > ijk) { val = -val; }
                }
              r4_cross(a, b, c, d);
              for (l=0; l<N; l++)
                { if (d[l] != (l==ijk ? val : 0.0))
                    { fprintf(stderr, "a = "); r4_print(stderr, a); fputc('\n', stderr);
                      fprintf(stderr, "b = "); r4_print(stderr, b); fputc('\n', stderr);
                      fprintf(stderr, "c = "); r4_print(stderr, c); fputc('\n', stderr);
                      fprintf(stderr, "d = "); r4_print(stderr, d); fputc('\n', stderr);
                      assert (0, "r4_cross error(4)");
                    }
                }
	    }
        }
    }
 
  /* r4_orthize: */
  r4_throw_cube(a);
  r4_throw_cube(b);
  r = r4_orthize(a, b, c); /* Now c == a - r * b */
  s = r4_dot(b, b);  
  t = r4_dot(c, b);
  assert (fabs(t) <= 0.000000001 * s, "r4_orthize error(1)");
  s = r4_dot(c, c);  
  t = r4_dot(c, a);
  assert (fabs(t - s) <= 0.000000001 * s, "r4_orthize error(2)");
  /* Should check coplanarity */

  /* r4_print: */
  r4_throw_cube (a);
  fprintf(stderr, "a = ");
  r4_print(stderr, a);
  fputc('\n', stderr);
}

void test_r4x4(void)
{
}  
    
