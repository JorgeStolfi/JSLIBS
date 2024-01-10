/* r2test --- test program for r2.h, r2x2.h  */
/* Last edited on 2001-10-21 22:00:11 by stolfi */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ioprotos.h>
#include <r2.h>
#include <r2x2.h>
#include <js.h>

#define N 2

/* Internal prototypes */

int main (int argc, char **argv);
void test_r2(void);
void test_r2x2(void);
double frac (double x);

int main (int argc, char **argv)
{
  int i;
  srand(1993);
  
  for (i=0; i<100; i++) test_r2();
  for (i=0; i<100; i++) test_r2x2();
  fclose(stderr);
  fclose(stdout);
  return (0);
}

double frac (double x)
{ int i = (int)x;
  double f = x - (double)i;
  return (f);
}

void test_r2(void)
{
  r2_t a, b, c;
  double r, s, t;
  int i, j, k;
  
  /* r2_zero: */
  r2_zero(a);
  for (i=0; i<N; i++)
    { assert(a[i] == 0.0, "r2_zero error"); }
    
  /* r2_axis: */
  for (k=0; k<N; k++)
    { r2_axis(k, a);
      for (i=0; i<N; i++)
        { assert(a[i] == (i==k ? 1.0 : 0.0), "r2_axis error"); }
    }

  /* r2_throw_cube */
  r2_throw_cube(a);
  for (i=0; i<N; i++)
    { assert(a[i] != a[(i+1)%N], "r2_throw error(1)"); 
      assert(frac(a[i]*256.0) != 0.0, "r2_throw error(1)"); 
      assert((a[i] > -1.0) && (a[i] < 1.0), "r2_throw error(2)"); 
    }

  /* r2_add: */
  r2_throw_cube(a);
  r2_throw_cube(b);
  r2_add(a, b, c);
  for (i=0; i<N; i++)
    { assert(c[i] == a[i] + b[i], "r2_add error"); }

  /* r2_sub: */
  r2_throw_cube(a);
  r2_throw_cube(b);
  r2_sub(a, b, c);
  for (i=0; i<N; i++)
    { assert(c[i] == a[i] - b[i], "r2_sub error"); }

  /* r2_scale: */
  s = frandom();
  r2_throw_cube(a);
  r2_scale(s, a, c);
  for (i=0; i<N; i++)
    { assert(c[i] == s * a[i], "r2_scale error"); }

  /* r2_dist: */
  r2_throw_cube(a);
  r2_throw_cube(b);
  r = r2_dist(a, b);
  s = sqrt(r2_dot(a,a)) + sqrt(r2_dot(b,b));
  t = 0.0;
  for (i=0; i<N; i++)
    { t += (a[i] - b[i])*(a[i] - b[i]); }
  t = sqrt(t);
  assert(fabs(t - r) <= 0.000000001 * s, "r2_dist error");
  
  /* r2_dot: */
  r2_throw_cube(a);
  r2_throw_cube(b);
  r = r2_dot(a, b);
  s = sqrt(r2_dot(a,a)*r2_dot(b,b));
  t = 0.0;
  for (i=0; i<N; i++)
    { t += a[i]*b[i]; }
  assert(fabs(t-r) <= 0.000000001 * s, "r2_dot error(1)");
  for (i=0; i<N; i++)
    { r2_axis(i, a);
      for (j=0; j<N; j++)
        { r2_axis(j, b);
          r = r2_dot(a, b);
          t = (i==j ? 1.0 : 0.0);
          assert(r == t, "r2_dot error(2)");
        }
    }
  
  /* r2_cross: */
  r2_throw_cube(a);
  r2_throw_cube(b);
  r = r2_cross(a, b);
  /* Should check length and sign */
  /* r2_cross on basis vectors: */
  for (i=0; i<N; i++)
    { r2_axis(i, a);
      for (j=0; j<N; j++)
        { double val = 0.0;
          if ((i == 0) && (j == 1)) { val = 1.0; }
          if ((i == 1) && (j == 0)) { val = -1.0; }
          
          r2_axis(j, b);
          r = r2_cross(a, b);
          assert (r == val, "r2_cross error");
        }
    }
 
  /* r2_orthize: */
  r2_throw_cube(a);
  r2_throw_cube(b);
  r = r2_orthize(a, b, c); /* Now c == a - r * b */
  s = r2_dot(b, b);  
  t = r2_dot(c, b);
  assert (fabs(t) <= 0.000000001 * s, "r2_orthize error(1)");
  s = r2_dot(c, c);  
  t = r2_dot(c, a);
  assert (fabs(t - s) <= 0.000000001 * s, "r2_orthize error(2)");
  /* Should check coplanarity */

  /* r2_print: */
  r2_throw_cube (a);
  fprintf(stderr, "a = ");
  r2_print(stderr, a);
  fputc('\n', stderr);
}

void test_r2x2(void)
{
}  
    
