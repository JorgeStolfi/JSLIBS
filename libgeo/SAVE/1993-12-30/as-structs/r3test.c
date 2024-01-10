/* r3test --- test program for r3.h, r3x3.h  */

#include <math.h>
#include <stdio.h>
#include <ioprotos.h>
#include <r3.h>
#include <r3x3.h>
#include <js.h>

#define N 3

/* Internal prototypes */

int main (...);
void test_r3(void);
void test_r3x3(void);
double frac (double x);

int main (...)
{
  int i;
  srandom(1993);
  
  for (i=0; i<100; i++) test_r3();
  for (i=0; i<100; i++) test_r3x3();
  fclose(stderr);
  fclose(stdout);
  return (0);
}

double frac (double x)
{ int i = (int)x;
  double f = x - (double)i;
  return (f);
}

void test_r3(void)
{
  r3_t a, b, c;
  double r, s, t;
  int i, j, k;
  
  /* r3_zero: */
  r3_zero(&a);
  for (i=0; i<N; i++)
    { assert(a.c[i] == 0.0, "r3_zero error"); }
    
  /* r3_axis: */
  for (k=0; k<N; k++)
    { r3_axis(k, &a);
      for (i=0; i<N; i++)
        { assert(a.c[i] == (i==k ? 1.0 : 0.0), "r3_axis error"); }
    }

  /* r3_throw_cube */
  r3_throw_cube(&a);
  for (i=0; i<N; i++)
    { assert(a.c[i] != a.c[(i+1)%N], "r3_throw error(1)"); 
      assert(frac(a.c[i]*256.0) != 0.0, "r3_throw error(1)"); 
      assert((a.c[i] > -1.0) && (a.c[i] < 1.0), "r3_throw error(2)"); 
    }

  /* r3_add: */
  r3_throw_cube(&a);
  r3_throw_cube(&b);
  r3_add(&a, &b, &c);
  for (i=0; i<N; i++)
    { assert(c.c[i] == a.c[i] + b.c[i], "r3_add error"); }

  /* r3_sub: */
  r3_throw_cube(&a);
  r3_throw_cube(&b);
  r3_sub(&a, &b, &c);
  for (i=0; i<N; i++)
    { assert(c.c[i] == a.c[i] - b.c[i], "r3_sub error"); }

  /* r3_scale: */
  s = frandom();
  r3_throw_cube(&a);
  r3_scale(s, &a, &c);
  for (i=0; i<N; i++)
    { assert(c.c[i] == s * a.c[i], "r3_scale error"); }

  /* r3_dist: */
  r3_throw_cube(&a);
  r3_throw_cube(&b);
  r = r3_dist(&a, &b);
  t = 0.0;
  for (i=0; i<N; i++)
    { t += (a.c[i] - b.c[i])*(a.c[i] - b.c[i]); }
  assert(sqrt(t) == r, "r3_dist error");
  
  /* r3_dot: */
  r3_throw_cube(&a);
  r3_throw_cube(&b);
  r = r3_dot(&a, &b);
  t = 0.0;
  for (i=0; i<N; i++)
    { t += a.c[i]*b.c[i]; }
  assert(t == r, "r3_dot error(1)");
  for (i=0; i<N; i++)
    { r3_axis(i, &a);
      for (j=0; j<N; j++)
        { r3_axis(j, &b);
          r = r3_dot(&a, &b);
          t = (i==j ? 1.0 : 0.0);
          assert(r == t, "r3_dot error(2)");
        }
    }
  
  /* r3_cross: */
  r3_throw_cube(&a);
  r3_throw_cube(&b);
  r3_cross(&a, &b, &c);
  s = r3_dot(&b, &b);  
  t = r3_dot(&c, &b);
  assert (fabs(t) <= 0.000000001 * s, "r3_cross error(1)");
  s = r3_dot(&a, &a);  
  t = r3_dot(&c, &a);
  assert (fabs(t) <= 0.000000001 * s, "r3_cross error(2)");
  /* Should check length and sign */
  /* r3_cross on basis vectors: */
  for (i=0; i<N; i++)
    { r3_axis(i, &a);
      for (j=0; j<N; j++)
        { int ij = i; double val = 0.0;
          if ((i+1)%N == j) { ij = (j+1)%N; val = 1.0; }
          if (i == (j+1)%N) { ij = (i+1)%N; val = -1.0; }
          
          r3_axis(j, &b);
          r3_cross(&a, &b, &c);
          for (k=0; k<N; k++)
            { assert (c.c[k] == (k==ij ? val : 0.0), "r3_cross error(3)"); }
        }
    }
 
  /* r3_orthize: */
  r3_throw_cube(&a);
  r3_throw_cube(&b);
  r = r3_orthize(&a, &b, &c); /* Now c == a - r * b */
  s = r3_dot(&b, &b);  
  t = r3_dot(&c, &b);
  assert (fabs(t) <= 0.000000001 * s, "r3_orthize error(1)");
  s = r3_dot(&c, &c);  
  t = r3_dot(&c, &a);
  assert (fabs(t - s) <= 0.000000001 * s, "r3_orthize error(2)");
  /* Should check coplanarity */

  /* r3_print: */
  r3_throw_cube (&a);
  fprintf(stderr, "a = ");
  r3_print(stderr, &a);
  fputc('\n', stderr);
}

void test_r3x3(void)
{
}  
    
