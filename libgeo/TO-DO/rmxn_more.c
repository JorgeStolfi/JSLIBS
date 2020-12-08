/* See {rmxn_more.h} */
/* Last edited on 2007-10-13 20:02:48 by stolfi */

  
void rmxn_regular_simplex(int n, double V[], double U[])  
  { /* We have {a+b = (sqrt(n+1)-1)/n} for regularity. */
    /* We have {b - n*a = 1} for zero-centering. */
    /* So {(n+1)*a = (sqrt(n+1)-1)/n - 1}. */
    double N = (double)n;
    double a = ((sqrt(N+1)-1)/N - 1)/(N+1);
    double b = 1 + N*a;
    int i, j;
    /* Set the matrix {V}: */
    for (i = 0; i <= n; i++) 
      { int ni = i*n;
        if (i < n)
          { /* Set row {i} to {u_i + (a,a,...,a)}, for {i} in {0..n-1}: */
            for (j = 0; j < n; j++) { V[ni + j] = (i == j ? 1 : 0) + a; }
          }
        else
          { /* Set the last row to {(-b,-b,...,-b)}: */
            for (j = 0; j < n; j++) { V[ni + j] = -b; }
          }
      }
    if (U != NULL)
      { /* Set {U} to the pseudo-inverse of {V} (special case of {rmxn_barycentric_matrix}). */
        double s = -a/(a*N + 1)
        for (i = 0; i < n; i++) 
          { int ni = i*(n+1);
            for (j = 0; j < n; j++)
              { U[ni + j] = (i == j ? 1 : 0) + s; }
            U[ni + n] = 
          }
      }
  }
