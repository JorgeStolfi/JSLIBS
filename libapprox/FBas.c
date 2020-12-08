/* See FBas.h
**
** Last edited on 2012-12-15 09:42:05 by stolfilocal
**
** Copyright © 2003 by Jorge Stolfi and Anamaria Gomide, the
** University of Campinas, Brazil. See the rights and conditions
** notice at the end of this file.
*/

#include <FBas.h>
#include <FBasPoly2.h>

#include <bool.h> 
#include <rn.h>
#include <gauss_elim.h>
#include <filefmt.h>
#include <nget.h>
#include <fget.h>
#include <affirm.h> 
#include <jsstring.h> 
    
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define T FBas

/* INTERNAL PROTOTYPES */

void FBas_WriteAsSubclass(FBas *bas, FILE *wr);
  /* Casts {bas} to its widest proper effective subclass, and 
    calls the corresponding {write} method (which must shadow,
    but must not override, the {write} method of {bas}. */

void FBas_PrintSystem(int ne, int nc, double *A, int nb);
  /* Prints the system matrix {A} which is assumed to have 
    {ne} rows and {nc} columns. Prints dividing lines after the 
    first {nb} equations and unknowns. */

/* DATA-FITTING SYSTEM CONSTRUCTION */

/* The procedures in this section set up a linear system {M X + U = 0}
  arising from the problem of fitting a function to a set of given
  data. The matrices {M} and {U} have sizes {ne×ne} and {ne×dv},
  respectively, where {dv} is the number of data values per sample,
  and {ne} is chosen by the procedure.   
  
  The input parameters {bas}, {nb}, {np}, {dp}, {P}, {dv}, and {V}
  have the same meaning as for {FBas_Fit}. 
  
  The matrices {M} and {U} are computed from the given samples and the
  basis functions, and are packaged as a single matrix {A} with size
  {ne×nc} where {nc = ne+dv}. The matrix {A} is allocated by the
  procedure and returned through {*AP}. The sizes {ne} and {nc} are
  returned through {*neP} and {*ncP}, respectively.
  
  The matrix {X}, with size {ne×dv}, contains the unknowns of the
  system, and is neither stored nor computed by these procedures. Once
  the system is solved, first {nb} rows of {X} will be the
  coefficients of the desired approximant in the basis {bas}. */

void FBas_BuildLeastSquaresSystem
  ( FBas *bas,      /* Basis to use. */  
    int nb,         /* Number of basis elements to use. */
    int np,         /* Number of samples. */
    int dp,         /* Coordinates per sample. */
    double *P,      /* Sample positions. */ 
    int dv,         /* Data values per sample. */
    double *V,      /* Sample values. */ 
    double bias,    /* Bias for the internal energy. */
    int *neP,       /* OUT: Number of equations (rows of system matrix). */
    int *ncP,       /* OUT: Number of variables and indep terms (cols). */
    double **AP     /* OUT: Matrix of linear system, stored by rows. */
  );
  /* Sets up a linear system {M X + U = 0} for data fitting by least
    squares (minimizing the sum of squared deviations at sample
    points), plus {bias} times the approximant's energy as defined by
    the given basis {bas}. */

void FBas_BuildInterpolationSystem
  ( FBas *bas,      /* Basis to use. */  
    int nb,         /* Number of basis elements to use. */
    int np,         /* Number of samples. */
    int dp,         /* Coordinates per sample. */
    double *P,      /* Sample positions. */ 
    int dv,         /* Data values per sample. */
    double *V,      /* Sample values. */ 
    int *neP,       /* OUT: Number of equations (rows of system matrix). */
    int *ncP,       /* OUT: Number of variables and indep terms (cols). */
    double **AP     /* OUT: Matrix of linear system, stored by rows. */
  );
  /* Sets up a linear system {M X + U = 0} for data fitting by
    minimizing the the approximant's energy, as defined by the given
    basis {bas}, subject to the constraint that the approximant
    matches the given data values at all sample points. */

/* IMPLEMENTATIONS */

FBas *FBas_Cast(OBJ *bas)
  { if ((bas != NULL) && isprefix(FBas_TypeId, ((FBas *)bas)->type))
      { return (FBas *)bas; }
    else
      { return NULL; }
  }

#define FBas_FileFormat "2002-11-12"

void FBas_M_Write(T *bas, FILE *wr)
  { char *t = bas->type;
    affirm(isprefix(FBas_TypeId, t), "type mismatch");
    filefmt_write_header(wr, "FBas", FBas_FileFormat);
    fprintf(wr, "type = %s\n", t);
    FBas_WriteAsSubclass(bas, wr);
    filefmt_write_footer(wr, "FBas");
    fflush(wr);
  }
  
void FBas_WriteAsSubclass(FBas *bas, FILE *wr)
  { FBasPoly2 *fp = FBasPoly2_Cast(bas);
    if (fp != NULL) { fp->m->write(fp, wr); return; } 
    fprintf(stderr, "unknown FBas subclass \"%s\"\n", bas->type);
    affirm(FALSE, "aborted");
  }

FBas *FBas_Read(FILE *rd)
  { char *t;
    FBas *bas = NULL;
    filefmt_read_header(rd, "FBas", FBas_FileFormat);
    t = nget_string(rd, "type"); fget_eol(rd);
    affirm(isprefix(FBas_TypeId, t), "incorrect FBas type");
    if (isprefix(FBasPoly2_TypeId, t))
      { bas = (FBas *)FBasPoly2_Read(rd); }
    else 
      { fprintf(stderr, "unknown FBas subclass \"%s\"\n", t);
        affirm(FALSE, "aborted");
      }
    fget_skip_formatting_chars(rd);
    filefmt_read_footer(rd, "FBas");
    free(t);
    return bas;
  }

/* DATA FITTING */

bool_t FBas_Debug = TRUE;

void FBas_Fit
  ( int np,         /* Number of samples. */
    int dp,         /* Number of coordinates per sample. */
    double *P,      /* Sample positions. */ 
    int dv,         /* Number of data values per sample. */
    double *V,      /* Sample values. */  
    FBas *bas,      /* Basis to use. */  
    int nb,         /* How many elements of {bas} to use. */
    double *W       /* OUT: coefficients of fit. */
  )
  { if(FBas_Debug)
      { fprintf(stderr, "-- samples --\n");
        int k;
        for (k = 0; k < np; k++)
          { double *Pk = &(P[dp*k]);
            double *Vk = &(V[dv*k]);
            fprintf(stderr, "P[%02d] =", k); 
            rn_gen_print(stderr, dp, Pk, "%9.5f", "( ", " ", " )"); 
            fprintf(stderr, "  V = "); 
            rn_gen_print(stderr, dv, Vk, "%9.5f", "( ", " ", " )"); 
            fprintf(stderr, "\n");
          }
      }
     
    /* Build linear system {A[0..ne-1,0..nc-1]}: */
    int ne, nc;
    double *A;
    if (nb < np)
      { /* Use least-squares with slight bias from energy. */
        double bias = 1.0e-5; /* Must do something smarter... */
        FBas_BuildLeastSquaresSystem
          ( bas, nb, np, dp, P, dv, V, bias, &ne, &nc, &A );
      }
    else
      { /* Use energy minimization with interpolation constraints. */
        FBas_BuildInterpolationSystem
          ( bas, nb, np, dp, P, dv, V, &ne, &nc, &A );
      }
    /* Solve system: */   
    if (FBas_Debug) { FBas_PrintSystem(ne, nc, A, nb); }
    gsel_triangularize(ne, nc, A, TRUE, 0.0);
    gsel_diagonalize(ne, nc, A, TRUE);
    gsel_normalize(ne, nc, A, TRUE);
    if (FBas_Debug) { FBas_PrintSystem(ne, nc, A, nb); }

    /* Extract solution: */
    int nx = nc - dv;
    double X[nx*dv];
    gsel_extract_solution(ne, nc, A, dv, X, TRUE);
    free(A);

    if(FBas_Debug)
      { fprintf(stderr, "-- system solution --\n");
        int i;
        for (i = 0; i < nx; i++)
          { double *Xi = &(X[dv*i]); 
            fprintf(stderr, "X[%02d] = ", i); 
            rn_gen_print(stderr, dv, Xi, "%9.5f", "( ", " ", " )");
            fprintf(stderr, "\n");
          }
      }
    /* Copy coeffs to {W}: */
    int i;
    for(i = 0; i < nb; i++)
      { int r;
        double *Wir = &(W[dv*i]);
        double *Xir = &(X[dv*i]);
        for (r = 0; r < dv; r++, Wir++, Xir++) 
          { (*Wir) = (*Xir); }
      }
  }
  
void FBas_BuildLeastSquaresSystem
  ( FBas *bas,      /* Basis to use. */  
    int nb,         /* Number of basis elements to use. */
    int np,         /* Number of samples. */
    int dp,         /* Coordinates per sample. */
    double *P,      /* Sample positions. */ 
    int dv,         /* Data values per sample. */
    double *V,      /* Sample values. */ 
    double bias,    /* Bias for the internal energy. */
    int *neP,       /* OUT: Number of equations (rows of system matrix). */
    int *ncP,       /* OUT: Number of variables and indep terms (cols). */
    double **AP     /* OUT: Matrix of linear system, stored by rows. */
  )
  { /* Use least squares-plus-energy minimization. */
    /* Compute matrix {B[i,k] = basis[i](P[k])}: */
    double B[nb*np];
    int i;
    for (i = 0; i < nb; i++)
      { int k;
        double *Bik = &(B[np*i]);
        for (k = 0; k < np; k++, Bik++)
          { double *Pk = &(P[dp*k]);
            (*Bik) = bas->m->eval(bas, i, dp, Pk);
          }
      }
            
    /* Now compute the system matrix {A}: */
    int nx = nb;       /* Total number of unknowns (coeffs). */ 
    int ne = nx;       /* Total number of equations. */
    int nc = nx + dv;  /* Total number of cols in the {A} matrix. */
    double *A = (double *)malloc(ne*nc*sizeof(double));
    affirm(A != NULL, "out of mem");
    
    /* The total energy {E} of a set of coefficients {W[0..nb-1]} is
      the sum of two terms: the sum {S(f)} of the squared deviations
      between the approximant {f(P[k])} and the given intensities
      {V[k]}; and the internal energy {Q(f)} of the approximant, which
      is some basis-specific quadratic function of the coefficients
      {W[i]}. */
      
    /* Build the matrix {A}: */
    for (i = 0; i < ne; i++)
      { /* Compute the main entries for equation {i}: */
        int j;
        double *Mij = &(A[nc*i]);  /* Scans row {i}. */
        double *Mji = &(A[i]);     /* Scans column {i} */
        for (j = 0; j <= i; j++, Mij++, Mji += nc)
          { /* Compute the internal energy term. */
            double Q = bas->m->energy(bas, i, j);
            /* Compute least-squares term: */
            double S = 0.0;
            int k;
            double *Bik = &(B[np*i]);
            double *Bjk = &(B[np*j]);
            for (k = 0; k < np; k++, Bik++, Bjk++)
              { S += (*Bik)*(*Bjk); }
            (*Mij) = (*Mji) = bias*Q + S;
          }
        /* Compute the independent terms: */
        int r;
        double *Uir = &(A[nc*i + nx]);
        for (r = 0; r < dv; r++, Uir++) 
          { double T = 0.0;
            int k;
            double *Bik = &(B[np*i]);
            double *Vkr = &(V[r]);
            for (k = 0; k < np; k++, Bik++, Vkr += dv)
              { T += (*Bik)*(*Vkr); }
            (*Uir) = -T;
          }
      }
    (*neP) = ne;
    (*ncP) = nc;
    (*AP) = A;
  }

void FBas_BuildInterpolationSystem
  ( FBas *bas,      /* Basis to use. */  
    int nb,         /* Number of basis elements to use. */
    int np,         /* Number of samples. */
    int dp,         /* Coordinates per sample. */
    double *P,      /* Sample positions. */ 
    int dv,         /* Data values per sample. */
    double *V,      /* Sample values. */ 
    int *neP,       /* OUT: Number of equations (rows of system matrix). */
    int *ncP,       /* OUT: Number of variables and indep terms (cols). */
    double **AP     /* OUT: Matrix of linear system, stored by rows. */
  )
  { int nm = np;          /* Number of Lagrange multipliers. */
    int nx = nb + nm;     /* Total number of unknowns (coeffs + lambdas). */
    int ne = nx;          /* Total number of equations. */
    int nc = nx + dv;   /* Total number of cols in the {A} matrix. */
    double *A = (double *)malloc(ne*nc*sizeof(double));
    affirm(A != NULL, "out of mem");

    /* The first {nb} equations are energy minimization conditions. */
    /* The remaining {np} equations are interpolation constraints. */
    /* The first {nb} variables are the coeffs of the basis. */
    /* The remaining {np} variables are Lagrange multipliers. */

    /* Build the matrix {A}: */
    int i;
    for (i = 0; i < ne; i++)
      { /* Compute equation {i}: */
        int j, r;
        double *Mij = &(A[nc*i]);
        double *Mji = &(A[i]);
        double *Uir = &(A[nc*i + nx]);
        if (i < nb) 
          { /* Equation is a minimum-energy constraint: */
            for (j = 0; j <= i; j++, Mij++, Mji += nc)
              { /* Compute metric. */
                (*Mij) = (*Mji) = bas->m->energy(bas, i, j);
              }
            /* The independent term is zero: */
            for (r = 0; r < dv; r++, Uir++) { (*Uir) = 0.0; }
          }
        else
          { int k = i - nb; 
            /* The equation is an interpolation constraint at {P[k]}: */
            double *Pk = &(P[dp*k]);
            for (j = 0; j < i; j++, Mij++, Mji += nc)
              { double B = (j < nb ? bas->m->eval(bas, j, dp, Pk) : 0.0);
                (*Mij) = (*Mji) = B;
              }
            /* The independent terms are the observed data at {P[k]}, negated: */
            double *Vkr = &(V[dv*k]);
            for (r = 0; r < dv; r++, Uir++, Vkr++) 
              { (*Uir) = -(*Vkr); }
          }
      }
    (*neP) = ne;
    (*ncP) = nc;
    (*AP) = A;
  }

void FBas_PrintSystem(int ne, int nc, double *A, int nb)
  { fprintf(stderr, "-- systems matrix --\n");
    int i, j, ij = 0;
    for (i = 0; i < ne; i++)
      { if (i == nb)
          { for (j = 0; j < nc; j++)
              { if (j > 0) { fprintf(stderr, "-"); }
                if ((j == nb) || (j == ne)) { fprintf(stderr, "+-"); }
                fprintf(stderr, "------------");
              }
            fprintf(stderr, "\n");
          }
        for (j = 0; j < nc; j++)
          { double Aij = A[ij]; ij++; 
            if (j > 0) { fprintf(stderr, " "); }
            if ((j == nb) || (j == ne)) { fprintf(stderr, "| "); }
            if (Aij != 0.0) 
              { fprintf(stderr, "%12.5f", Aij); }
            else
              { fprintf(stderr, "%12s", "0.0    "); }
          }
        fprintf(stderr, "\n");
      }
  }

