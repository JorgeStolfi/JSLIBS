/* See FBasPoly2.h.
**
** Last edited on 2007-01-04 03:42:51 by stolfi
**
** Copyright © 2003 by Jorge Stolfi and Anamaria Gomide, the
** University of Campinas, Brazil. See the rights and conditions
** notice at the end of this file.
*/

#include <FBasPoly2.h>
#include <FBas.h>

#include <filefmt.h>
#include <nget.h>
#include <fget.h>
#include <affirm.h>
#include <jsstring.h>
    
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define T FBasPoly2

double FBasPoly2_M_Eval(T *bas, int index, int dp, double *p);
double FBasPoly2_M_Energy(T *bas, int i, int j);
int FBasPoly2_M_RoundSize(T *bas, int nb, bool_t up);
T *FBasPoly2_M_Copy(T *bas);
void FBasPoly2_M_Free(T *bas);
void FBasPoly2_M_Write(T *bas, FILE *wr);

FBasPoly2_Mths *FBasPoly2_Mths_New(void);
FBasPoly2_Data *FBasPoly2_Data_New(void);
FBasPoly2 *FBasPoly2_T_New(void);
  /* Allocates new records of the specified types.
    All fields are left undefined. */

FBasPoly2 *FBasPoly2_New(void);
  /* Allocates a new {FBasPoly2} object {bas}, plus all its owned
    internal storage, with appropriate {type} field and methods
    record. The contents of the owned storage is left undefined. */

FBasPoly2 *FBasPoly2_Cast(OBJ *o)
  { if ((o != NULL) && isprefix(FBasPoly2_TypeId, ((FBas *)o)->type))
      { return (FBasPoly2 *)o; }
    else
      { return NULL; }
  }

#define FBasPoly2_FileFormat "2002-11-18"

/* OVERRIDES FOR PARENT CLASS METHODS */

double FBasPoly2_M_Eval(T *bas, int index, int dp, double *p)
  { affirm(dp == 2, "wrong point dimension");
    /* Find degree {d} of element {index}: */
    int d = 0, imin = 0;
    while (index > imin + d) { d++; imin += d; }
    /* Compute expoents of element: */
    int ix = index - imin, iy = d - ix;
    /* Evaluate the monomial {z = x^ix*y^iy} */
    /* (don't care aboout speed for now): */
    double x = p[0], y = p[1], z = 1.0;
    int r;
    for (r = 1; r <= ix; r++) { z *= x; }
    for (r = 1; r <= iy; r++) { z *= y; }
    return z;
  }
  
#define f2baspoly_TB_DEG 9
#define f2baspoly_TB_SZ ((f2baspoly_TB_DEG+1)*(f2baspoly_TB_DEG+2)/2)
static int FBasPoly2_EnergyTable[f2baspoly_TB_SZ] = 
  { 1,
    1,1,
    3,1,3,
    10,2,2,10,
    35,5,3,5,35,
    126,14,6,6,14,126,
    462,42,14,10,14,42,462,
    1716,132,36,20,20,36,132,1716,
    6435,429,99,45,35,45,99,429,6435,
    24310,1430,286,110,70,70,110,286,1430,24310
  };
   
double FBasPoly2_M_Energy(T *bas, int i, int j)
  { /* Find degree {d} of element {k := max(i,j)}: */
    int k = (i > j ? i : j);
    int d = 0, imin = 0;
    while (k > imin + d) { d++; imin += d; }
    affirm(k < (d+1)*(d+2)/2, "degree bug");
    affirm(imin == d*(d+1)/2, "imin bug");

    /* The quadratic form {Q} is proportional to the integral
      over a disk of radius {1} of the squared deviation between
      the polynomial {f} and the polynomial of total degree {d-1}
      that osculates {f} to order {d-1} at the center of the disk. 
      That is the integral of the square of the degree-{d} terms
      only. */
      
    if ((i < imin) || (j < imin)) { return 0.0; }
    /* Find exponents {ix} and {jx} of {x} in the monomials {i,j}: */
    int ix = (i - imin);
    int jx = (j - imin);
    if ((ix+jx) % 2 != 0) { return 0.0; }
    /* Look up energy coeff in table: */
    affirm(d <= f2baspoly_TB_DEG, "unimplemented");
    int r = (ix + jx)/2;
    return (double)FBasPoly2_EnergyTable[imin + r];
  }

int FBasPoly2_M_RoundSize(T *bas, int n, bool_t up)
  { /* Use fact that {deg_down(nb) = deg_up(nb+1)-1}: */
    if (! up) { n++; }
    /* Find degree {d} that gives at {nb >= nb} elements: */
    int d = -1, nb = 0; 
    while (nb < n) { d++; nb += d+1; }
    if (! up) { nb -= d+1; d--; affirm(nb >= 0, "bad size"); }
    return nb;
  }
  
T *FBasPoly2_M_Copy(T *bas)
  { affirm(isprefix(FBasPoly2_TypeId, bas->type), "type/method bug");
    FBasPoly2 *cpy = FBasPoly2_New();
    return cpy;
  }

void FBasPoly2_M_Free(T *bas)
  { affirm(isprefix(FBasPoly2_TypeId, bas->type), "type/method bug");
    free(bas->d);
    free(bas);
  }

/* CLASS-SPECIFIC METHODS */
  
void FBasPoly2_M_Write(T *bas, FILE *wr)
  { filefmt_write_header(wr, "FBasPoly2", FBasPoly2_FileFormat);
    filefmt_write_footer(wr, "FBasPoly2");
    fflush(wr);
  }
    
/* OTHER PROCS */
  
FBasPoly2_Mths *FBasPoly2_Mths_New(void)
  { void *v = malloc(sizeof(FBasPoly2_Mths));
    return (FBasPoly2_Mths *)notnull(v, "no mem for FBasPoly2_Mths");
  }

FBasPoly2_Data *FBasPoly2_Data_New(void)
  { affirm(sizeof(FBasPoly2_Data) == 0, "empty record size");
    return NULL;
    /* void *v = malloc(sizeof(FBasPoly2_Data)); */
    /* return (FBasPoly2_Data *)notnull(v, "no mem for FBasPoly2_Data"); */
  }

FBasPoly2 *FBasPoly2_T_New(void)
  { void *v = malloc(sizeof(FBasPoly2));
    return (FBasPoly2 *)notnull(v, "no mem for FBasPoly2");
  }

static FBasPoly2_Mths *Mths = NULL;

FBasPoly2 *FBasPoly2_New(void)
  { FBasPoly2 *bas = FBasPoly2_T_New();
    bas->type = FBasPoly2_TypeId;
    bas->d = FBasPoly2_Data_New();
    if (Mths == NULL)
      { Mths = FBasPoly2_Mths_New();
        Mths->pm.eval = (FBas_EvalMth *)&FBasPoly2_M_Eval;
        Mths->pm.energy = (FBas_EnergyMth *)&FBasPoly2_M_Energy;
        Mths->pm.roundsz = (FBas_RoundSizeMth *)&FBasPoly2_M_RoundSize;
        Mths->pm.copy = (FBas_CopyMth *)&FBasPoly2_M_Copy;
        Mths->pm.free = (FBas_FreeMth *)&FBasPoly2_M_Free;
        /* Note: the {fn.write} method is inherited from {FBas}! */
        Mths->pm.write = (FBas_WriteMth *)&FBas_M_Write;
        /* Class-specific methods */
        Mths->write = (FBas_WriteMth *)&FBasPoly2_M_Write;
      }
    bas->m = Mths;
    return bas;
  }
    
FBasPoly2 *FBasPoly2_Read(FILE *rd)
  { FBasPoly2 *bas = FBasPoly2_New();
    filefmt_read_header(rd, "FBasPoly2", FBasPoly2_FileFormat);
    filefmt_read_footer(rd, "FBasPoly2");
    return bas;
  }
