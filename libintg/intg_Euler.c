/* See intg_Euler.h. */
/* Last edited on 2003-09-01 17:15:44 by stolfi */

#include <intg_Euler.h>
#include <intg_gen.h>

#include <vec.h>
#include <affirm.h>

#include <math.h>
#include <stdlib.h>

/* INTERNAL PROTOTYPES */

Time Intg_Euler_adjust
  ( OBJ *self,
    Time dt, 
    Dist error, 
    Dist tol,
    Time dtMin, 
    Time dtMax
  );

bool_t Intg_Euler_step
  ( OBJ *self,
    Intg_RHS *rhs,
    Time ta,
    State *sa,
    Velocity *va,
    Time tb,
    State *sb,
    Error *eb
  );
  
/* Maximum scale factor allowed in {adjust}: */
#define MaxScale (10.0)

/* Minimum scale factor allowed in {adjust}: */
#define MinScale (1.0/MaxScale)

/* Reduce computed step by this much, just in case: */
#define Safety   (0.95)

Time Intg_Euler_adjust
  ( OBJ *self,
    Time dt, 
    Dist error, 
    Dist tol,
    Time dtMin, 
    Time dtMax
  )
  { double scale;
    if (error > 0.0)
      { scale = Safety * sqrt(tol/error);
        if (scale < MinScale) { scale = MinScale; }
        if (scale > MaxScale) { scale = MaxScale; }
      }
    else
      { scale = MaxScale; }
    dt = scale * dt;
    if (dt < dtMin) { dt = dtMin; }
    if (dt > dtMax) { dt = dtMax; }
    return dt;
  }

bool_t Intg_Euler_step
  ( OBJ *self,
    Intg_RHS *rhs,
    Time ta,
    State *sa,
    Velocity *va,
    Time tb,
    State *sb,
    Error *eb
  )
  { int n = sa->ne;
    affirm(va->ne == n, "inconsistent sizes"); 
    affirm(sb->ne == n, "inconsistent sizes"); 
    affirm(eb->ne == n, "inconsistent sizes"); 
    int i;
    
    /* Handy names for argument vectors: */
    Coord *Sa = sa->e;
    Coord *Sb = sb->e;
    Speed *Va = va->e;
    Dist *Eb = eb->e;

    /* Constants for Euler step formula: */
    Time dt = tb - ta;

    double Cba = dt;
    
    double Cea = dt;

    /* Error estimation (hack --- should check the theory and do it right!): */
    for (i = 0; i < n; i++)
      { Sb[i] = Sa[i] + Cba*Va[i];
        Eb[i] = Cea*Va[i];
      }

    return FALSE;
  }

Intg_Euler_T *Intg_Euler_new(void)
  { Intg_Euler_T *e = (Intg_Euler_T *)malloc(sizeof(Intg_Euler_T));
    affirm(e != NULL, "out of mem");
    e->intg = (Intg_T)
      { Intg_Euler_TypeId,
        "Intg_Euler_T", 
        &Intg_Euler_step, 
        &Intg_Euler_adjust
      };
    return e;
  }
 
