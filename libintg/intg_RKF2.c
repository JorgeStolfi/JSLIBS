/* See intg_RKF2.h. */
/* Last edited on 2024-12-05 10:31:28 by stolfi */

#include <intg_RKF2.h>
#include <intg_gen.h>

#include <vec.h>
#include <affirm.h>

#include <math.h>
#include <stdlib.h>

/* INTERNAL PROTOTYPES */

Time Intg_RKF2_adjust
  ( OBJ *self,
    Time dt, 
    Dist error, 
    Dist tol,
    Time dtMin, 
    Time dtMax
  );

bool_t Intg_RKF2_step
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

Time Intg_RKF2_adjust
  ( OBJ *self,
    Time dt, 
    Dist error, 
    Dist tol,
    Time dtMin, 
    Time dtMax
  )
  { double scale;
    if (error > 0.0)
      { scale = Safety * cbrt(tol/error);
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

bool_t Intg_RKF2_step
  ( OBJ *self,
    Intg_RHS *rhs,
    Time ta,
    State *sa,
    Velocity *va,
    Time tb,
    State *sb,
    Error *eb
  )
  { Intg_RKF2_T *t = (Intg_RKF2_T *)self;
    int n = sa->ne;
    affirm(va->ne == n, "inconsistent sizes"); 
    affirm(sb->ne == n, "inconsistent sizes"); 
    affirm(eb->ne == n, "inconsistent sizes"); 
    int i;
    
    /* Handy names for argument vectors: */
    Coord *Sa = sa->e;
    Coord *Sb = sb->e;
    Speed *Va = va->e;
    Dist *Eb = eb->e;

    /* Adjust size of work areas as needed: */
    double_vec_expand(&(t->s), n-1);
    double_vec_expand(&(t->v), n-1);
    
    /* Handy names for work areas: */
    State *sm = &(t->s);      Coord *Sm = sm->e;
    State *v1 = &(t->v);    Speed *V1 = v1->e;
   
    /* Constants for RKF2 step formula: */
    Time dt = tb - ta;

    double C1a = 0.5*dt;
    double C1t = C1a;
    
    double Cba = dt;
    
    double Cea = 0.5*dt;
    double Ce1 = -0.5*dt;

    /* Degree 1 extrapolation: */
    for (i = 0; i < n; i++)
      { Sm[i] = Sa[i] + C1a*Va[i]; }
    if (rhs(ta + C1t, sm, v1)) { return TRUE; }

    /* Degree 2 extrapolation and error estimate: */
    for (i = 0; i < n; i++)
      { Sb[i] = Sa[i] + Cba*Va[i];
        Eb[i] = Cea*Va[i] + Ce1*V1[i];
      }

    return FALSE;
  }

Intg_RKF2_T *Intg_RKF2_new(void)
  { Intg_RKF2_T *e = (Intg_RKF2_T *)malloc(sizeof(Intg_RKF2_T));
    affirm(e != NULL, "out of mem");
    e->intg = (Intg_T)
      { Intg_RKF2_TypeId,
        "Intg_RKF2_T", 
        &Intg_RKF2_step, 
        &Intg_RKF2_adjust
      };
    e->s = double_vec_new(0);
    e->v = double_vec_new(0);
    return e;
  }
 
