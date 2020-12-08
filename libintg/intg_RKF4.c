/* See intg_RKF4.h. */
/* Last edited on 2003-09-01 11:28:02 by stolfi */

#include <intg_RKF4.h>
#include <intg_gen.h>

#include <vec.h>
#include <affirm.h>

#include <math.h>
#include <stdlib.h>

/* INTERNAL PROTOTYPES */

Time Intg_RKF4_adjust
  ( OBJ *self,
    Time dt, 
    Dist error, 
    Dist tol,
    Time dtMin, 
    Time dtMax
  );

bool_t Intg_RKF4_step
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
#define Safety   (0.840896415253)

Time Intg_RKF4_adjust
  ( OBJ *self,
    Time dt, 
    Dist error, 
    Dist tol,
    Time dtMin, 
    Time dtMax
  )
  { double scale;
    if (error > 0.0)
      { scale = Safety * sqrt(sqrt(tol/error));
        if (scale <= MinScale) { scale = MinScale; } 
        if (scale >= MaxScale) { scale = MaxScale; } 
      }
    else
      { scale = MaxScale; }
    dt = scale * dt;
    if (dt < dtMin) { dt = dtMin; }
    if (dt > dtMax) { dt = dtMax; }
    return dt;
  }

bool_t Intg_RKF4_step
  ( OBJ *self,
    Intg_RHS *rhs,
    Time ta,
    State *sa,
    Velocity *va,
    Time tb,
    State *sb,
    Error *eb
  )
  { Intg_RKF4_T *t = (Intg_RKF4_T *)self;
    int n = sa->ne;
    affirm(va->ne == n, "inconsistent sizes"); 
    affirm(sb->ne == n, "inconsistent sizes"); 
    affirm(eb->ne == n, "inconsistent sizes"); 
    int i;
    
    /* Handy names for argument vectors: */
    Coord *Sa = sa->e;
    Speed *Va = va->e;
    Coord *Sb = sb->e;
    Dist *Eb = eb->e;

    /* Adjust size of work areas as needed: */
    double_vec_expand(&(t->s), n-1);
    for (i = 0; i < 5; i++) { double_vec_expand(&(t->v[i]), n-1); }
  
    /* Handy names for work areas: */
    State *sm = &(t->s);       Coord *Sm = sm->e;
    Velocity *v1 = &(t->v[0]); Speed *V1 = v1->e;
    Velocity *v2 = &(t->v[1]); Speed *V2 = v2->e;
    Velocity *v3 = &(t->v[2]); Speed *V3 = v3->e;
    Velocity *v4 = &(t->v[3]); Speed *V4 = v4->e;
    Velocity *v5 = &(t->v[4]); Speed *V5 = v5->e;
    
    /* Constants for Runge-Kutta-Fehlberg step formulas: */

    Time dt = tb - ta;

    double C1a = dt*(1.0/4.0);
    double C1t = (C1a);
    
    double C2a = dt*(3.0/32.0);
    double C21 = dt*(9.0/32.0);
    double C2t = (C2a + C21);
    
    double C3a = dt*(1932.0/2197.0);
    double C31 = dt*(-7200.0/2197.0);
    double C32 = dt*(7296.0/2197.0);
    double C3t = (C3a + C31 + C32);
    
    double C4a = dt*(439.0/216.0);
    double C41 = dt*(-8.0);
    double C42 = dt*(3680.0/513.0);
    double C43 = dt*(-845.0/4104.0);
    double C4t = (C4a + C41 + C42 + C43);
    
    double C5a = dt*(-8.0/27.0);
    double C51 = dt*(2.0);
    double C52 = dt*(-3544.0/2565.0);
    double C53 = dt*(1859.0/4104.0);
    double C54 = dt*(-11.0/40.0);
    double C5t = (C5a + C51 + C52 + C53 + C54);
    
    double Cba = dt*(25.0/216.0);
    double Cb2 = dt*(1408.0/2565.0);
    double Cb3 = dt*(2197.0/4104.0);
    double Cb4 = dt*(-1.0/5.0);

    double Cea = (1.0/360.0);
    double Ce2 = (-128.0/4275.0);
    double Ce3 = (-2197.0/75240.0);
    double Ce4 = (1.0/50.0);
    double Ce5 = (2.0/55.0);
    
    /* We could do without the internal work vector {sm} by using {sb}
      instead. However, if the client gave us the same vector for {sa}
      and {sb}, the result would be garbage.

      We could also use {v5} instead of {sm} in all steps excet the
      degree 5 extrapolation, where we could use {v1}. But that would
      be too obscure...

      It is safer this way... */
      
    /* Degree 1 extrapolation: */
    for (i = 0; i < n; i++)
      { Sm[i] = Sa[i] + C1a*Va[i]; }
    if (rhs(ta + C1t, sm, v1)) { return TRUE; }

    /* Degree 2 extrapolation: */
    for (i = 0; i < n; i++)
      { Sm[i] = Sa[i] + (C2a*Va[i] + C21*V1[i]); }
    if (rhs(ta + C2t, sm, v2)) { return TRUE; }

    /* Degree 3 extrapolation: */
    for (i = 0; i < n; i++)
      { Sm[i] = Sa[i] + (C3a*Va[i] + C31*V1[i] + C32*V2[i]); }
    if (rhs(ta + C3t, sm, v3)) { return TRUE; }

    /* Degree 4 extrapolation: */
    for (i = 0; i < n; i++)
      { Sm[i] = Sa[i] + (C4a*Va[i] + C41*V1[i] + C42*V2[i] + C43*V3[i]); }
    if (rhs(ta + C4t, sm, v4)) { return TRUE; }

    /* Degree 5 extrapolation: */
    for (i = 0; i < n; i++)
      { Sm[i] = Sa[i] + (C5a*Va[i] + C51*V1[i] + C52*V2[i] + C53*V3[i] + C54*V4[i]); }
    if (rhs(ta + C5t, sm, v5)) { return TRUE; }

    /* Final extrapolation and error estimate: */
    for (i = 0; i < n; i++)
      { Sb[i] = Sa[i] + (Cba*Va[i] + Cb2*V2[i] + Cb3*V3[i] + Cb4*V4[i]);
        Eb[i] = Cea*Va[i] + Ce2*V2[i] + Ce3*V3[i] + Ce4*V4[i] + Ce5*V5[i];
      }

    return FALSE;
  }

Intg_RKF4_T *Intg_RKF4_new(void)
  { Intg_RKF4_T *t = (Intg_RKF4_T *)malloc(sizeof(Intg_RKF4_T));
    affirm(t != NULL, "out of mem");
    t->intg = (Intg_T)
      { Intg_RKF4_TypeId,
        "Intg_RKF4_T", 
        &Intg_RKF4_step, 
        &Intg_RKF4_adjust
      };
    t->s = double_vec_new(0);
    int i;
    for(i = 0; i < 5; i++) { t->v[i] = double_vec_new(0); }
    return t;
  }
