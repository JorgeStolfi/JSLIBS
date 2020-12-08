/* intg_problem - a menagerie of ODE integration tests. */
/* Last edited on 2007-01-04 00:17:30 by stolfi */

#ifndef intg_problem_H
#define intg_problem_H

#include <intg_gen.h>

#include <vec.h>

typedef void intg_sol_t(Time t, State *s);
  /* A procedure that returns in {s} the value at time {t}
    of the correct solution of some ODE problem. */

typedef struct intg_prob_t 
  { char *tag;        /* A tag that identifies the problem. */
    unsigned n;       /* Dimension of state vector. */
    Intg_RHS *rhs;    /* Right-hand side of the equation. */
    Time t0;          /* Integraton start time. */
    Time t1;          /* Integraton stopping time. */
    Time dt;          /* Initial time step. */
    Time dtMin;       /* Minimum time step for adaptive integration. */
    Time dtMax;       /* Maximum time step for adaptive integration. */
    intg_sol_t *sol;  /* The correct solution. */
  } intg_prob_t;
    
vec_typedef(intg_prob_vec_t,intg_prob_vec,intg_prob_t);
  /* A vector of iteegration problems ({intg_prob_t}s). */
  
intg_prob_vec_t intg_prob_t_Sample(void);
 /* Returns a sample of ODE integration problems. */

#endif
