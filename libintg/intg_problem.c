/* See intg_problem.h */   
/* Last edited on 2009-02-10 20:37:10 by stolfi */


#include <intg_gen.h>
#include <intg_problem.h>

#include <vec.h>
#include <math.h>

/* INTERNAL PROTOTYPES */

void Quadr_Sol(Time t, State *s);
bool_t Quadr_RHS(Time t, State *s, Velocity *v);

void Quart_Sol(Time t, State *s);
bool_t Quart_RHS(Time t, State *s, Velocity *v);

void Quint_Sol(Time t, State *s);
bool_t Quint_RHS(Time t, State *s, Velocity *v);

void Sinus_Sol(Time t, State *s);
bool_t Sinus_RHS(Time t, State *s, Velocity *v);

void Twang_Sol(Time t, State *s);
bool_t Twang_RHS(Time t, State *s, Velocity *v);

void Gauss_Sol(Time t, State *s);
bool_t Gauss_RHS(Time t, State *s, Velocity *v);

void Hyper_Sol(Time t, State *s);
bool_t Hyper_RHS(Time t, State *s, Velocity *v);

void Notch_Sol(Time t, State *s);
bool_t Notch_RHS(Time t, State *s, Velocity *v);

/* PACK THEM ALL: */

intg_prob_vec_t intg_prob_t_Sample(void)
  {
    intg_prob_vec_t p = intg_prob_vec_new(0);
    int np = 0;
    
    auto void mkproblem(intg_sol_t *sol, Intg_RHS *rhs, char *tag, int n);
    void mkproblem(intg_sol_t *sol, Intg_RHS *rhs, char *tag, int n)
      {
        intg_prob_vec_expand(&p, np);
        intg_prob_t *q = &(p.e[np]);
        q->tag = tag;
        q->n = n;
        q->rhs = rhs;
        q->sol = sol;
        q->t0 = -1.0;
        q->t1 = +1.0;
        q->dt = 0.1;
        q->dtMin = 1.0e-6;
        q->dtMax = 1.0;
        q->sol = sol;
        np++;
      }
    
    mkproblem(&Quadr_Sol, &Quadr_RHS, "quadr", 2);
    mkproblem(&Quart_Sol, &Quart_RHS, "quart", 2);
    mkproblem(&Quint_Sol, &Quint_RHS, "quint", 2);
    mkproblem(&Sinus_Sol, &Sinus_RHS, "sinus", 2);
    mkproblem(&Gauss_Sol, &Gauss_RHS, "gauss", 2);
    mkproblem(&Twang_Sol, &Twang_RHS, "twang", 2);
    mkproblem(&Hyper_Sol, &Hyper_RHS, "hyper", 2);
    mkproblem(&Notch_Sol, &Notch_RHS, "notch", 1);
    
    intg_prob_vec_trim(&p, np);
    return p;
  }

vec_typeimpl(intg_prob_vec_t,intg_prob_vec,intg_prob_t);
  /* Manipulation of integration problem vectors. */

/* THE PROBLEMS */

/* Quadratic */

void Quadr_Sol(Time t, State *s)
  { s->e[0] = t*t;
    s->e[1] = 2.0 * t;
  }

bool_t Quadr_RHS(Time t, State *s, Velocity *v)
  { v->e[0] = s->e[1];
    v->e[1] = 2.0;
    return FALSE;
  }

/* Quartic */

void Quart_Sol(Time t, State *s)
  { s->e[0] = t*(t*t*t - 1.0);
    s->e[1] = 4.0 * t*t*t - 1.0;
  }

bool_t Quart_RHS(Time t, State *s, Velocity *v)
  { v->e[0] = s->e[1];
    v->e[1] = 12.0*t*t;
    return FALSE;
  }

/* Quintic */

void Quint_Sol(Time t, State *s)
  { s->e[0] = t*(t*t*t*t - 1.0);
    s->e[1] = 5.0 * t*t*t*t - 1.0;
  }

bool_t Quint_RHS(Time t, State *s, Velocity *v)
  { v->e[0] = s->e[1];
    v->e[1] = 20.0*t*t*t;
    return FALSE;
  }

/* Sinusoid */

#define SinusW (20.0)

void Sinus_Sol(Time t, State *s)
  { double W = SinusW;
    double ct = cos(W*t);
    double st = sin(W*t);
    s->e[0] = ct;
    s->e[1] = -W*st;
  }

bool_t Sinus_RHS(Time t, State *s, Velocity *v)
  { double W = SinusW;
    v->e[0] = s->e[1];
    v->e[1] = - W*W*s->e[0];
    return FALSE;
  }

/* Decaying sinusoid */

#define TwangW (20.0)
#define TwangK (2.0)

void Twang_Sol(Time t, State *s)
  { double W = TwangW;
    double K = TwangK;
    double et = exp(-K*t);
    double ct = et * cos(W*t);
    double st = et * sin(W*t);
    s->e[0] = ct;
    s->e[1] = st;
  }

bool_t Twang_RHS(Time t, State *s, Velocity *v)
  { double W = TwangW;
    double K = TwangK;
    v->e[0] = - K * s->e[0] - W * s->e[1];
    v->e[1] = - K * s->e[1] + W * s->e[0];
    return FALSE;
  }

/* Gaussian bell */

#define GaussK (3.0)

void Gauss_Sol(Time t, State *s)
  { double K = GaussK;
    double kt = K*t;
    double et = exp(-kt*kt);
    s->e[0] = et;
    s->e[1] = -2.0*K*K*t*et;
  }

bool_t Gauss_RHS(Time t, State *s, Velocity *v)
  { double K = GaussK;
    double r = s->e[1]/s->e[0];
    v->e[0] = s->e[1];
    v->e[1] = s->e[0]*(r*r - 2.0*K*K);
    return FALSE;
  }

/* Hyperbola */

#define HyperE (0.05)

void Hyper_Sol(Time t, State *s)
  { double E = HyperE;
    double et = sqrt(t*t + E*E);
    s->e[0] = et;
    s->e[1] = t/et;
  }

bool_t Hyper_RHS(Time t, State *s, Velocity *v)
  { double E = HyperE;
    double r = 1.0/s->e[0];
    v->e[0] = s->e[1];
    v->e[1] = E*E * r*r*r;
    return FALSE;
  }

/* Notches: log(2 + sin(t)) */

#define NotchK (1.41421356237309504880)
#define NotchW (4.0)

void Notch_Sol(Time t, State *s)
  { double K = NotchK;
    double W = NotchW;
    double et = log(K + sin(W*t));
    s->e[0] = et;
  }

bool_t Notch_RHS(Time t, State *s, Velocity *v)
  { double W = NotchW;
    v->e[0] = W*cos(W*t)/exp(s->e[0]);
    return FALSE;
  }
