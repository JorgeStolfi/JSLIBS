/* See {intgprod.h} */
/* _*_ coding: iso-8859-1 _*_ */
/* Last edited on 2009-03-05 10:37:25 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <interval.h>
#include <affirm.h>
    
#include <intgprod.h>

intgprod_node_t *intgprod_node_new(void)
  { intgprod_node_t *t = notnull(malloc(sizeof(intgprod_node_t)), "no mem");
    *t = (intgprod_node_t){
      .xm = NAN,
      .ch = {NULL,NULL},
      .IoF = (interval_t){{ 0, +INF }},
      .fP = (interval_t){{ 0, +INF }},
      .nv = 0,
      .iv = NULL,
      .fY = NULL
    };
    return t;
  }

void intgprod_node_free(intgprod_node_t *t)
  {
    if (t->fY != NULL) { free(t->fY); }
    if (t->iv != NULL) { free(t->iv); }
    free(t);
  }
      
void intgprod_condense_factor_list(intgprod_node_t *t, double minVar)
  {
    int old_nv = t->nv; /* Number of existing explicit factors. */
    int new_nv = 0;     /* Number of surviving explicit factors. */
    int j;
    for (j = 0; j < old_nv; j++)
      { interval_t *fYi = &(t->fY[j]);
        assert(fYi->end[0] >= 0);
        if (interval_rad(fYi) >= minVar * fYi->end[0] / 2)
          { t->iv[new_nv] = t->iv[j];
            t->fY[new_nv] = *fYi;
            new_nv++;
          }
        else 
          { t->fP.end[0] *= fYi->end[0];
            t->fP.end[1] *= fYi->end[1];
          }
      }
    if (new_nv != old_nv)
      { /* Assuming that {t->iv} and {t->fY} were allocated with {malloc}: */
        t->iv = realloc(t->iv, new_nv*sizeof(int));
        t->fY = realloc(t->iv, new_nv*sizeof(interval_t));
        t->nv = new_nv;      
      }
  }

interval_t intgprod_compute_IoF(interval_t *X, int n, interval_t Y[])
  { assert(X->end[1] >= X->end[0]);       /* The {X} interval must be non-empty. */
    interval_t F = (interval_t){{ 1,1 }}; /* Enclosure for {F(X[])}. */
    interval_t *Y = t->[]->Y;
    for (i = 0; i < n; i++)
      { /* Multiply {Y[i]} into {F}: */
        interval_t *Yi = &(Y[i]);
        assert(Yi->end[0] >= 0);          /* Must be non-negative. */
        assert(Yi->end[1] >= Yi->end[0]); /* Must be non-empty. */
        F.end[0] *= Yi->end[0];
        F.end[1] *= Yi->end[1];
      }
    double dX = X->end[1] - X->end[0];
    return (interval_t){{ F.end[0]*dX, F.end[1]*dX }};
  }

void intgprod_node_split
  ( intgprod_node_t *t, 
    interval_t *X, 
    intgprod_est_func_t *est_f
  )
  {
    auto void intgprod_node_leaf_split(t, X, n, est_f);

    if ((t->ch[0] == NULL) || (t->ch[1] == NULL))
      { /* Turn {t} into a split node with two leaf children: */
        void intgprod_split_node_leaf(t, X, n, est_f);
      }
    else
      { /* Compute the contributions of children to the uncertainty of {t->IoF}: */
        double v0 = interval_rad(t->ch[0]->IoF);
        double v1 = interval_rad(t->ch[1]->IoF);
        /* Split some leaf of the subtree with largest contribution: */
        interval_t chX;
        intgprod_node_t *ch;
        if (v0 >= v1)
          { chX = (interval_t){{ X.end[0], t->xm }};  ch = t->ch[0]; }
        else
          { chX = (interval_t){{ t->xm, X.end[1] }};  ch = t->ch[1]; }
        intgprod_node_split(ch, &chX, n, est_f);
      }

    /* Update {t->IoF}: */
    t->IoF = interval_add(&(t->ch[0]->IoF), &(t->ch[1]->IoF));

    return;
      
    void intgprod_node_leaf_split(t, X, n, est_f)
      { /* Node {t} is a leaf node. */
        assert((t->ch[0] == NULL) && (t->ch[1] == NULL)); /* Paranoia. */
        assert(isnan(t->xm)); /* Paranoia. */

        /* Split {*X} into {chX[0..1]}: */
        t->xm = interval_mid(X);
        interval_t chX[2];
        chX[0] = (interval_t){{ X.end[0], xm }};
        chX[1] = (interval_t){{ xm, X.end[1] }};

        /* Create the two children of {t}: */
        int i, k;
        for (k = 0; k < 2; k++) { t->ch[k] = intgprod_node_new(); }

        /* Compute the graph boxes for each half: */
        { interval_t *Y = t->Y;
          t->Y = NULL; 
          t->ch[0]->Y = Y; /* Reuse old {t->Y} vector as the Y-range list of child 0. */
          t->ch[1]->Y = notnull(malloc(n * sizeof(interval_t)), "no mem");
          for (i = 0; i < t->nv; i++)
            { /* Split the box {X × Y[i]} into two boxes at {xm}: */
              interval_t Yi0, Yi1; /* Local temps, in case {est_f} is careless. */
              est_f(t->iv[i], X, &(Y[i]), t->xm, &Yi0, &Yi1); 
              t->ch[0]->Y[i] = Yi0;
              t->ch[1]->Y[i] = Yi1;
            }
        }

        /* Compute the children's {IoF} fields: */
        for (k = 0; k < 2; k++) 
          { t->ch[k].IoF = intgprod_compute_IoF(&(chX[k]), n, t->ch[k]->Y); }
      }
  }

interval_t intgprod_integrate(double a[])
  {
    demand(FALSE, "!!! not implemented yet !!!");
  }
