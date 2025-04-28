/* See {psgr_shrink_stacy.h} */
/* Last edited on 2025-04-27 02:57:12 by stolfi */

/* Created by Rafael F. V. Saracchini */

/* See "A topological multi-grid method for integration of two-dimensional meshes"
by Rafael F. V. Saracchini and Jorge Stolfi. 
Unpublished draft, article-d.tex, 2016-01-08.  From page 11:

Table 2. Formulas for the weight w'_{01} of the new edge e'_{01} = (v0,v1) created by the star-cycle swap,
for each degree k. The same formulas are valid for each other edge e'_{i,i+1} except by the arc indices which
are increased by i modulo k.

k   w_{01}
2   w_0*w_1/w_t
3   0.5*(w_0*w_1 + w1*w2)/w_t
4   (w_0*w_1 + 0.5*(w_0*w_2 + w_1*w_3))/w_t
5   (w_0*w_1 + 1.169*(w_2*w_4 + w_0*w_2 + w_1*w_4))/w_t
6   (w_0*w_1 + 2*w_5*w_2 + 1.5*(w_5*w_1 + w_0*w_2))/w_t

*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <float_image.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>
#include <rn.h>

#include <psgr_types.h>
#include <psgr.h>
#include <psgr_copy.h>
#include <psgr_shrink.h>

#include <psgr_shrink_stacy.h>

/* SHORTER NAMES */  

#define VNONE psgr_vertex_NONE
#define ENONE psgr_star_NONE
#define ANONE psgr_arc_NONE

#define CLEAR   psgr_mark_CLEAR
#define KILL    psgr_mark_KILL   
#define KEEP    psgr_mark_KEEP   
#define DELETED psgr_mark_DELETED
    
void psgr_shrink_stacy_debug_star_data(psgr_t *gr, psgr_arc_t a);

double psgr_shrink_stacy_weight_3
  ( psgr_t* gr, 
    uint32_t n,
    uint32_t k0,
    double d_star[],
    double w_star[]
  );

double psgr_shrink_stacy_weight_4
  ( psgr_t* gr, 
    uint32_t n,
    uint32_t k0,
    double d_star[],
    double w_star[]
  );

double psgr_shrink_stacy_weight_5
  ( psgr_t* gr, 
    uint32_t n,
    uint32_t k0,
    double d_star[],
    double w_star[] 
  );

double psgr_shrink_stacy_weight_6
  ( psgr_t* gr, 
    uint32_t n,
    uint32_t k0,
    double d_star[],
    double w_star[]
  );
  
/* IMPLEMENTATIONS */
    
void psgr_shrink_stacy_get_star_data
  ( psgr_t *gr,
    psgr_arc_t a,
    uint32_t n,
    double d_star[],
    double w_star[],
    double *wtot_P
  )
  { 
    psgr_arc_t b = a;
    double wtot = 0;
    for (int32_t k = 0; k < n; k++)
      { double delta = psgr_arc_delta(gr, b);
        d_star[k] = delta;

        double weight = psgr_arc_weight(gr, b);
        demand(isfinite(weight) && (weight > 0), "arc with weight");
        w_star[k] = weight;
        wtot += weight;
        
        b = psgr_arc_onext(gr, b);
      }

    demand(isfinite(wtot), "overflow in {wtot}");
    assert(wtot > 0);
    (*wtot_P) = wtot;
  }

void psgr_shrink_stacy_compute_cycle_data
  ( psgr_t *gr,
    uint32_t n,
    double d_star[],
    double w_star[],
    double wtot,
    double d_cycle[],
    double w_cycle[]
  )
  {
    demand(isfinite(wtot) && (wtot > 0), "invalid {wtot}");
      
    for (uint32_t k0 = 0; k0 < n; k0++)
      { uint32_t k1 = (k0 + 1) % n;
        d_cycle[k0] = d_star[k1] - d_star[k0];
        
        demand(isfinite(w_star[k0]) && (w_star[k0] > 0), "invalid star edge weight");
        switch (n)
          { case 3: w_cycle[k0] = psgr_shrink_stacy_weight_3(gr, k0, n, w_star, d_star); break;
            case 4: w_cycle[k0] = psgr_shrink_stacy_weight_4(gr, k0, n, w_star, d_star); break;
            case 5: w_cycle[k0] = psgr_shrink_stacy_weight_5(gr, k0, n, w_star, d_star); break;
            case 6: w_cycle[k0] = psgr_shrink_stacy_weight_6(gr, k0, n, w_star, d_star); break;
            default: assert(FALSE);
          }
        w_cycle[k0] = w_cycle[k0]/wtot;
        /* The computed weight may be zero because of underflow: */
        assert(isfinite(w_cycle[k0]) && (w_cycle[k0] >= 0));
      }
  }

double psgr_shrink_stacy_weight_3
  ( psgr_t* gr, 
    uint32_t n,
    uint32_t k0,
    double d_star[],
    double w_star[]
  )
  { assert(n == 3);
    
    uint32_t k1 = (k0 + 1) % n;
    uint32_t k2 = (k0 + 2) % n;

    double weight = 0.5*(w_star[k0]*w_star[k1] + w_star[k1]*w_star[k2]);
    return weight;
  }

double psgr_shrink_stacy_weight_4
  ( psgr_t* gr, 
    uint32_t n,
    uint32_t k0,
    double d_star[],
    double w_star[]
  )
  { assert(n == 4);
  
    uint32_t k1 = (k0 + 1) % n;
    uint32_t k2 = (k0 + 2) % n;
    uint32_t k3 = (k0 + 3) % n;

    double weight = w_star[k0]*w_star[k1] + 0.5*(w_star[k0]*w_star[k2] + w_star[k1]*w_star[k3]);

    return weight;
  }

double psgr_shrink_stacy_weight_5
  ( psgr_t* gr, 
    uint32_t n,
    uint32_t k0,
    double d_star[],
    double w_star[]
  )
  { assert(n == 5);
  
    double A = 1.1690;
  
    uint32_t k1 = (k0 + 1) % n;
    uint32_t k2 = (k0 + 2) % n;
    uint32_t k4 = (k0 + 4) % n;

    double weight = w_star[k0]*w_star[k1];
    weight += A*(w_star[k2]*w_star[k4] + w_star[k0]*w_star[k2] + w_star[k1]*w_star[k4]);

    return weight;
  }

double psgr_shrink_stacy_weight_6
  ( psgr_t* gr, 
    uint32_t n,
    uint32_t k0,
    double d_star[],
    double w_star[]
  )
  {
    double A = 2.0;
    double C = 1.5;

    uint32_t k1 = (k0 + 1) % n;
    uint32_t k2 = (k0 + 2) % n;
    uint32_t k5 = (k0 + 5) % n;

    double weight = w_star[k0]*w_star[k1];
    weight += A*(w_star[k5]*w_star[k2]);
    weight += C*(w_star[k5]*w_star[k1] + w_star[k0]*w_star[k2]);

    return weight;
  }
