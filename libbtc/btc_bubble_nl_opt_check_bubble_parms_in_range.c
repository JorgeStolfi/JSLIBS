/* See {btc_bubble_nl_opt_check_bubble_parms_in_range.h} */
/* Last edited on 2024-12-21 11:58:24 by stolfi */

#include <stdlib.h>

#include <bool.h>
#include <affirm.h>

#include <btc_bubble_t.h>
#include <btc_is_trivial_int_range.h>
#include <btc_is_in_int_range.h>
#include <btc_is_in_double_range.h>

#include <btc_bubble_nl_opt_check_bubble_parms_in_range.h>

void btc_bubble_nl_opt_check_bubble_parms_in_range
  ( int nb, 
    btc_bubble_t bp_lo[], 
    btc_bubble_t bp[], 
    btc_bubble_t bp_hi[]
  )
  { 
    demand((bp_lo != NULL) && (bp_hi != NULL), "min and max parameter sets must be given");

    int ib;
    for (ib = 0; ib < nb; ib++)
      { btc_bubble_t* blo = &(bp_lo[ib]);
        btc_bubble_t* b   = &(bp[ib]);
        btc_bubble_t* bhi = &(bp_hi[ib]);
        
        (void)btc_is_trivial_int_range(blo->id_ini_sg, b->id_ini_sg, bhi->id_ini_sg, ib, ".id_ini_sg", TRUE);
        (void)btc_is_in_int_range(blo->id_fin_up, b->id_fin_up, bhi->id_fin_up, ib, ".id_fin_up", TRUE);
        (void)btc_is_in_double_range(blo->rt_up, b->rt_up, bhi->rt_up, ib, ".rt_up", TRUE);
        (void)btc_is_in_int_range(blo->wd_plat, b->wd_plat, bhi->wd_plat, ib, ".wd_plat", TRUE);
        (void)btc_is_in_double_range(blo->rt_dn, b->rt_dn, bhi->rt_dn, ib, ".rt_dn", TRUE);
        (void)btc_is_trivial_int_range(blo->id_fin_sg, b->id_fin_sg, bhi->id_fin_sg, ib, ".id_fin_sg", TRUE);
      }
  }
    


