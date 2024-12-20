/* See {btc_bubble_parms_validate.h} */
/* Last edited on 2024-12-05 10:23:13 by stolfi */

#include <stdio.h>
#include <bool.h>

#include <btc_bubble_t.h>


#include <btc_bubble_parms_validate.h>

bool_t btc_bubble_parms_validate(char* fName, int nlin, int nd, char* dt[], btc_bubble_t* bpj)
  {
    bool_t ok = TRUE;

    int id_fin_up = bpj->id_fin_up; /* Date of end of rally, start of plateau: */
    int id_ini_dn = bpj->id_fin_up + bpj->wd_plat; /* Date of end of plateau, start of decay. */
    
    if (bpj->id_ini_sg >= id_fin_up)
      { fprintf(stderr, "%s:%d: !! warning: start-relevant date after peak %s %s", fName, nlin, dt[bpj->id_ini_sg], dt[id_fin_up]); }
    if (id_fin_up > id_ini_dn)
      { fprintf(stderr, "%s:%d: ** invalid date interval %s %s", fName, nlin, dt[id_fin_up], dt[id_ini_dn]);
        ok = FALSE;
      }
    if (id_ini_dn >= bpj->id_fin_sg)
      { fprintf(stderr, "%s:%d: !! warning: end-relevant date before peak %s %s", fName, nlin, dt[id_ini_dn], dt[bpj->id_fin_sg]); }


    if (bpj->rt_up < 1.0)
      { fprintf(stderr, "%s:%d: ** growth rate %8.6f must be at least 1.0", fName, nlin, bpj->rt_up);
        ok = FALSE;
      }
    if (bpj->rt_up > 10.0)
      { fprintf(stderr, "%s:%d: ** growth rate %8.6f too big", fName, nlin, bpj->rt_up);
        ok = FALSE;
      }
    if (bpj->rt_dn > 1.1)
      { fprintf(stderr, "%s:%d: ** decay rate %8.6f must be at most 1.1", fName, nlin, bpj->rt_dn);
        ok = FALSE;
      }
    if (bpj->rt_dn < 0.1)
      { fprintf(stderr, "%s:%d: ** decay rate %8.6f too small", fName, nlin, bpj->rt_dn);
        ok = FALSE;
      }
    if (((bpj->rt_up == 1.0) || (bpj->rt_dn == 1.0)) && (id_fin_up != id_ini_dn))
      { fprintf(stderr, "%s:%d: ** invalid dates for degenerate bubble %s %s", fName, nlin, dt[id_fin_up], dt[id_ini_dn]);
        ok = FALSE;
      }
      
    return ok;
  }

