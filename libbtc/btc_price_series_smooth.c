/* See {btc_price_series_smooth.h} */
/* Last edited on 2024-12-05 10:23:57 by stolfi */

#include <btc_price_series_smooth.h>
#include <btc_price_series_local_avg.h>


void btc_price_series_smooth(int nd, double vi[], int hrad, double vo[])
  {
    int id;
    for (id = 0; id <= nd-1; id++)
      { vo[id] = btc_price_series_local_avg(nd,vi,id,hrad); }
  }
