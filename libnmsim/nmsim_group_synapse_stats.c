/* See {nmsim_group_synapse_stats.h} */
/* Last edited on 2020-12-14 10:05:31 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <sign.h>

#include <nmsim_basic.h>

#include <nmsim_group_synapse_stats.h>

void nmsim_group_synapse_stats_print_gen(FILE *wr, nmsim_group_synapse_stats_t *st, int32_t hdr)
  {
    char *fmt_hd1 = "%5s %8s %11s %17s";
    char *fmt_hd2 = "%5s %8s %11s %17s";
    char *fmt_hd3 = "%5s %8s %5s %5s %8s %8s";
    char *fmt_val = "%5d %8d %5d %5d %8d %8d";
    switch (hdr) {
      case 0: 
        { fprintf
            ( wr, fmt_val, 
              st->isg, st->nse, 
              st->ing_pre, st->ing_pos, 
              st->nne_pre, st->nne_pos
            );
        }
        break;
      case 1:
        fprintf(wr, fmt_hd1, "", "", "neur groups", "num neurons");
        break;
      case 2:
        fprintf(wr, fmt_hd2, "syn", "total", "-----------", "-----------------");
        break;
      case 3:
        fprintf(wr, fmt_hd3, "group", "synapses", "pre", "post", "pre", "post");
        break;
      case 4:
        fprintf(wr, fmt_hd3, "-----", "--------", "-----", "-----", "--------", "--------");
        break;
      default:
        assert(FALSE);
    }
  }

void nmsim_group_synapse_stats_print_sgn(FILE *wr, nmsim_group_synapse_stats_t *st, sign_t sgn, int32_t hdr)
  {
    /* Format fields mea: */
    char *fmt_hd1 = "%57s";
    char *fmt_hd2 = "%57s";
    char *fmt_hd3 = "%8s %8s %8s %8s %10s %10s";
    char *fmt_val = "%8d %+8.4f %8.4f %8.1f %+10.4f %10.4f";
    switch (hdr) {
      case 0: 
        { double navg = ((double)st->nse_cons)/((double)st->nne_pos); /* Avg good inputs per post-syn neuron. */
          fprintf(wr, fmt_val, st->nse_cons, st->Wavg, st->Wdev, navg, st->Wina, st->Wind); 
        }
        break;
      case 1:
        { char *stype = (sgn < 0 ? "inhibitory synapses" : "excitatory synapses");
          fprintf(wr, fmt_hd1, stype);
        }
        break;
      case 2:
        fprintf(wr, fmt_hd2, "---------------------------------------------------------");
        break;
      case 3:
        fprintf(wr, fmt_hd3, "count", "Wavg", "Wdev", "navg", "Wina", "Wind");
        break;
      case 4:
        fprintf(wr, fmt_hd3, "--------", "--------", "--------", "--------", "----------", "----------");
        break;
      default:
        assert(FALSE);
    }
  }

void nmsim_group_synapse_stats_print_uns(FILE *wr, nmsim_group_synapse_stats_t *st, int32_t hdr)
  {
    char *fmt_hd1 = "%37s";
    char *fmt_hd2 = "%37s";
    char *fmt_hd3 = "%8s %8s %8s %10s";
    char *fmt_val = "%8d %+8.4f %8.1f %+10.4f";
    
    switch (hdr) {
      case 0: 
        fprintf(wr, fmt_val, st->nse_cons, st->Wavg, st->navg, st->Wina);
        break;
      case 1:
        fprintf(wr, fmt_hd1, "nonzero synapses");
        break;
      case 2:
        fprintf(wr, fmt_hd2, "-------------------------------------");
        break;
      case 3:
        fprintf(wr, fmt_hd3, "count", "Wavg", "navg", "Wina");
        break;
      case 4:
        fprintf(wr, fmt_hd3, "--------", "--------", "--------", "----------");
        break;
      default:
        assert(FALSE);
    }
  }
