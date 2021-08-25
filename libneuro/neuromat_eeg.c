/* See {neuromat_eeg.h}. */
/* Last edited on 2021-08-23 23:19:22 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>

#include <fget.h>
#include <affirm.h>
#include <jsstring.h>
#include <neuromat_eeg.h>

void neuromat_eeg_get_channel_names(int32_t ne, int32_t nv, char **evname, int32_t *ncP, char ***chnamesP)
  { int32_t nc = 0; /* Number of electrodes plust trigger/reference channel. */
    char **chnames = NULL; /* Electrode names and trigger/reference channel name. */
    if (ne == 20)
      { neuromat_eeg_get_20_channel_names(&nc, &chnames); }
    else if (ne == 128)
      { neuromat_eeg_get_128_channel_names(&nc, &chnames); }
    else
      { demand(FALSE, "invalid electrode count"); }
    if (nv > 0)
      { /* Append the event channel names: */
        chnames = notnull(realloc(chnames, (nc+nv)*sizeof(char *)), "no mem");
        for (int32_t i = 0; i < nv; i++) 
          { chnames[nc+i] = txtcat(evname[i], ""); }
        nc = nc + nv;
      }
    /* Return results: */
    (*ncP) = nc; 
    (*chnamesP) = chnames;
  }

void neuromat_eeg_get_20_channel_names(int32_t *ncP, char ***chnamesP)
  { 
    int32_t nc = 21;
    char **chnames = notnull(malloc(nc*sizeof(char*)), "no mem");
    chnames[ 0] = "F7";
    chnames[ 1] = "T3";
    chnames[ 2] = "T5";
    chnames[ 3] = "Fp1";
    chnames[ 4] = "F3";
    chnames[ 5] = "C3";
    chnames[ 6] = "P3";
    chnames[ 7] = "O1";
    chnames[ 8] = "F8";
    chnames[ 9] = "T4";
    chnames[10] = "T6";
    chnames[11] = "Fp2";
    chnames[12] = "F4";
    chnames[13] = "C4";
    chnames[14] = "P4";
    chnames[15] = "O2";
    chnames[16] = "Fz";
    chnames[17] = "Cz";
    chnames[18] = "Pz";
    chnames[19] = "Oz"; 
    chnames[20] = "TR"; /* Trigger signal. */
    /* Make all strings be newly allocated: */
    for (int32_t ic = 0; ic < nc; ic++) { chnames[ic] = txtcat(chnames[ic], ""); }
    /* Return results: */
    (*ncP) = nc;
    (*chnamesP) = chnames;
  }
  
void neuromat_eeg_get_128_channel_names(int32_t *ncP, char ***chnamesP)  
  {
    int32_t ne = 128;
    int32_t nc = ne + 1;
    char **chnames = notnull(malloc(nc*sizeof(char *)), "no mem");
    /* Electrode channels: */
    for (int32_t i = 0; i < ne; i++) 
      { char *name = NULL;
        asprintf(&name, "C%03d", i+1);
        chnames[i] = name;
      }
    /* Reference electrode: */
    chnames[ne] = txtcat("CZ", ""); /* To get a newly allocated copy. */
    /* Return results: */
    (*ncP) = nc;
    (*chnamesP) = chnames;
  }

int32_t neuromat_eeg_find_channel_by_name(char *name, int32_t ic_start, int32_t ic_end, char *chnames[], bool_t die)
  {
    int32_t ict = ic_start;
    while ((ict <= ic_end) && (strcmp(name, chnames[ict]) != 0)) { ict++; }
    if (ict > ic_end)
      { if (die)
          { fprintf(stderr, "looking for channel \"%s\" in channels [%d..%d] = {", name, ic_start, ic_end);
            for (int32_t jc = ic_start; jc <= ic_end; jc++) { fprintf(stderr, " %s", chnames[jc]); }
            fprintf(stderr, " }\n");
            demand(FALSE, "no such trigger channel");
          }
        else
          { ict = -1; }
      }
    return ict;
  }
    
int32_t neuromat_eeg_locate_pulses
  ( int32_t nt,              /* Number of data frames. */
    int32_t nc,              /* Number of data channels. */
    double **val,        /* The EEG dataset ({nt} frames with {nc} channels each). */
    int32_t ict,             /* Index of trigger channel. */
    double vmin,         /* "Low" trigger channel value. */
    double vmax,         /* "High" trigger channel value. */
    int32_t np_max,          /* Max expected number of pulses in file. */
    int32_t it_pulse_ini[],  /* (OUT) Index of first frame in each pulse. */
    int32_t it_pulse_fin[]   /* (OUT) Index of last frame in each pulse. */
  )
  {
    demand((ict) && (ict < nc), "invalid trigger channel index");
    double vtr_prev = 0; /* Trigger channel value in previous frame. */
    int32_t np = 0; /* Number of pulses found. */
    for (int32_t it = 0; it <= nt; it++)
      { double vtr = (it < nt ? val[it][ict] : 0.0);
        demand((vtr == vmin) || (vtr == vmax), "invalid trigger value");
        if ((vtr_prev == vmin) && (vtr == vmax))
          { /* Trigger up-event: */
            demand(it > 0, "incomplete trigger pulse at start of file");
            demand(np < np_max, "too many trigger pulses");
            it_pulse_ini[np] = it;
            np++;
          }
        else if ((vtr_prev == vmax) && (vtr == vmin))
          { /* Trigger down-event: */
            demand(it < nt, "incomplete trigger pulse at end of file");
            assert(np > 0);
            assert(it > 0);
            it_pulse_fin[np - 1] = it - 1;
          }
        vtr_prev = vtr;
      }
    return np;
  }

void neuromat_eeg_report_pulse(FILE *wr, char *pre, int32_t ic, char *name, int32_t it_ini, int32_t it_fin, double fsmp, char *suf)
  {
    if (pre != NULL) { fprintf(wr, "%s", pre); }
    fprintf(wr, "pulse in channel %d = \"%s\"", ic, name);
    fprintf(wr, " spanning frames %d..%d", it_ini, it_fin);
    fprintf(wr, " (%.4f _ %.4f sec)", ((double)it_ini)/fsmp, ((double)it_fin)/fsmp);
    if (suf != NULL) { fprintf(wr, "%s", suf); }
    fflush(wr);
  }
