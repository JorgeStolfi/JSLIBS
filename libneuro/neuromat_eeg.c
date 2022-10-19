/* See {neuromat_eeg.h}. */
/* Last edited on 2021-08-28 02:41:22 by stolfi */

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

void neuromat_eeg_get_20_channel_names(int32_t ne, char *chname[]);
  /* Returns in {chname[0..ne-1]} the names of the raw electrode channels in the 20-electrode experiment.
    They are assumed to be the 20 eelectrodes followed by {nv} extra marker channels, called {evname[0..nv-1]}. 
    The size of {channels} must be at least {ne}.
    
    The 20 electrode names are as follows.
    Note that the C-language indices ("C" in the table below) range from
    0 to 19 while the Fortran/Matlab indices ("F") range from 1 to 20.
    
      |    C  F  Name    C  F  Name    C  F  Name    C  F  Name  
      |   -- -- -----   -- -- -----   -- -- -----   -- -- -----
      |    0  1  "F7"    5  6 "C3"    10 11 "T6"    15 16 "O2"
      |    1  2  "T3"    6  7 "P3"    11 12 "Fp2"   16 17 "Fz"
      |    2  3  "T5"    7  8 "O1"    12 13 "F4"    17 18 "Cz"
      |    3  4  "Fp1"   8  9 "F8"    13 14 "C4"    18 19 "Pz"
      |    4  5  "F3"    9 10 "T4"    14 15 "P4"    19 20 "Oz"
      
    A 20-electrode data file usually has at least one marker channel, 
    a trigger channel with C-index 20 (F-index 21) and name "TR". 
    
    All the strings {chname[0..ne-1]} are newly allocated by the procedure.
  */

void neuromat_eeg_get_128_channel_names(int32_t ne, char *chname[]);
  /* Returns in {*chname[0..ne-1]} the names of the electrode channels channels for {ne}-electrode EEG datasets, 
    where {ne} is 128 or 129.
    The size of {channels} must be at least {ne}.
    
    Specifically, names {chname[0..127]} are "C001" to "C128" (ie.
    {chname[i]} is "C{i+1}"), for {i} in {0..127}. If {ne} is 129, the
    extra channel is assumed to be the reference (ground) electrode, so
    {chname[128] is set to "CZ" which is (name used by the Neuromat
    team; labeled "VREF" in the electrode diagram).
    
    All the strings {chname[0..ne-1]} are newly allocated by the procedure. */

char **neuromat_eeg_get_channel_names(int32_t ne, int32_t nv, char *evname[])
  { int32_t nc = ne + nv; /* Number of electrodes plus marker channels. */
    char **chname = notnull(malloc(nc*sizeof(char*)), "no mem"); /* Electrode names and trigger/reference channel name. */
    if (ne == 20)
      { neuromat_eeg_get_20_channel_names(ne, chname); }
    else if ((ne == 128) || (ne == 129))
      { neuromat_eeg_get_128_channel_names(ne, chname); }
    else
      { demand(FALSE, "invalid electrode count"); }
    if (nv > 0)
      { /* Append the trigger channel names: */
        for (int32_t iv = 0; iv < nv; iv++) { chname[ne+iv] = txtcat(evname[iv], ""); }
      }
    /* Return results: */
    return chname;
  }

void neuromat_eeg_get_20_channel_names(int32_t ne, char *chname[])
  { 
    demand(ne == 20, "invalid {ne}");
    chname[ 0] = "F7";
    chname[ 1] = "T3";
    chname[ 2] = "T5";
    chname[ 3] = "Fp1";
    chname[ 4] = "F3";
    chname[ 5] = "C3";
    chname[ 6] = "P3";
    chname[ 7] = "O1";
    chname[ 8] = "F8";
    chname[ 9] = "T4";
    chname[10] = "T6";
    chname[11] = "Fp2";
    chname[12] = "F4";
    chname[13] = "C4";
    chname[14] = "P4";
    chname[15] = "O2";
    chname[16] = "Fz";
    chname[17] = "Cz";
    chname[18] = "Pz";
    chname[19] = "Oz"; 
    /* Make all strings be newly allocated: */
    for (int32_t ie = 0; ie < ne; ie++) { chname[ie] = txtcat(chname[ie], ""); }
  }
  
void neuromat_eeg_get_128_channel_names(int32_t ne, char *chname[])  
  {
    demand((ne == 128) || (ne == 129), "invalid {ne}");
    /* Electrode channels: */
    for (int32_t ie = 0; ie < 128; ie++) 
      { char *name = NULL;
        asprintf(&name, "C%03d", ie+1);
        chname[ie] = name;
      }
    if (ne == 129)
      { /* Reference electrode: */
        chname[128] = txtcat("CZ", ""); /* To get a newly allocated copy. */
      }
  }

int32_t neuromat_eeg_find_channel_by_name(char *name, int32_t ic_start, int32_t ic_end, char *chname[], bool_t die)
  {
    int32_t ict = ic_start;
    while ((ict <= ic_end) && (strcmp(name, chname[ict]) != 0)) { ict++; }
    if (ict > ic_end)
      { if (die)
          { fprintf(stderr, "looking for channel \"%s\" in channels [%d..%d] = {", name, ic_start, ic_end);
            for (int32_t jc = ic_start; jc <= ic_end; jc++) { fprintf(stderr, " %s", chname[jc]); }
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
