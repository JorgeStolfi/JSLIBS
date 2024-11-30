/* See {neuromat_eeg.h}. */
/* Last edited on 2023-11-22 10:22:32 by stolfi */

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

double **neuromat_eeg_new(int32_t nt, int32_t nc)
  { double **val = talloc(nt, double*);
    for (uint32_t it = 0;  it < nt; it++) { val[it] = talloc(nc, double); }
    return val;
  }
  
void neuromat_eeg_free(double **val, int32_t nt, int32_t nc)
  { for (uint32_t it = 0;  it < nt; it++) { free(val[it]); }
    free(val);
  }

void neuromat_eeg_get_R20_channel_names(int32_t ne, char *chname[]);
  /* Returns in {chname[0..ne-1]} the names of the raw electrode names in the 20-electrode cap
    used at INDC sometime before 2013 by Ghislain et al. 
    The size of {channels} must be at least {ne} which must be 20.
    
    Set {chname[0..19]} as follows.
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
    
    All the strings {chname[0..ne-1]} are newly allocated by the procedure. */

void neuromat_eeg_get_R128_channel_names(int32_t ne, char *chname[]);
  /* Returns in {*chname[0..ne-1]} the electrode names for
    {ne}-electrode EEG cap used at INDC sometime before 2013,
    by Ghislain et al, Andressa et al.
    
    The size of {channels} must be at least {ne}, which must be 128 or 129.
    
    Sets {chname[0..127]} to "C001" to "C128" (ie. {chname[i]} is "C{i+1}"), 
    for {i} in {0..127}. If {ne} is 129, {chname[128]} is set to "CZ" which 
    is name used by the Neuromat team for the reference (ground) electrode 
    (labeled "VREF" in the electrode diagram).
    
    All the strings {chname[0..ne-1]} are newly allocated by the procedure. */

void neuromat_eeg_get_FN3_channel_names(int32_t ne, char *chname[]);
  /* Returns in {*chname[0..ne-1]} the electrode names for
    virtual electrodes used in the Renewal Points paper by Fernando
    Najman et al in 2023.
    
    The size of {channels} must be at least {ne}, which must be 3.
    The three channel names are "RPF", "LPF", and "OCC".
    
    All the strings {chname[0..ne-1]} are newly allocated by the procedure. */

void neuromat_eeg_get_FN128_channel_names(int32_t ne, char *chname[]);
  /* Returns in {*chname[0..ne-1]} the electrode names for
    the real raw electrodes used in the Renewal Points paper by Fernando
    Najman et al in 2023.
    
    The size of {channels} must be at least {ne}, which must be 128.
    The channel names will be "C001",.. "C128". 
    
    All the strings {chname[0..ne-1]} are newly allocated by the procedure. */

void neuromat_eeg_get_channel_names(char *capType, int32_t nv, char *evname[], int32_t *ne_P, char ***chname_P)
  { int32_t ne = 0;
    char **chname = NULL;
    for (uint32_t pass = 0;  pass < 2; pass++)
      { if (strcmp(capType, "R20") == 0)
          { if (pass == 0) { ne = 20; }  else { neuromat_eeg_get_R20_channel_names(ne, chname); } }
        else if (strcmp(capType, "R128") == 0)
          { if (pass == 0) { ne = 128; }  else { neuromat_eeg_get_R128_channel_names(ne, chname); } }
        else if (strcmp(capType, "R129") == 0)
          { if (pass == 0) { ne = 129; }  else { neuromat_eeg_get_R128_channel_names(ne, chname); } }
        else if (strcmp(capType, "FN3") == 0)
          { if (pass == 0) { ne = 3; }  else { neuromat_eeg_get_FN3_channel_names(ne, chname); } }
        else if (strcmp(capType, "FN128") == 0)
          { if (pass == 0) { ne = 128; }  else { neuromat_eeg_get_FN128_channel_names(ne, chname); } }
        else
          { demand(FALSE, "unknown cap type"); }
        if (pass == 0)
          { int32_t nc = ne + nv; /* Number of electrodes plus marker channels. */
            chname = talloc(nc, char*); /* Electrode names and trigger/reference channel name. */
          }
      }
    if (nv > 0)
      { /* Append the trigger channel names: */
        for (uint32_t iv = 0;  iv < nv; iv++) { chname[ne+iv] = txtcat(evname[iv], ""); }
      }
    /* Return results: */
    (*ne_P) = ne;
    (*chname_P) = chname;
  }

void neuromat_eeg_get_R20_channel_names(int32_t ne, char *chname[])
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
    for (uint32_t ie = 0;  ie < ne; ie++) { chname[ie] = txtcat(chname[ie], ""); }
  }
  
void neuromat_eeg_get_R128_channel_names(int32_t ne, char *chname[])  
  {
    demand((ne == 128) || (ne == 129), "invalid {ne}");
    /* Electrode channels: */
    for (uint32_t ie = 0;  ie < 128; ie++) 
      { char *name = NULL;
        char *name = jsprintf("C%03d", ie+1);
        chname[ie] = name;
      }
    if (ne == 129)
      { /* Reference electrode: */
        chname[128] = txtcat("CZ", ""); /* To get a newly allocated copy. */
      }
  }

void neuromat_eeg_get_FN3_channel_names(int32_t ne, char *chname[]) 
  {
    demand(ne == 3, "invalid {ne}");
    chname[0] = txtcat("RPF","");
    chname[1] = txtcat("LPF","");
    chname[2] = txtcat("OCC","");
  }

void neuromat_eeg_get_FN128_channel_names(int32_t ne, char *chname[])
  { demand(ne == 128, "invalid {ne}");
    for (uint32_t ke = 0;  ke < ne; ke++)
      { char *name = NULL;
        char *name = jsprintf("C%03d", ke+1);
        chname[ke] = name;
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
    for (uint32_t it = 0;  it <= nt; it++)
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
