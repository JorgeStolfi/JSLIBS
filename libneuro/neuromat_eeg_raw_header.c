/* See {neuromat_eeg_raw_header.h}. */
/* Last edited on 2023-11-02 05:30:51 by stolfi */
  
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <jsmath.h>
#include <jsstring.h>

#include <neuromat_eeg_header.h>
#include <neuromat_eeg_raw_io.h>
#include <neuromat_eeg_raw_header.h>

neuromat_eeg_raw_header_t *neuromat_eeg_raw_header_read(FILE *rd)
  {
    neuromat_eeg_raw_header_t *hr = notnull(malloc(sizeof(neuromat_eeg_raw_header_t)), "no mem");
    
    /* Parse the fixed partof the header: */
    (*hr) = (neuromat_eeg_raw_header_t)
      { .version = neuromat_eeg_raw_read_int32(rd, "version", 2, 4),        /* Sample data type (2 = {int16_t}, 4 = {float}). */
        .year =    neuromat_eeg_raw_read_int16(rd, "year",  2000, 2100),    /* Recording time: year (4-digit). */
        .month =   neuromat_eeg_raw_read_int16(rd, "month", 1, 12),         /* Recording time: month. */
        .day =     neuromat_eeg_raw_read_int16(rd, "day",   1, 31),         /* Recording time: day. */
        .hour =    neuromat_eeg_raw_read_int16(rd, "hour",  0, 23),         /* Recording time: hour. */
        .min =     neuromat_eeg_raw_read_int16(rd, "min",   0, 59),         /* Recording time: minute. */
        .sec =     neuromat_eeg_raw_read_int16(rd, "sec",   0, 59),         /* Recording time: second. */
        .msec =    neuromat_eeg_raw_read_int32(rd, "msec",  0, 999),        /* Recording time: millisecond. */
        .rate =    neuromat_eeg_raw_read_int16(rd, "rate",  1, 10000),      /* Sampling rate (samples per second). */
        .nc =      neuromat_eeg_raw_read_int16(rd, "nc",    1, 999),        /* Number of channels. */
        .gain =    neuromat_eeg_raw_read_int16(rd, "gain",  1, 8),          /* Board gain (1,2,4,8). */
        .bits =    neuromat_eeg_raw_read_int16(rd, "bits",  0, 15),         /* Number of conversion bits. */
        .uVmax =   neuromat_eeg_raw_read_int16(rd, "uVmax", 0, 5000),       /* Full-scale range of amplifier in uV. */
        .nt =      neuromat_eeg_raw_read_int32(rd, "nt",    0, 1000000000), /* Number of samples. */
        .nv =      neuromat_eeg_raw_read_int16(rd, "nv",    0, 10000),      /* Number of event channels. */
        .evnames =  NULL                                                     /* Event channel names (for now) */
      };
    
    /* Parse the event channel names: */
    hr->evnames = notnull(malloc(hr->nv*sizeof(char *)), "no mem");
    for(int32_t ie = 0; ie < hr->nv; ie++)
      { hr->evnames[ie] = neuromat_eeg_raw_read_event_code(rd, "event"); }
    
    return hr;
  }

void neuromat_eeg_raw_header_print(FILE *wr, neuromat_eeg_raw_header_t *hr)      
  {
    fprintf(wr, "sample data type = %d\n", hr->version);
    fprintf(wr, "recording date = %04d-%02d-%02d\n", hr->year, hr->month, hr->day);
    fprintf(wr, "recording time = %02d:%02d:%02d.%03d\n", hr->hour, hr->min, hr->sec, hr->msec);
    fprintf(wr, "sampling rate = %d Hz\n", hr->rate);
    fprintf(wr, "number of channels = %d\n", hr->nc);
    fprintf(wr, "board gain = %d\n", hr->gain);
    fprintf(wr, "conversion bits = %d\n", hr->bits);    
    fprintf(wr, "amplifier range = %d uV\n", hr->uVmax);
    fprintf(wr, "number of data frames = %d\n", hr->nt);
    fprintf(wr, "number of event channels = %d\n", hr->nv);
    
    fprintf(wr, "event channel names =");
    for(int32_t ie = 0; ie < hr->nv; ie++) 
      { if ((ie > 0) && ((ie % 10) == 0)) { fprintf(wr, "\n "); }
        fprintf(wr, " %s", hr->evnames[ie]);
      }
    fprintf(wr, "\n");
    fflush(wr);
  }
  
neuromat_eeg_header_t *neuromat_eeg_raw_header_to_plain_header
  ( neuromat_eeg_raw_header_t *hr,
    char *file,
    int32_t skip,
    int32_t copy,
    int32_t subject,
    int32_t run
  )
  {
    int32_t nc = hr->nc;
    int32_t ne = nc - hr->nv;  /* ??? MAY BE WRONG ??? */
    /* Hack: guess the cap type from number of electrodes. */
    /* !!! Cap type should be a parameter. !!! */
    char *capType = NULL;
    if (ne == 20)
      { capType = "R20"; }
    else if ((ne == 128) || (ne == 129))
      { capType = "R128"; }
    else if (ne == 3)
      { capType = "FN3"; }
    else
      { demand(FALSE, "invalid electrode count"); }
    int32_t ntOrig = hr->nt; /* Number of frames in original file. */
    
    /* Get the complete channel count and names (including event channels): */
    int32_t ne_full = -1;
    char **chname = NULL;
    neuromat_eeg_get_channel_names(capType, hr->nv, hr->evnames, &ne_full, &chname);
    assert((ne == ne_full) || ((ne == 128) && (ne_full == 129)));
    
    /* Allocate the plain-format header: */
    neuromat_eeg_header_t *ht = neuromat_eeg_header_new();
    neuromat_eeg_source_t *orig = ht->orig; /* Source file and segment description. */
    
    (*ht) = (neuromat_eeg_header_t)
      {
        /* Dataset size and content: */
        .nt = copy,           /* Num of data frames. */
        .nc = nc,             /* Num of channels (incl. trigger, reference, event, etc.) */

        /* Nature of data: */
        .fsmp = (double)(hr->rate), /* Sampling rate. */
        .ne = ne,             /* Num of electrode channels. */
        .subject = subject,   /* Subject ID number. */
        .run = run,           /* Run rumber. */
        .type = NULL,         /* Type of run. */
        .chname = chname,     /* Channel names. */
        .kfmax = INT32_MIN,   /* Max frequency index (only if spectrum data). */
        .component = NULL,    /* Principal component name. */

        /* Filtering data: */
        .tdeg = INT32_MIN,     /* Degree of polynomial trend. */
        .tkeep = INT32_MIN,    /* Whether the polynomial trend was kept (1) or suppressed (0). */
        .flo0 = NAN,           /* Lower cutoff frequency. */
        .flo1 = NAN,           /* Nominal lowest preserved frequency. */
        .fhi1 = NAN,           /* Nominal highest preserved frequency. */
        .fhi0 = NAN,           /* Higher cutoff frequency. */
        .finvert = INT32_MIN,  /* Filter complementation: 0 for bandpass, 1 for bandkill. */
        
        /* Source file information: */
        .orig = orig           /* Source parameters. */
      };
      
    (*orig) = (neuromat_eeg_source_t)
      { .file = txtcat(file, ""),
        .nt = ntOrig,
        .it_ini = skip,
        .it_fin = skip + copy - 1,
        .fsmp = (double)hr->rate,
        .subject = subject,
        .run = run
      };

    return ht;
  }
