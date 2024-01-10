/* See {neuromat_eeg_frame_buffer.oh}. */
/* Last edited on 2013-11-13 19:28:08 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_io.h>
#include <neuromat_eeg_raw_io.h>

#include <neuromat_eeg_frame_buffer.h>
      
neuromat_eeg_frame_buffer_t *neuromat_eeg_frame_buffer_new(int size, int nc)
  {
    neuromat_eeg_frame_buffer_t *buf = notnull(malloc(sizeof(neuromat_eeg_frame_buffer_t)), "no mem");
    
    (*buf) = (neuromat_eeg_frame_buffer_t)
      { .size = size,
        .val = notnull(malloc(size*sizeof(double *)), "no mem"),
        .nc = nc,    /* Number of channels */
        .it_ini = 0, /* Index of first frame in buffer. */
        .it_fin = -1 /* Index of last frame in buffer. */
      };
    int itb;
    for (itb = 0; itb < size; itb++) { buf->val[itb] = NULL; }
    return buf;
  }

void neuromat_eeg_frame_buffer_free(neuromat_eeg_frame_buffer_t *buf)
  {
    int itb;
    for (itb = 0; itb < buf->size; itb++) { free(buf->val[itb]); }
    free(buf->val);
    free(buf);
  }

int neuromat_eeg_frame_buffer_get_frame
  ( neuromat_eeg_frame_buffer_t *buf,
    int it,
    neuromat_eeg_frame_buffer_read_proc_t read_frame
  )
  {
    bool_t debug = FALSE;
    demand (it >= buf->it_ini, "requested data frame was flushed and cannot be re-read");
    assert(buf->it_ini <= buf->it_fin + 1);
    while (it > buf->it_fin)
      { /* Flush the first frame, if needed, to make space for one more frame: */
        int nb = buf->it_fin + 1 - buf->it_ini;
        if (nb >= buf->size) { (buf->it_ini)++; }
        /* Get the index of the next frame in file {jt} and in the buffer {jt_buf}: */
        int jt = buf->it_fin + 1;
        int jt_buf = (jt % buf->size);
        /* Get or allocate storage for the next frame: */
        if (debug) { fprintf(stderr, "  reading into row %d of buffer\n", jt_buf); }
        double *frm = buf->val[jt_buf];
        if (frm == NULL) { frm = notnull(malloc(buf->nc*sizeof(double)), "no mem"); }
        /* Read the next frame into buffer: */
        int nrd = read_frame(buf->nc, frm);
        if (nrd == 0) 
          { /* Hit end of file: */
            if (frm != NULL) { free(frm); frm = NULL; }
            break;
          }
        /* Read succedded: */
        assert(nrd == 1);
        buf->val[jt_buf] = frm;
        (buf->it_fin)++;
      }
    assert(buf->it_ini <= buf->it_fin + 1); /* Paranoia. */
    assert(buf->it_ini <= it); /* Start of queue cannot have surpassed {it}. */

    if (it <= buf->it_fin)
      { /* We succeeded: */
        return (it % buf->size);
      }
    else
      { /* File must have ended before the desired record: */
        return -1;
      }
  }
      
void neuromat_eeg_frame_buffer_find_next_pulse
  ( neuromat_eeg_frame_buffer_t *buf, 
    neuromat_eeg_frame_buffer_read_proc_t *read_frame,
    int it_start, 
    int ns,
    int ichs[], 
    int *it_iniP, 
    int *it_finP, 
    int *sP
  )
  {
    /* Data for a pulse: */
    int it_ini = -1, it_fin = -1; /* First and last frame of pulse. */
    int s = -1;  /* Phase type: stimulus type index {0..ns-1} or -1 if fixation. */ 
    
    int it_next = it_start;
    neuromat_eeg_frame_buffer_find_next_pulse_start(buf, read_frame, it_next, ns, ichs, &it_ini, &s);
    if (it_ini >= 0) 
      { assert(it_ini >= it_next);
        assert((s >= 0) && (s < ns));
        int ic = ichs[s]; /* Channel where pulse actually occurred. */
        it_next = it_ini;
        /* Skip rest of pulse: */
        neuromat_eeg_frame_buffer_find_pulse_end(buf, read_frame, it_next, ic, &it_fin);
        assert(it_fin >= it_ini);
      }
    (*it_iniP) = it_ini;
    (*it_finP) = it_fin;
    (*sP) = s;
  }

void neuromat_eeg_frame_buffer_find_next_pulse_start
  ( neuromat_eeg_frame_buffer_t *buf, 
    neuromat_eeg_frame_buffer_read_proc_t *read_frame,
    int it_start, 
    int ns, 
    int ichs[], 
    int *itP,
    int *sP
  )
  {
    int it = it_start;
    int s = -1; /* Channel that was up. */
    while (TRUE)
      { /* Make sure that frame {it} is in the buffer, grab it: */
        int itb = neuromat_eeg_frame_buffer_get_frame(buf, it, read_frame);
        if (itb < 0) 
          { /* End of file before next frame. */
            it = -1; 
            break;
          }
        double *vt = buf->val[itb];
        /* Check if any of the channels are up, save in {*sP}: */
        int nup = 0; /* Number of marker channels that are up. */ 
        int r;
        for (r = 0; r < ns; r++)
          { /* Check channel {ichs[r]}: */
            int ic = ichs[r];
            demand(vt[ic] >= 0, "stimulus phase marker channel is negative");
            if (vt[ic] > 0) 
              { if (s == -1) { s = r; }
                nup++;
              }
          }
        if (nup > 1) { fprintf(stderr, "!! warning: overlapping pulses"); }
        if (nup > 0) { break; }
        it++;
      }
    (*itP) = it; 
    (*sP) = s;
    
    return;
  }
    
void neuromat_eeg_frame_buffer_find_pulse_end
  ( neuromat_eeg_frame_buffer_t *buf, 
    neuromat_eeg_frame_buffer_read_proc_t *read_frame, 
    int it_start, 
    int ic, 
    int *itP
  )
  {
    int it = it_start;
    while(TRUE)
      { /* Make sure that frame {it} is in the buffer: */
        int itb = neuromat_eeg_frame_buffer_get_frame(buf, it, read_frame);
        if (itb < 0) 
          { /* End of file before next frame. */
            demand(it > it_start, "start frame does not exist");
            break;
          }
        double *vt = buf->val[itb];
        demand(vt[ic] >= 0, "phase marker channel is negative");
        if (vt[ic] == 0) 
          { demand(it > it_start, "pulse is off in initial frame");
            break;
          }
        it++;
      }
    (*itP) = it-1; 
  }

