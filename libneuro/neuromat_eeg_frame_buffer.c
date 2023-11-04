/* See {neuromat_eeg_frame_buffer.oh}. */
/* Last edited on 2023-10-21 21:48:08 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_io.h>
#include <neuromat_eeg_raw_io.h>

#include <neuromat_eeg_frame_buffer.h>
      
neuromat_eeg_frame_buffer_t *neuromat_eeg_frame_buffer_new(int32_t size, int32_t nc)
  {
    neuromat_eeg_frame_buffer_t *buf = notnull(malloc(sizeof(neuromat_eeg_frame_buffer_t)), "no mem");
    
    (*buf) = (neuromat_eeg_frame_buffer_t)
      { .size = size,
        .val = notnull(malloc(size*sizeof(double *)), "no mem"),
        .nc = nc,    /* Number of channels */
        .it_ini = 0, /* Index of first frame in buffer. */
        .it_fin = -1 /* Index of last frame in buffer. */
      };
    for (int32_t itb = 0; itb < size; itb++) { buf->val[itb] = NULL; }
    return buf;
  }

void neuromat_eeg_frame_buffer_free(neuromat_eeg_frame_buffer_t *buf)
  {
    for (int32_t itb = 0; itb < buf->size; itb++) { free(buf->val[itb]); }
    free(buf->val);
    free(buf);
  }

int32_t neuromat_eeg_frame_buffer_get_frame
  ( neuromat_eeg_frame_buffer_t *buf,
    int32_t it,
    neuromat_eeg_frame_buffer_read_proc_t read_frame,
    bool_t mirror
  )
  {
    bool_t debug = FALSE;
    
    int32_t it_raw = it; /* For messages. */
    int32_t it_end = -1; /* Index of last frame in file, or {-1} if not known. */

    /* Mirror negative indices about frame zero: */
    if (it < 0) 
      { if (mirror) 
          { it = -it; }
        else
          { demand(FALSE, "invalid negative frame index"); }
      }
    
    /* Make sure that the file has been read at least up to frame {it}: */
    /* Note that {it} may get reduced by mirroring about the last frame. */
    while (it > buf->it_fin)
      { /* Flush the first frame, if needed, to make space for one more frame: */
        assert(buf->it_ini <= buf->it_fin + 1);
        int32_t nb = buf->it_fin + 1 - buf->it_ini;
        if (nb >= buf->size) { (buf->it_ini)++; }
        /* Get the index of the next frame in file {jt} and in the buffer {jt_buf}: */
        int32_t jt = buf->it_fin + 1;
        int32_t jt_buf = (jt % buf->size);
        /* Get or allocate storage for the next frame: */
        if (debug) { fprintf(stderr, "  reading into row %d of buffer\n", jt_buf); }
        double *frm = buf->val[jt_buf];
        if (frm == NULL) { frm = notnull(malloc(buf->nc*sizeof(double)), "no mem"); }
        /* Read the next frame into buffer: */
        int32_t nrd = (it_end >= 0 ? 0 : read_frame(buf->nc, frm));
        if (nrd == 1) 
          { assert(it_end < 0); /* We must not have hit EOF before. */
            /* Read succedded, store in buffer: */
            buf->val[jt_buf] = frm;
            (buf->it_fin)++;
            /* Loop back and test if {it} is in the buffered range. */
          }
        else
          { assert(nrd == 0);
            /* Hit end of file: */
            if (buf->it_fin < 0) 
              { /* File is empty, fail: */
                if ((frm != NULL) && (frm != buf->val[jt_buf])) { free(frm); frm = NULL; }
                fprintf(stderr, "!! file is empty\n");
                return -1;
              }
            else
              { /* The last valid frame index is {buf->it_fin}. */
                it_end = buf->it_fin;
                if (mirror)
                  { /* Mirror {it} about last frame: */
                    int32_t nt = it_end + 1; /* Number of frames in file. */
                    /* Execute multiple mirrors about end and zero: */
                    it = it % (2*nt); 
                    /* Mirror one last time about end if needed: */
                    if (it >= nt) { it = 2*it_end - it; }
                    assert((it >= 0) && (it <= it_end));
                    /* Loop back and test whether the new {it} has been read. */
                  }
                else
                  { /* No such frame: */
                    if ((frm != NULL) && (frm != buf->val[jt_buf])) { free(frm); frm = NULL; }
                    return -1;
                  }
              }
          }
      }
    
    /* If we got here, we must have at least one frame in the buffer: */
    assert(buf->it_ini <= buf->it_fin);
    
    if (it != it_raw)
      { /* Frame index was mirrored, warn user: */
        fprintf(stderr, "!! frame index %d folded to %d", it_raw, it);
        if (it_end >= 0) { fprintf(stderr, " to be in 0..%d", it_end); }
        fprintf(stderr, "\n");
        assert(mirror);
      }

    /* Check if desired frame (after mirroring) has been flushed: */
    demand (it >= buf->it_ini, "requested data frame was flushed and cannot be re-read");
    assert(it <= buf->it_fin); /* Because of mirroring. */

    return (it % buf->size);
  }
      
void neuromat_eeg_frame_buffer_find_next_pulse
  ( neuromat_eeg_frame_buffer_t *buf, 
    neuromat_eeg_frame_buffer_read_proc_t *read_frame,
    int32_t it_start, 
    int32_t ns,
    int32_t ichs[], 
    int32_t *it_iniP, 
    int32_t *it_finP, 
    int32_t *sP
  )
  {
    bool_t debug = TRUE;
    if (debug) { fprintf(stderr, "    %s: starting at buf.it_ini = %d\n", __FUNCTION__, it_start); }
    
    /* Data for a pulse: */
    int32_t it_ini = -1, it_fin = -1; /* First and last frame of pulse. */
    int32_t s = -1;  /* Phase type: stimulus type index {0..ns-1} or -1 if fixation. */ 
    
    int32_t it_next = it_start;
    neuromat_eeg_frame_buffer_find_next_pulse_start(buf, read_frame, it_next, ns, ichs, &it_ini, &s);
    if (it_ini != INT32_MIN) 
      { assert(it_ini >= it_next);
        assert((s >= 0) && (s < ns));
        int32_t ic = ichs[s]; /* Channel where pulse actually occurred. */
        it_next = it_ini;
        /* Skip rest of pulse: */
        neuromat_eeg_frame_buffer_find_pulse_end(buf, read_frame, it_next, ic, &it_fin);
        if (debug) { fprintf(stderr, "    %s: found pulse = %d..%d type %d\n", __FUNCTION__, it_ini, it_fin, s); }
        assert(it_fin >= it_ini);
      }
   else
     { if (debug) { fprintf(stderr, "    %s: no pulse found\n", __FUNCTION__); }
     }
    (*it_iniP) = it_ini;
    (*it_finP) = it_fin;
    (*sP) = s;
  }

void neuromat_eeg_frame_buffer_find_next_pulse_start
  ( neuromat_eeg_frame_buffer_t *buf, 
    neuromat_eeg_frame_buffer_read_proc_t *read_frame,
    int32_t it_start, 
    int32_t ns, 
    int32_t ichs[], 
    int32_t *itP,
    int32_t *sP
  )
  {
    int32_t it = it_start;
    int32_t s = -1; /* Channel that was up. */
    while (TRUE)
      { /* Make sure that frame {it} is in the buffer, grab it: */
        int32_t itb = neuromat_eeg_frame_buffer_get_frame(buf, it, read_frame, FALSE);
        if (itb < 0) 
          { /* End of file before next frame. */
            it = INT32_MIN; 
            break;
          }
        double *vt = buf->val[itb];
        /* Check if any of the channels are up, save in {*sP}: */
        int32_t nup = 0; /* Number of marker channels that are up. */ 
        for (int32_t r = 0; r < ns; r++)
          { /* Check channel {ichs[r]}: */
            int32_t ic = ichs[r];
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
    int32_t it_start, 
    int32_t ic, 
    int32_t *itP
  )
  {
    int32_t it = it_start;
    while(TRUE)
      { /* Make sure that frame {it} is in the buffer: */
        int32_t itb = neuromat_eeg_frame_buffer_get_frame(buf, it, read_frame, FALSE);
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

void neuromat_eeg_frame_buffer_get_next_pulse_pair
  ( neuromat_eeg_frame_buffer_t *buf,
    neuromat_eeg_frame_buffer_read_proc_t *read_frame,
    int32_t it_start, 
    int32_t ns, 
    int32_t ichs[], 
    int32_t nt_fx_default,
    bool_t verbose,
    char *chname[], 
    double fsmp,
    int32_t *it_fx_iniP,
    int32_t *it_fx_finP, 
    int32_t *it_st_iniP, 
    int32_t *it_st_finP,
    int32_t *s_stP
  )
  {
    int32_t s_fx = ns-1;  /* Phase type of the fixation phase: */

    /* Local frame indices: */
    int32_t it_fx_ini, it_fx_fin; /* First and last frame of start-of-fixation pulse. */
    int32_t it_st_ini, it_st_fin; /* First and last frame of start-of-stimulus pulse. */
    int32_t s_st;  /* Index of start-of-stimulus channel in {ichs[0..ns-1]}. */
    
    int32_t it_next = it_start;
    
    /* Loop until we get a well-formed run or the file is exhausted: */
    while(TRUE)
      { 
        /* Look for a fixation pulse, fake one if needed: */
        /* Note that it may be negative. */
        neuromat_eeg_frame_buffer_get_next_fixation_pulse
          ( buf, read_frame, it_next, 
            ns, ichs, nt_fx_default,
            verbose, chname, fsmp, 
            &it_fx_ini, &it_fx_fin
          );

        if (it_fx_ini == INT32_MIN)
          { /* No more runs: */
            it_st_ini = INT32_MIN;
            it_st_fin = INT32_MIN;
            s_st = INT32_MIN;
            break;
          }
        else
          { /* Found one more run, with real or faked start-of-fixation pulse. */
            /* Look for the start-of-stimulus pulse.  Note that it may be faked negative. */
            /* Allow overlapping pulses for now. */
            it_next = it_fx_ini + 1; 
            neuromat_eeg_frame_buffer_get_next_stimulus_pulse
              ( buf, read_frame, it_next, 
                ns, ichs, 
                verbose, chname, fsmp, 
                &it_st_ini, &it_st_fin, &s_st
              );
            if (it_st_ini == INT32_MIN)
              { /* Missing start-of-stimulus pulse, due to eof or another fixation pulse. */
                /* Discard fixation pulse: */
                fprintf(stderr, "!! start of stimulus pulse not found\n");
                fprintf(stderr, "!! ignoring start-of-fixation pulse at frame %d\n", it_fx_ini);
                it_next = it_fx_fin + 1;
                it_fx_ini = INT32_MIN;
                it_fx_fin = INT32_MIN;
                /* Look again for the start-of-fixation pulse: */
                continue;
              }
            else
              { /* Got a start-of-stimulus pulse. */
                break;
              }
          }
      }
      
    /* Either we found a start-of-stimulus pulse in the actual file, or eof: */
    if (it_st_ini != INT32_MIN)
      { /* Paranoia: */
        assert(it_fx_ini <= it_fx_fin);
        assert(it_fx_ini < it_st_ini);
        assert(it_st_ini <= it_st_fin);
        assert((s_st >= 0) && (s_st < ns));
        assert(s_st != s_fx);
       }
     else
       { assert(it_fx_fin == INT32_MIN);
         assert(it_st_ini == INT32_MIN);
         assert(it_st_fin == INT32_MIN);
         assert(s_st == INT32_MIN);
       }

    /* Return results: */
    (*it_fx_iniP) = it_fx_ini;
    (*it_fx_finP) = it_fx_fin;
    (*it_st_iniP) = it_st_ini;
    (*it_st_finP) = it_st_fin;
    (*s_stP) = s_st;
  }

void neuromat_eeg_frame_buffer_get_next_fixation_pulse
  ( neuromat_eeg_frame_buffer_t *buf,
    neuromat_eeg_frame_buffer_read_proc_t *read_frame,
    int32_t it_start,
    int32_t ns, 
    int32_t ichs[],
    int32_t nt_fx_default,
    bool_t verbose,
    char *chname[], 
    double fsmp,
    int32_t *it_fx_iniP,
    int32_t *it_fx_finP
  )
  {  
    int32_t s_fx = ns-1;  /* Phase type of the fixation phase: */

    /* Look for the next start-of-phase pulse in the fixation or stimulus trigger channels: */
    int32_t it_pu_ini, it_pu_fin; /* First and last frame of pulse. */
    int32_t s_pu;  /* Phase type in {0..ns-1}. */ 
    neuromat_eeg_frame_buffer_find_next_pulse
      ( buf, read_frame, it_start,
        ns, ichs, 
        &it_pu_ini, &it_pu_fin, &s_pu
      );

    if (it_pu_ini == INT32_MIN) 
      { /* No more start-of-something pulses in file: */
        (*it_fx_iniP) = INT32_MIN;
        (*it_fx_finP) = INT32_MIN;
        return;
      }
    else
      { /* Found a start-of-something pulse: */
        assert((it_pu_ini >= it_start) && (it_pu_ini <= it_pu_fin));
        assert((s_pu >= 0) && (s_pu < ns));

        if (s_pu == s_fx)
          { /* We got a start-of-fixation pulse: */
            int32_t ic_pu = ichs[s_pu];
            if (verbose)
              { neuromat_eeg_report_pulse
                  ( stderr, " [fx] ", ic_pu, chname[ic_pu], it_pu_ini, it_pu_fin, fsmp, "\n" );
              }
            (*it_fx_iniP) = it_pu_ini;
            (*it_fx_finP) = it_pu_fin;
            return;
          }
        else
          { /* We got an unexpected start-of-stimulus pulse. */
            int32_t ic_pu = ichs[s_pu];
            fprintf(stderr, "!! warning: missing start-of-fixation pulse - start-of-stimulus pulse in channel %d\n", ic_pu);
            /* Fake the start-of-fixation pulse: */
            if (it_pu_ini < it_start + nt_fx_default)
              { (*it_fx_iniP) = it_start; }
            else
              { (*it_fx_iniP) = it_pu_ini - nt_fx_default; }
            (*it_fx_finP) = (*it_fx_iniP);
            return;
          }
      }
  }
  
void neuromat_eeg_frame_buffer_get_next_stimulus_pulse
  ( neuromat_eeg_frame_buffer_t *buf,
    neuromat_eeg_frame_buffer_read_proc_t *read_frame,
    int32_t it_start, 
    int32_t ns, 
    int32_t ichs[], 
    bool_t verbose,
    char *chname[], 
    double fsmp,
    int32_t *it_st_iniP, 
    int32_t *it_st_finP,
    int32_t *s_stP
  )
  { 
    int32_t s_fx = ns-1;  /* Phase type of the fixation phase: */

    /* Look for the next start-of-phase pulse in the fixation or stimulus trigger channels: */
    int32_t it_pu_ini, it_pu_fin; /* First and last frame of pulse. */
    int32_t s_pu;  /* Phase type in {0..ns-1}. */ 
    neuromat_eeg_frame_buffer_find_next_pulse
      ( buf, read_frame, it_start, 
        ns, ichs, 
        &it_pu_ini, &it_pu_fin, &s_pu
      );

    if (it_pu_ini == INT32_MIN) 
      { /* No more start-of-something pulses in file: */
        (*it_st_iniP) = INT32_MIN;
        (*it_st_finP) = INT32_MIN;
        (*s_stP) = INT32_MIN;
        return;
      }
    else
      { /* Found a start-of-something pulse: */
        assert((it_pu_ini >= it_start) && (it_pu_ini <= it_pu_fin));
        assert((s_pu >= 0) && (s_pu < ns));

        if (s_pu != s_fx)
          { /* We got a start-of-stimulus pulse: */
            int32_t ic_pu = ichs[s_pu];
            if (verbose)
              { neuromat_eeg_report_pulse
                  ( stderr, " [st] ", ic_pu, chname[ic_pu], it_pu_ini, it_pu_fin, fsmp, "\n" );
              }
            (*it_st_iniP) = it_pu_ini;
            (*it_st_finP) = it_pu_fin;
            (*s_stP) = s_pu;
          }
        else
          { /* We got an unexpected start-of-fixation pulse: */
            int32_t ic_pu = ichs[s_pu];
            fprintf(stderr, "!! warning: missing start-of-stimulus pulse - start-of-fixation pulse in channel %d\n", ic_pu);
            /* Signal "start-of-stimulus pulse not found": */
            (*it_st_iniP) = INT32_MIN;
            (*it_st_finP) = INT32_MIN;
            (*s_stP) = INT32_MIN;
          }
        return;
      }
  }

