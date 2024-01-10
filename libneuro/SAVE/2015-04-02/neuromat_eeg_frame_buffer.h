#ifndef neuromat_eeg_frame_buffer_H
#define neuromat_eeg_frame_buffer_H

/* Buffer for readin consecutive dataframes from plain text NeuroMat EEG datasets. */
/* Last edited on 2013-11-13 00:08:40 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>

#include <bool.h>
 
typedef struct neuromat_eeg_frame_buffer_t 
  {
    int size;         /* Max number of frames in buffer. */
    int nc;       /* Total number of channels. */
    double **val; /* Cyclic queue of frame pointers. */
    int it_ini;   /* Index of first frame in buffer. */
    int it_fin;   /* Index of last frame in buffer. */
  } neuromat_eeg_frame_buffer_t;
  /* A buffer that stores up to {size} consecutive 
    data frames from a text NeuroMat EEG dataset, each with {nc}
    data channels (including electrodes, triggers, reference, and other event marker channels).
    A frame with time index {it} in the range {it_ini..it_fin} is stored 
    into {val[it % size][0..nc]}; other frames are not stored.
    Entries of {val} that are not used in this way are either {NULL} or 
    point to an unused data frame record.
    
    The index {it_ini} is always non-negative, and {it_fin} is {it_ini-1} iff
    the buffer is empty, at least {it_ini} if it has at least one frame loaded. */

neuromat_eeg_frame_buffer_t *neuromat_eeg_frame_buffer_new(int size, int nc);
  /* Allocate EEG data buffer, with space for up to {size} frames of {nc}
    channels each. */

void neuromat_eeg_frame_buffer_free(neuromat_eeg_frame_buffer_t *buf);
  /* De-allocates all the internal subsidiary records of {buf}, and the 
    record {*buf} itself. */

typedef int neuromat_eeg_frame_buffer_read_proc_t(int nc, double val[]);
  /* Type of a client procedure that reads a frame from somewhere and stores it into {val[0..nc-1]}.
    Should return 1 if success, 0 if there are no more frames. */

int neuromat_eeg_frame_buffer_get_frame
  ( neuromat_eeg_frame_buffer_t *buf,
    int it,
    neuromat_eeg_frame_buffer_read_proc_t read_frame
  );
  /* Tries to ensure that the requested data frame {it} is in {buf}.
  
    Namely, if {it} is already in {buf.it_ini..buf.it_fin}, does nothing.
    If {it} is greater than {buf.it_fin}, reads one or more frames from {rd}
    using the client-given {read_frame} procedure, stores them into
    their places in the buffer, and increments {buf.it_fin};
    until {buf.it_fin} is equal to {it}. Allocates frame storage as needed. If necessary,
    reuses one or more buffer frames, starting at index {buf.it_ini} and
    incrementing the latter. 
    
    If it succeeds, returns a non-negative index {itb = (it%buf.size)}
    such that {buf.val[itb]} points to the frame with index {it} in the
    whole dataset. Returns {-1} if the frames were exhausted before the frame
    with index {it}. */

void neuromat_eeg_frame_buffer_find_next_pulse
  ( neuromat_eeg_frame_buffer_t *buf, 
    neuromat_eeg_frame_buffer_read_proc_t *read_frame,
    int it_start, 
    int ns,
    int ichs[], 
    int *it_iniP, 
    int *it_finP, 
    int *sP
  );
  /* Looks for the first and last frame of the next pulse in certain
    channels of the input, starting with frame {it_start}. Assumes that
    some frames are in {buf}, reads more frames with {read_frame} as
    needed. The frame {it_start} must not have been flushed from buffer
    {buf} yet.

    Specifically, the procedure assumes that the pulses of interest
    pulses are in channels {ichs[0..ns-1]}. It looks for the first frame
    at or after frame {it_start} where one of those channels is
    positive,and then for the last frame where that same channel is
    positive, extending possibly to the end of the dataset.

    If it succeeds, returns in {*it_iniP,*it_finP} the first and
    last frames of the pulse. Also returns in {*sP} the the index
    in {0..ns-1} such that the the pulse was in channel {ichs[*sP]}.

    If all channels {chs[0..ns-1]} are zero in all frames at or after frame
    {it_start}, the procedure returns {-1} in {*it_fx_iniP}. In that case
    {*it_fx_finP,*sP} are undefined. */

void neuromat_eeg_frame_buffer_find_next_pulse_start
  ( neuromat_eeg_frame_buffer_t *buf, 
    neuromat_eeg_frame_buffer_read_proc_t *read_frame,
    int it_start, 
    int ns, 
    int ichs[], 
    int *itP,
    int *sP
  );
  /* Finds the next up-transition in one or more channels of the input,
    starting at or after frame {it_start}.  Assumes that some frames are in {buf}, reads more frames
    with {read_frame} as needed. The frame {it_start} must not have been flushed from buffer {buf} yet.
    
    Specifically, the procedure assumes that the pulses of interest are
    in channels {ichs[0..ns-1]}. It looks for the first frame, beginning
    with {it_start}, where one of those channels is positive. 
    
    If successful, returns in {*itP} the index of that frame, and in {*sP} the
    channel index in {0..ns-1} among the specified -- that is, an integer
    such that the pulse occurred in channel {ichs[*sP]}.  If there are
    two or more pulses starting at the same time, returns the one with smallest
    {*sP}.
    
    Returns {-1} in {*itP} if the dataset ends before finding
    the required pulse.  In that case, {*sP} will be undefined.  */

void neuromat_eeg_frame_buffer_find_pulse_end
  ( neuromat_eeg_frame_buffer_t *buf, 
    neuromat_eeg_frame_buffer_read_proc_t *read_frame, 
    int it_start, 
    int ic, 
    int *itP
  );
  /* Finds the frame at or after frame {it_start} where the current
    pulse in channel {ic} ends. Returns its index in {*itP}. Assumes
    that some frames are in {buf}, reads more frames with {read_frame}
    as needed. The frame {it_start} must not have been flushed from
    buffer {buf} yet.
    
    Specifically, the procedure requires that channel {ic} of frame
    {t_start} is positive. The procedure then looks for the first frame
    after that where that channel is zero. 
    
    If successful, returns in {*itP} the index pf the *previous* frame
    (the last one where the channel was nonzero).
    
    If the dataset ends before finding that frame, it returns in
    {*itP} the index of the last existing frame. */

#endif
