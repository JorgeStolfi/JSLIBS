#ifndef neuromat_eeg_frame_buffer_H
#define neuromat_eeg_frame_buffer_H

/* Buffer for readin consecutive dataframes from plain text NeuroMat EEG datasets. */
/* Last edited on 2017-10-18 13:15:40 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <limits.h>

#include <bool.h>
 
typedef struct neuromat_eeg_frame_buffer_t 
  {
    int size;     /* Max number of frames in buffer. */
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
    neuromat_eeg_frame_buffer_read_proc_t read_frame,
    bool_t mirror
  );
  /* Tries to ensure that the requested data frame {it} is in {buf}.
  
    Namely, if {it} is already in {buf.it_ini..buf.it_fin}, does nothing.
    If {it} is greater than {buf.it_fin}, reads one or more frames from {rd}
    using the client-given {read_frame} procedure, stores them into
    their places in the buffer, and increments {buf.it_fin};
    until {buf.it_fin} is equal to {it}. Allocates frame storage as needed. If necessary,
    reuses ("flushes") one or more buffer frame slots, starting at index {buf.it_ini} and
    incrementing the latter. 
    
    If it succeeds, returns a non-negative index {itb = (it%buf.size)}
    such that {buf.val[itb]} points to the frame with index {it} in the
    whole dataset.
    
    Fails, returning {-1}, if the file contains
    no frames, or if {mirror} is false and
    {it} is outside the range of frame indices in the file,
    namely {0..nt-1} for a file with {nt} frames.
    
    If {mirror} is true and the file is non-empty, any index {it}
    outside the range {0..nt-1} is mapped back into that range by
    folding about the extremes. Namely, if {it < 0} is negative,
    replaces it internally by its absolute value, thus implicitly
    mirroring the file about frame 0. Likewise, if {it >= nt}, replaces
    {it} internally by {2*(nt-1) - it}, thus implicitly mirroring the
    file about frame {nt-1}. These mirrorings are repeated as often as
    neccessary to get {it} in the range {0..nt-1}. */

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
    are in channels {ichs[0..ns-1]}. It looks for the first frame
    at or after frame {it_start} where one of those channels is
    positive, and then for the last frame where that same channel is
    positive, extending possibly to the end of the dataset.

    If it succeeds, returns in {*it_iniP,*it_finP} the first and
    last frames of the pulse. Also returns in {*sP} the the index
    in {0..ns-1} such that the the pulse was in channel {ichs[*sP]}.

    If all channels {chs[0..ns-1]} are zero in all frames at or after frame
    {it_start}, the procedure returns {INT_MIN} in {*it_fx_iniP}. In that case
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
    
    Returns {INT_MIN} in {*itP} if the dataset ends before finding
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

/* TWO-PHASE EXPERIMENT RUNS 
  
  The following procedures are intended to locate
  experimental runs consisting of a "fixation" phase
  followed by a "stimulus" phase.  They assume that the start
  of each phase is indicated by the rising edge of a pulse
  in some marker channel.
  
  For now, the procedures assume that the start-of-fixation pulse
  occcurs in only one marker channel, while the start-of-stimulus
  pulse may occur in one of several channels, depending on the 
  run type. */

void neuromat_eeg_frame_buffer_get_next_pulse_pair
  ( neuromat_eeg_frame_buffer_t *buf,
    neuromat_eeg_frame_buffer_read_proc_t *read_frame,
    int it_start, 
    int ns, 
    int ichs[], 
    int nt_fx_default,
    bool_t verbose,
    char *chnames[], 
    double fsmp,
    int *it_fx_iniP,
    int *it_fx_finP, 
    int *it_st_iniP, 
    int *it_st_finP,
    int *s_stP
  );
  /* Looks for the next experimental run in the input file, defined by a
    start-of-fixation pulse (possibly missing) followed by a start-of-stimulus pulse,
    starting at or after frame {it_start}. Assumes that some frames are
    in {buf}; reads more with {read_frame} if needed. The frame
    {it_start} must not have been flushed from {buf} yet.

    Specifically, the procedure assumes that the valid start-of-stimulus pulses are
    in channels {ichs[0..ns-2]}, and that the start-of-fixation
    pulses are in channel {ic_fx = ichs[ns-1]}. It looks for the first frame,
    beginning with {it_start}, where channel {ic_fx} is positive, and
    considers that the first frame of the start-of-fixation pulse. Then
    looks for the first frame after that where some channel {ic_st = ichs[s_st]}
    becomes positive, for some {s_st} in {0..ns-2}; and considers
    that the first frame of the start-of-stimulus pulse.

    If it succeeds, returns in {*it_fx_iniP,*it_fx_finP} the first and
    last frames of the start-of-fixation pulse, and in
    {*it_st_iniP,*it_st_finP} the first and last frames of the
    start-of-stimulus pulse. Also returns in {*s_stP} the the stimulus
    type in {0..ns-2}. Reports the pulse to
    {stderr}, assuming {chnames[o..buf->nc-1]} are the channel names and
    {fsmp} is the sampling rate (in Hz).

    If there are no more runs in the input, returns {INT_MIN} in {*it_fx_iniP}.
    In that case {*it_fx_finP,*it_st_iniP,*it_st_finP,*s_stP} are undefined.
    
    If the first pulse found after {it_start} is a start-of-stimulus pulse, tries 
    to get around the situation by assuming the the lost start-of-fixation 
    pulse came {nt_fx_default} frames before the start-of-stimulus pulse
    (or at the first frame, if there is not enough frames for that).
    Note that, as a result, {*it_fx_ini_P} and possibly {*it_fx_finP}
    may be negative.
    
    If two or more start-of-fixation pulses are found with no intervening
    start-of-stimulus pulse, ignores all but the last one. */
    
void neuromat_eeg_frame_buffer_get_next_fixation_pulse
  ( neuromat_eeg_frame_buffer_t *buf,
    neuromat_eeg_frame_buffer_read_proc_t *read_frame,
    int it_start, 
    int ns, 
    int ichs[],
    int nt_fx_default,
    bool_t verbose,
    char *chnames[], 
    double fsmp,
    int *it_fx_iniP,
    int *it_fx_finP
  );
  /* Looks for the next start-of-fixation pulse in the input file,
    faking one if it is missing, starting at or after frame {it_start}.
    Assumes that some frames are in {buf}; reads more with {read_frame}
    if needed. The frame {it_start} must not have been flushed from
    {buf} yet.
    
    Specifically, the procedure assumes that the valid start-of-stimulus
    pulses are in channels {ichs[0..ns-2]}, and that the
    start-of-fixation pulses are in channel {ic_fx = ichs[ns-1]}. It
    looks for the first frame, beginning with {it_start}, where channel
    {ic_fx} is positive, and considers that the first frame of the
    start-of-fixation pulse.

    If it succeeds, returns in {*it_fx_iniP,*it_fx_finP} the first and
    last frames of the start-of-fixation pulse. Reports the pulse to
    {stderr}, assuming {chnames[o..buf->nc-1]} are the channel names and
    {fsmp} is the sampling rate (in Hz).

    If there are no more runs in the input, returns {INT_MIN} in {*it_fx_iniP}.
    In that case {*it_fx_finP} is undefined.
    
    If the first pulse found after {it_start} is a start-of-stimulus pulse, tries 
    to get around the situation by assuming the the start-of-fixation 
    pulse came {nt_fx_default} frames before the start-of-stimulus pulse.
    Note that this may result in a negative value of {*it_fx_iniP}
    and maybe also {*it_fx_finP}. */
  
void neuromat_eeg_frame_buffer_get_next_stimulus_pulse
  ( neuromat_eeg_frame_buffer_t *buf,
    neuromat_eeg_frame_buffer_read_proc_t *read_frame,
    int it_start, 
    int ns, 
    int ichs[], 
    bool_t verbose,
    char *chnames[], 
    double fsmp,
    int *it_st_iniP, 
    int *it_st_finP,
    int *s_stP
  );
  /* Looks for the next start-of-stimulus pulse in the input file,
    faking one if it is missing, starting at or after frame {it_start}.
    Assumes that some frames are in {buf}; reads more with {read_frame}
    if needed. The frame {it_start} must not have been flushed from
    {buf} yet.
    
    Specifically, the procedure assumes that the valid start-of-stimulus
    pulses are in channels {ichs[0..ns-2]}, and that the
    start-of-fixation pulses are in channel {ic_fx = ichs[ns-1]}. It
    looks for the first frame, beginning with {it_start}, where any channel
    {ichs[s_st]} is positive, for any {s_st} in {0..ns-2},
    and considers that the first frame of the
    start-of-stimulus pulse.

    If it succeeds, returns in {*it_st_iniP,*it_st_finP} the first and
    last frames of the start-of-stimulus pulse, and in {*s_stP} the 
    phase {s_st}. Reports the pulse to
    {stderr}, assuming {chnames[o..buf->nc-1]} are the channel names and
    {fsmp} is the sampling rate (in Hz).

    If reaches the end-of-file before finding such a frame,
    returns {INT_MIN} in {*it_st_iniP}. In that case {*it_st_finP} 
    and {*s_stP} are undefined.
    
    If the first pulse in channels {ichs[0..ns-1]} 
    found after {it_start} is a start-of-fixation pulse, 
    rather than a start-of-stimulus, also returns  {INT_MIN} in {*it_st_iniP}. */

#endif
