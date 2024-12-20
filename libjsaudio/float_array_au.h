/* {float_array_au.h} - read/write Sun AU audio files as {float_array_t}s */
/* Last edited on 2024-12-05 10:32:03 by stolfi */

/* Derived from {rusound.h}, created by Rumiko Oishi Stolfi
  and Jorge Stolfi on sep/2006.
*/

#ifndef float_array_au_H
#define float_array_au_H

#include <stdint.h>
#include <stdio.h>

#include <float_array.h>
#include <jsaudio.h>

/* CONVERSION BETWEEN AU SAMPLE ARRAYS AND FLOAT ARRAYS

  The procedures in this interface perform conversion between
  a {float_array_t} {A} (with float samples) and a {sound_t} 
  structure {snd} (with double-precision samples).  

  The indices of the {float_array_t} will have the following meaning:

    | {ix[0]} = sound channel.
    | {ix[1]} = sampling time.

  The remaining indices {ix[2..5]} are trivial (always 0). */

float_array_t *float_array_from_sound(sound_t *snd, bool_t verbose);
  /* Converts a sound array {snd} to a {float_array_t}.
    Samples are just converted from {double} to {float}.

    The returned array has two indices, /channel/ and /frame/.
    
    If {snd} is mono, the channel index {ix[0]} has only one
    legal value, namely 0. If {snd} is stereo, then {ix[0] = 0,1} 
    select the left or right channel, respectively.  The 
    frame index {ix[1]} selects the sampling instant.

    If {verbose} is TRUE, the procedure prints statistics of the
    conversion to {stderr}. */

sound_t *float_array_to_sound(float_array_t *A, int32_t chns, int32_t ch[], double freq, bool_t verbose);
  /* Converts a {float_array_t} {A} to a sound file {snd}.
  
    The number of channels of {snd} will be {chns} --- usually 1 for
    mono, 2 for stereo, etc.. The samples of channel {k} of {snd} will
    be taken from the subarray of {A} with {ix[0] = ch[k]}, for each
    {k} in {0..chns-1}. If {ch[k]} is invalid (negative, or {>=
    A.ds.sz[0]}), the samples in channel {k} of {snd} are all set to
    zero. If {ch} is NULL, it defaults to the identity vector
    {(0,1,2,...chns-1)}.
    
    The sampling frame index is {ix[1]}, so the number of samples per channel
    is {A.ds.sz[1]}. 
    
    The parameter {freq} is the sampling frequency in Hz.  It is stored as 
    {snd.fsmp} but is not used for anything else.
    
    
    
    

    If {verbose} is TRUE, the procedure prints statistics of the
    conversion to {stderr}. */
  
#endif

