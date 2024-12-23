/* jsaudio.h - audio data representation, I/O, and manipulation */
/* Last edited on 2024-12-21 03:21:23 by stolfi */

/* 
  Derived from {rusound.h}, created by Rumiko Oishi Stolfi
  and Jorge Stolfi on sep/2006.
*/

#ifndef jsaudio_H
#define jsaudio_H

#include <stdio.h>
#include <stdint.h>

/* GENERIC SOUND CLIPS */

typedef struct jsaudio_t /* A sound clip in memory. */
  { double fsmp;    /* Samples per second per channel. */
    uint32_t nc;    /* Number of channels. */
    uint32_t ns;    /* Number of samples per channel. */
    double **sv;    /* {sv[ic][k]} is value of sample {k} in channel {ic}. */
  } jsaudio_t;
  
jsaudio_t jsaudio_allocate_sound(uint32_t nc, uint32_t ns);
  /* Allocates memory for a {jsaudio_t} with {nc} channels and {ns}
    samples per channel.  The samples are all set to 0.0. */

void jsaudio_free_sound(jsaudio_t *s);
  /* Frees the memory used by the sound clip {s} (except the struct
    {s} itself). */

jsaudio_t jsaudio_copy_sound(jsaudio_t *s, uint32_t ini, uint32_t ns);
  /* Creates a copy of a portion of sound clip {s}. 
    The portion starts at sample number {ini} and has {ns} samples,. */

void jsaudio_add_sound(jsaudio_t *s, uint32_t sskip, jsaudio_t *r, uint32_t rskip, uint32_t ns);
  /* Adds {ns} samples of sound {r}, starting at sample number {rskip},
    to sound {s}, starting at sample number {sskip}.  The operation
    is limited to samples and channels that exist in both clips. */

#endif
