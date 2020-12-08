/* jsaudio.h - audio data representation, I/O, and manipulation */
/* Last edited on 2006-10-29 00:36:29 by stolfi */

/* 
  Derived from {rusound.h}, created by Rumiko Oishi Stolfi
  and Jorge Stolfi on sep/2006.
*/

#ifndef jsaudio_H
#define jsaudio_H

#include <stdio.h>
#include <stdint.h>

/* GENERIC SOUND CLIPS */

typedef struct sound_t /* A sound clip in memory. */
  { double fsmp;    /* Samples per second per channel. */
    int nc;         /* Number of channels. */
    int ns;         /* Number of samples per channel. */
    double **sv;    /* {sv[c][k]} is value of sample {k} in channel {c}. */
  } sound_t;
  
sound_t jsa_allocate_sound(int nc, int ns);
  /* Allocates memory for a {sound_t} with {nc} channels and {ns}
    samples per channel.  The samples are all set to 0.0. */

void jsa_free_sound(sound_t *s);
  /* Frees the memory used by the sound clip {s} (except the struct
    {s} itself). */

sound_t jsa_copy_sound(sound_t *s, int ini, int ns);
  /* Creates a copy of a portion of sound clip {s}. 
    The portion starts at sample number {ini} and has {ns} samples,. */

void jsa_add_sound(sound_t *s, int sskip, sound_t *r, int rskip, int ns);
  /* Adds {ns} samples of sound {r}, starting at sample number {rskip},
    to sound {s}, starting at sample number {sskip}.  The operation
    is limited to samples and channels that exist in both clips. */

#endif
