/* jsaudio_au.h - Sun AU audio file I/O */
/* Last edited on 2024-12-21 03:46:05 by stolfi */

/* 
  Derived from {rusound.h}, created by Rumiko Oishi Stolfi
  and Jorge Stolfi on sep/2006.
*/

#ifndef jsaudio_au_H
#define jsaudio_au_H

#include <stdio.h>
#include <stdint.h>

#include <jsaudio.h>

/* HIGH-LEVEL AUDIO FILE I/O */

jsaudio_t jsaudio_au_read_file(FILE *rd);
  /*  Reads from file {rd} a sound file in Sun's ".au" format. 
    Returns the samples as a {jsaudio_t} structure. */
     
void jsaudio_au_write_file(FILE *wr, jsaudio_t *s);
  /* Writes to file {wr} the sound clip {s},
    in Sun's ".au" format (16 bits signed complement, 
    linear encoding). */

/* LOW-LEVEL AUDIO FILE I/O */

/* 
  THE SUN ".au" AUDIO FILE FORMAT
  
  The .au file format, originally by SUN, is (fortunately) a very
  straightforward audio format, unfortunately it isn't widely
  supported outside the UNIX community. The file format is split into
  three parts....

   * A header containing basic information such as the length, number
     of channels, sample frequency, and data format.

   * A variable length informational field. This is designed for
     copyright information, authors name, etc.

   * The audio data which may be stored in a number of formats.

Header

  The header is a total of six 4 byte quantities. The description of
  each 4 byte word is described below.

    * A "magic" identification word of ".snd", otherwise known as 0x2e736e64.

    * A data offset (bytes) to the audio data. If there is no information 
      section then this would be 24, the size of the remainder of this header.

    * The data size (bytes) of the audio data. Since this can be worked 
      out knowing the file length and the data offset above, it is
      permitted to set this to 0xffffffff.

    * The encoding type used for the data, a number 
      from 1, 2, 3, 4, 5, 6, 7, 23, 24, 25, 26, 27. 
      See later.

    * The sampling frequency in Hz. The most commonly used 
      frequencies are 11025, 16000, 22050, 32000, 44100, and 48000 Hz.

    * The number of channels. For more than one channel the 
      samples are interleaved.

Information

  This can be any information at all, it automatically starts 24 bytes
  from the start of the file. If the data offset is 0 then this
  section is empty. The length of this section is calculated by
  subtracting 24 from the data offset field in the header.
  
Audio samples

  The audio samples can be encoded in a number of formats, the exact
  format is described in the 4th word of the header.  The encodings
  are:
     
        1 = 8 bit ISDN u-law
        2 = 8 bit linear PCM
        3 = 16 bit linear PCM
        4 = 24 bit linear PCM
        5 = 32 bit linear PCM
        6 = 32 bit IEEE floating point
        7 = 64 bit IEEE floating point
       23 = 4 bit CCITT G721 ADPCM
       24 = CCITT G722 ADPCM
       25 = CCITT G723 ADPCM
       26 = 5 bit CCITT G723 ADPCM
       27 = 8 bit ISDN A-law

*/

/* Fields from the Sun AU file header: */
typedef struct jsaudio_au_file_header_t 
  { uint32_t magic;       /* Magic number. */
    uint32_t hdr_size;    /* Size of this header. */
    uint32_t data_size;   /* Length of data (optional). */
    uint32_t encoding;    /* Data encoding format. */
    uint32_t sample_rate; /* Samples per second. */
    uint32_t channels;    /* Number of interleaved channels. */
  } jsaudio_au_file_header_t;

uint32_t jsaudio_au_file_bytes_per_sample(uint32_t enc);
  /* Returns the number of bytes per sample for the encoding type {enc}.
     Works only for encodings where each sample has a fixed number 
     of bytes; namely {enc = 1,2,3,4,5,6,7,27}.  Bombs out for 
     other encodings. */

jsaudio_au_file_header_t jsaudio_au_read_file_header(FILE *rd);
  /* Reads the header of a Sun ".au" audio file from stream {rd}. 
    Leaves the stream positioned just before the first sample data byte. */

void jsaudio_au_skip_file_samples(FILE *rd, jsaudio_au_file_header_t *h, uint32_t ns);
  /* Skips {ns} samples from each channel of a Sun ".au" file opened
    as stream {rd}.  Assumes that the file header data is {h} 
    (so that the total samples skipped is {ns*h.nc}) and that {rd} is
    positioned just before the first byte of a sample of 
    channel 0. Fails if EOF is encountered prematurely. */

void jsaudio_au_read_file_samples(FILE *rd, jsaudio_au_file_header_t *h, jsaudio_t *s, uint32_t skip, uint32_t ns);
  /* Reads {ns} samples from each channel of a Sun ".au" file opened as
    stream {rd}. Assumes that the file header data is {h} (so that the
    total samples skipped is {ns*h.nc}) and that {rd} is positioned just
    before the first byte of a sample of channel 0. The samples are
    stored into the sound clip {s}, starting at sample number {skip} of
    each channel. Fails if {skip..skip+ns-1} is not a subset of
    {0..s.ns-1}, or if EOF is encountered prematurely. */

void jsaudio_au_write_file_header(FILE *wr, jsaudio_au_file_header_t *h);
  /* Writes a Sun ".au" file header, with data taken from {h}, to stream {wr}. */
  
void jsaudio_au_write_file_samples(FILE *wr, jsaudio_au_file_header_t *h, jsaudio_t *s, uint32_t skip, uint32_t ns);
  /* Writes samples {s[skip..skip+ns-1]} to stream {wr}, assumed to contain
    a Sun ".au" file with header data {h}. */
  
#endif

