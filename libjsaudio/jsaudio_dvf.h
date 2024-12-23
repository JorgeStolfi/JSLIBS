/* jsaudio_dvf.h - Sony Recorder DVF audio file I/O */
/* Last edited on 2024-12-21 03:31:40 by stolfi */

/* 
  Derived from {rusound.h}, created by Rumiko Oishi Stolfi
  and Jorge Stolfi on sep/2006.
*/

#ifndef jsaudio_dvf_H
#define jsaudio_dvf_H

#include <stdio.h>
#include <stdint.h>

#include <jsaudio.h>

/* LOW-LEVEL AUDIO FILE I/O */

/* 
  THE SONY RECORDER ".dvf" AUDIO FILE FORMAT
  
  The .dvf file format used by Sony recorders is proprietary and
  undocumented. Therefore, this interface is limited to extracting a
  few basic header parameters, such as the date and time of the
  recording.
  
  Here is the deduced header format.  All multi-byte numbers are big-endian.
  
    Bytes      Contents
    ---------  -----------------------------
    0..7:      "MS_VOICE"
    8..11:     ??? = 80 
    12         ??? = 1
    13         ??? = 2
    14..15     zero?
    16..31     "SONY CORPORATION"
    32         80 'P'
    33..41     "Voice.dll"
    42..47     zero?
    48         1
    49..51     zero?
    52..58     Date I (start of recording?):
                 52..53  Year (4 digit) 2006 (7,214)
                 54      Month (1..12)
                 55      Day (1..31)
                 56      Hour (0..11 or 1..12?)
                 57      Minute (0..59)
                 58      Second (0..59)
    59         ??? 1 | 3 | 0 | 2 | 2 | 2 | 0
    60         zero?
    61         53
    62         0
    63         1
    64         0
    65         0
    66         28
    67         32
    68         0
    69         0
    70         3
    71         132
    72         4
    73         0
    74         0
    75         1
    76..79     zero?
    80         1
    81..88     20482? (80,2)
    89..96     24581? (96,5)
    97..104    65542? (1,0,6)
    105..112   4103? (16,7)
    113..120   8201? (32,9)
    121..128   4106? (16,10)
    129..136   4107? (16,11)
    137..144   4109? (16,13)
    145..152   4110? (16,14)
    153..160   262147? (4,0,3)
    161..168   122884? (1,224,4)
    169        0
    170        ???  1 | 0 |   1 |  3 |   3 |  3 |  1
    171        ??? 80 | 6 | 224 | 66 | 120 | 18 | 20
    172        zero?
    173        ???  68 |  56 | 12 | 161 |  39 |  0 | 79
    174        ??? 192 | 116 | 80 | 140 | 160 | 44 | 24
    175        zero?
    176        144
    177..431   zero?
    432        1
    433        0
    434..447   "SONY CORP.    "
    448        144
    449..479   zero?
    480..483   Recording time in seconds times 900. A multiple of 216.
    484..495   zero?
    496..499   255 '\377'
    500..544   zeros?
    545        ??? 158 | 59 | 119 | 174 | 255 | 150 | 87
    546        0
    547        255
    548..554   Date II:
               548..549   Year (4 digit) 2006 (7,214)
               550        Month (1..12)     
               551        Day (1..31)        
               552        Hour (0..23 or 1..24?)
               553        Minute (0..59)
               554        Second (0..59)
    555        ??? 1 | 3 | 0 | 2 | 2 | 2 | 0
    556..565   "SONY ICD-P"
    566..579   zero?
    580        34
    581        5
    582        7
    583        27
    584        204
    585        166
    586        172
    587        197
    588..595   zero?
    596..599   255 '\377'
    
  The two time/date fields are nearly identical, but field II is 
  often 1 second later.
  
  The field in bytes 480..483 is precisely the number of seconds of
  the OGG file, times 900. Let's therefore define a "Sony time unit"
  as 1/900 sec. Moreover, the recording time appears to be always a
  multiple of 8*27 = 216 Sony time units, that is, 0.24 seconds.
  Assuming 8000 samples per second, that means 1920 samples.
  
  The field in bytes 173..174 is roughly the number of samples divided
  by 2253.5, or the number of Sony time units divided by 253.5, or
  3.55 times the number of seconds. It may be combined with the field
  in bytes 170..171, in some fashion.

*/

typedef struct jsaudio_dvf_time_t
  { uint16_t year;      /*  52  53 Year (4 digit). */
    uint8_t  month;     /*  54     Month (1..12). */
    uint8_t  day;       /*  55     Day (1..31). */
    uint8_t  hour;      /*  56     Hour (0..23 or 1..24?). */
    uint8_t  minute;    /*  57     Minute (0..59). */
    uint8_t  second;    /*  58     Second (0..59). */
  } jsaudio_dvf_time_t;

/* Variable fields from the Sony Recorder DVF file header: */
typedef struct jsaudio_dvf_file_header_t 
  { 
    jsaudio_dvf_time_t time_1;    /*  52  58 Recording date and time (1). */
    uint8_t  unk_059_1;   /*  59     Unknown (0..3?). */
    uint8_t  unk_170_1;   /* 170     Unknown (0..3?). */
    uint8_t  unk_171_1;   /* 171     Unknown. */
    uint16_t unk_173_2;   /* 173 174 Unknown (related to length). */
    uint32_t length_su;   /* 480 483 Recording time in seconds times 900; a multiple of 216. */
    uint8_t  unk_545_1;   /* 545     Unknown. */           
    jsaudio_dvf_time_t time_2;    /* 548 554 Recording date and time (2). */
    uint8_t  unk_555_1;   /* 555     Unknown. */           
    char model[24];       /* 556 579 Recorder model, e.g. "SONY ICD-P", zero-padded. */
  } jsaudio_dvf_file_header_t;

jsaudio_dvf_file_header_t jsaudio_dvf_read_file_header(FILE *rd);
  /* Reads the header of a Sony ".dvf" audio file from stream {rd}. 
    Leaves the stream positioned just before the first sample data byte. */
 
#endif

