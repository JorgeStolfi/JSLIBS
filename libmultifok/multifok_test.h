/* Test tools for {multifok_focus_op} and related funcs. */
/* Last edited on 2024-10-10 22:33:13 by stolfi */

#ifndef multifok_test_H
#define multifok_test_H

#define _GNU_SOURCE
#include <stdint.h>

#include <interval.h>
#include <i2.h>
#include <r2.h>
#include <r3.h>
#include <bool.h>
#include <frgb.h>
#include <float_image.h>

#include <multifok_scene.h>
#include <multifok_term.h>

/* BASIS AND TERM DATA I/O */

FILE* multifok_test_open_text_file(char *outPrefix, char *tag);
  /* Opens a text file called "{outPrefix}{tag}.txt" for writing,
    returns the handle. */


#endif
