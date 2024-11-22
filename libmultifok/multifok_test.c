/* See {multifok_test.h}. */
/* Last edited on 2024-10-10 22:37:27 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <wt_table.h>
#include <wt_table_hann.h>
#include <affirm.h>
#include <interval.h>
#include <fget.h>
#include <i2.h>
#include <r2.h>
#include <r3.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <jsqroots.h>
#include <rn.h>

#include <float_image.h>
#include <float_image_paint.h>
#include <float_image_read_pnm.h>
#include <float_image_write_pnm.h>

#include <multifok_window.h>
#include <multifok_window.h>
#include <multifok_basis.h>
#include <multifok_score.h>
#include <multifok_term.h>
#include <multifok_scene.h>

#include <multifok_test.h>
    
/* IMPLEMENTATIONS */

FILE* multifok_test_open_text_file(char *outPrefix, char *tag)
  { 
    char *fname = jsprintf("%s%s.txt", outPrefix, tag);
    FILE *wr = open_write(fname, TRUE);
    free(fname);
    return wr;
  }

#define multifok_test_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

