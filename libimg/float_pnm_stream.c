/* See {float_pnm_stream.h}. */
/* Last edited on 2017-06-22 02:39:44 by stolfilocal */

#include <assert.h>
#include <affirm.h>
#include <bool.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <float_image_buffer.h>

#include <float_pnm_stream.h>

float_pnm_stream_t *float_pnm_stream_new(bool_t isMask, uint32_t badval)
  { /* Allocate top record: */
    float_pnm_stream_t *str = (float_pnm_stream_t *)notnull(malloc(sizeof(float_pnm_stream_t)), "no mem");
    /* Fileds that can be set now: */
    str->isMask = isMask;
    str->badval = badval; 
    /* Fields that are set in different ways depending on stream's direction: */
    str->chns = 0;
    str->rows = 0;
    str->cols = 0;
    str->maxval = 0;
    str->format = ('?' << 8) | '?';
    str->raw = FALSE;
    str->bits = FALSE;
    str->ftb = NULL;
    str->smp = NULL;
    str->buf = NULL;
    return str;
  }


void float_pnm_stream_free(float_pnm_stream_t *str)
  { if (str != NULL)
      { if (str->smp != NULL) { free(str->smp); str->smp = NULL; }
        if (str->ftb != NULL) { free(str->ftb); str->ftb = NULL; }
        if (str->buf != NULL) { float_image_buffer_free(str->buf); }
        free(str); str = NULL; 
      }
  }
