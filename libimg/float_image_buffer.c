/* See {float_image_buffer.h}. */
/* Last edited on 2024-12-05 00:44:12 by stolfi */

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <affirm.h>
#include <bool.h>
#include <rmxn.h>

#include <float_image_buffer.h>

/* Internal prototypes: */

double **float_image_buffer_alloc_rows(int32_t NC, int32_t NX, int32_t NB, bool_t clear);
  /* Allocates NB rows of {double} samples, each  with {NC*NX} samples. */

/* Implementations: */

float_image_buffer_t *float_image_buffer_new(int32_t NC, int32_t NX, int32_t NY, int32_t NB)
  {
    demand(NC >= 0, "invalid channel count");
    demand(NX >= 0, "invalid column count"); 
    demand(NY >= 0, "invalid row count");   
    demand(NB >= 0, "invalid buffer row count");
    /* Allocate top record: */
    float_image_buffer_t *buf = (float_image_buffer_t *)notnull(malloc(sizeof(float_image_buffer_t)), "no mem");
    buf->sz[0] = NC;
    buf->sz[1] = NX;
    buf->sz[2] = NY;
    /* Allocate row buffer vectors: */
    buf->NB = NB;
    buf->rowptr = float_image_buffer_alloc_rows(NC, NX, NB, FALSE);
    /* Set the buffer to empty band before row 0: */
    buf->yini = 0; buf->ylim = 0;
    return buf;
  }
  
void float_image_buffer_free(float_image_buffer_t *buf)
  { if (buf != NULL)
      { if (buf->rowptr != NULL)
          { int32_t k;
            for (k = 0; k < buf->NB; k++)
              { double **p = &(buf->rowptr[k]);
                if (*p != NULL) { free(*p); *p = NULL; }
              }
            free(buf->rowptr); buf->rowptr = NULL; 
          }
        free(buf); buf = NULL; 
      }
  }

double **float_image_buffer_alloc_rows(int32_t NC, int32_t NX, int32_t NB, bool_t clear)
  { double **rowptr = talloc(NB, double*);
    int32_t NCX = NC*NX; /* Samples per row. */
    int32_t yb;
    for (yb = 0; yb < NB; yb++) 
      { double *p = rmxn_alloc((uint32_t)NC, (uint32_t)NX);
        rowptr[yb] = p;
        if (clear) { int32_t k; for (k = 0; k < NCX; k++, p++) { (*p) = 0.0; } }
      }
    return rowptr;
  }

int32_t float_image_buffer_row_pos(float_image_buffer_t *buf, int32_t y)
  { 
    if (y < buf->yini) 
      { return -1; }
    else if (y >= buf->ylim) 
      { return +1; }
    else
      { return 00; }
  }

double *float_image_buffer_get_row(float_image_buffer_t *buf, int32_t y)
  { 
    if(float_image_buffer_row_pos(buf,y) == 00)
      { return buf->rowptr[y % (buf->NB)]; }
    else
      { return NULL; }
  }
  
void float_image_buffer_fill_row(float_image_buffer_t *buf, int32_t y, double val)
  { 
    double *row = float_image_buffer_get_row(buf, y);
    int32_t NCX = buf->sz[0]*buf->sz[1];
    int32_t k;
    for (k = 0; k < NCX; k++) { row[k] = val; }
  }

void float_image_buffer_advance(float_image_buffer_t *buf)
  { 
    int32_t y = buf->ylim;
    demand(y < buf->sz[2], "no next row");
    buf->ylim++;
    if (buf->ylim - buf->yini > buf->NB) { buf->yini++; }
    assert(buf->ylim - buf->yini <= buf->NB);
  }
