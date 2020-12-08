/* See {float_image_buffer.h}. */
/* Last edited on 2010-08-04 17:29:22 by stolfi */

#include <stdlib.h>
#include <assert.h>
#include <affirm.h>
#include <bool.h>

#include <float_image_buffer.h>

/* Internal prototypes: */

double **float_image_buffer_alloc_rows(int NC, int NX, int NB, bool_t clear);
  /* Allocates NB rows of {double} samples, each  with {NC*NX} samples. */

/* Implementations: */

float_image_buffer_t *float_image_buffer_new(int NC, int NX, int NY, int NB)
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
          { int k;
            for (k = 0; k < buf->NB; k++)
              { double **p = &(buf->rowptr[k]);
                if (*p != NULL) { free(*p); *p = NULL; }
              }
            free(buf->rowptr); buf->rowptr = NULL; 
          }
        free(buf); buf = NULL; 
      }
  }

double **float_image_buffer_alloc_rows(int NC, int NX, int NB, bool_t clear)
  { double **rowptr = (double **)notnull(malloc(NB*sizeof(double *)), "no mem");
    int NCX = NC*NX; /* Samples per row. */
    int yb;
    for (yb = 0; yb < NB; yb++) 
      { double *p = (double *)notnull(malloc(NCX*sizeof(double)), "no mem");
        rowptr[yb] = p;
        if (clear) { int k; for (k = 0; k < NCX; k++, p++) { (*p) = 0.0; } }
      }
    return rowptr;
  }

int float_image_buffer_row_pos(float_image_buffer_t *buf, int y)
  { 
    if (y < buf->yini) 
      { return -1; }
    else if (y >= buf->ylim) 
      { return +1; }
    else
      { return 00; }
  }

double *float_image_buffer_get_row(float_image_buffer_t *buf, int y)
  { 
    if(float_image_buffer_row_pos(buf,y) == 00)
      { return buf->rowptr[y % (buf->NB)]; }
    else
      { return NULL; }
  }
  
void float_image_buffer_fill_row(float_image_buffer_t *buf, int y, double val)
  { 
    double *row = float_image_buffer_get_row(buf, y);
    int NCX = buf->sz[0]*buf->sz[1];
    int k;
    for (k = 0; k < NCX; k++) { row[k] = val; }
  }

void float_image_buffer_advance(float_image_buffer_t *buf)
  { 
    int y = buf->ylim;
    demand(y < buf->sz[2], "no next row");
    buf->ylim++;
    if (buf->ylim - buf->yini > buf->NB) { buf->yini++; }
    assert(buf->ylim - buf->yini <= buf->NB);
  }
