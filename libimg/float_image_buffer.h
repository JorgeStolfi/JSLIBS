#ifndef float_image_buffer_H
#define float_image_buffer_H

/* Multipurpose buffer for multiple rows of an arbitrary image. */
/* Last edited on 2010-08-04 15:34:43 by stolfi */

#include <bool.h>

typedef struct float_image_buffer_t
  { int sz[3];          /* Number of image. */
    int NB;             /* Number of rows in buffer. */
    int yini;           /* Index of first image row currently in buffer. */
    int ylim;           /* Index of first image row after those in buffer. */
    double **rowptr;    /* Pointers to rows of pixels. */
  } float_image_buffer_t;
  /* A {float_image_buffer_t} stores a horizontal band from an image
    with {sz[0]} channels, {sz[1]} columns, and {sz[2]} rows. The
    buffer stores a band of consecutive rows with indices
    {yini..ylim-1}, a subset of {0..sz[2]-1}, for a maximum of
    {NB} rows. The buffer is empty when {yini==ylim}.
    
    The address {p} of row {y}, for {y} in that range, is stored in
    {rowptr[y%NB]}. The pixel with channel {c} and abscissa {x}
    in that row is stored in {p[sz[0]*x+c]} for {c} in {0..sz[0]-1}
    and {x} in {0..sz[1]-1}. */

float_image_buffer_t *float_image_buffer_new(int NC, int NX, int NY, int NB);
  /* Creates a new buffer structure {buf} for an image with {NC}
    channels, {NX} columns, and {NY} rows with space for at most
    {NB} rows. Allocates the the array {buf->rowptr}.  Initially
    {buf->yini==buf->ylim==0}, i.e. the buffer is empty and positioned
    just before row 0. */
    
int float_image_buffer_row_pos(float_image_buffer_t *buf, int y);
  /* Returns the position of row {y} relative to the buffered rows.
    Namely returns {-1} if {y < buf->yini}, {+1} if {y >= buf->ylim},
    and 0 if {y} is in {buf->yini..buf->ylim-1}. */
    
double *float_image_buffer_get_row(float_image_buffer_t *buf, int y);
  /* If {y} is in the range {buf->yini..buf->ylim-1}, returns 
    the address of the sample vector for row {y}, namely
    {buf->rowptr[y%buf->NB]}.  Otherwise returns NULL. */

void float_image_buffer_fill_row(float_image_buffer_t *buf, int y, double val);
  /* If {y} is in the range {buf->yini..buf->ylim-1}, stores {val} into all
    elements of the sample vector for row {y}, namely
    {buf->rowptr[y%buf->NB]}.  Fails otherwise. */

void float_image_buffer_advance(float_image_buffer_t *buf);
  /* Makes space in the buffer for the next image row, namely row 
    {y=buf->ylim}, and increments {buf->ylim}.
    
    If the buffer already contains {buf->NB} rows, increments
    {buf->yini}, thus implicitly deallocating row {buf->yini}. Fails
    if {buf->ylim >= buf->sz[2]}.
    
    The indices {buf->yini} and {buf->ylim} are incremented even when
    {buf->NB} is zero, but in that case (and only that case) row {y}
    will not be in the buffer. */

void float_image_buffer_free(float_image_buffer_t *buf);
  /* Frees all storage associated with {buf}, icluding {*buf} itself. */

#endif
