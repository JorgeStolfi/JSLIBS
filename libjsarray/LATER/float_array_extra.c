/* Last edited on 2023-03-18 11:03:56 by stolfi */

/* INTERNAL PROTOTYPES */

void float_array_get_pixel(float_array_t *A, int32_t x, int32_t y, float v[])
  { float *sp = float_array_get_elem_address(A, 0, x, y);
    int32_t NC = A->sz[0];
    int32_t c;
    for (c = 0; c < NC; c++) { v[c] = (*sp); sp += A->st[0]; }
  }

void float_array_set_pixel(float_array_t *A, int32_t x, int32_t y, float v[])
  { float *sp = float_array_get_elem_address(A, 0, x, y);
    int32_t NC = A->sz[0];
    int32_t c;
    for (c = 0; c < NC; c++) { (*sp) = v[c]; sp += A->st[0]; }
  }

void float_array_fill_pixel(float_array_t *A, int32_t x, int32_t y, float v)
  { float *sp = float_array_get_elem_address(A, 0, x, y);
    int32_t NC = A->sz[0];
    int32_t c;
    for (c = 0; c < NC; c++) { (*sp) = v; sp += A->st[0]; }
  }


void float_array_get_elem_row(float_array_t *A, int32_t c, int32_t y, float v[])
  { float *sp = float_array_get_elem_address(A, c, 0, y);
    int32_t NX = A->sz[1];
    int32_t x;
    for (x = 0; x < NX; x++) { v[x] = (*sp); sp += A->st[1]; }
  }

void float_array_set_elem_row(float_array_t *A, int32_t c, int32_t y, float v[])
  { float *sp = float_array_get_elem_address(A, c, 0, y);
    int32_t NX = A->sz[1];
    int32_t x;
    for (x = 0; x < NX; x++) { (*sp) = v[x]; sp += A->st[1]; }
  }

void float_array_fill_elem_row(float_array_t *A, int32_t c, int32_t y, float v)
  { float *sp = float_array_get_elem_address(A, c, 0, y);
    int32_t NX = A->sz[1];
    int32_t x;
    for (x = 0; x < NX; x++) { (*sp) = v; sp += A->st[1]; }
  }


void float_array_get_pixel_row(float_array_t *A, int32_t y, float v[])
  { float *pxp = float_array_get_elem_address(A, 0, 0, y);
    int32_t NC = A->sz[0];
    int32_t NX = A->sz[1];
    int32_t c, x; int32_t k = 0;
    for (x = 0; x < NX; x++)
      { float *sp = pxp;
        for (c = 0; c < NC; c++) 
          { v[k] = (*sp); sp += A->st[0]; k++; }
        pxp += A->st[1];
      }
  }

void float_array_set_pixel_row(float_array_t *A, int32_t y, float v[])
  { float *pxp = float_array_get_elem_address(A, 0, 0, y);
    int32_t NC = A->sz[0];
    int32_t NX = A->sz[1];
    int32_t c, x; int32_t k = 0;
    for (x = 0; x < NX; x++)
      { float *sp = pxp;
        for (c = 0; c < NC; c++) 
          { (*sp) = v[k]; sp += A->st[0]; k++; }
        pxp += A->st[1];
      }
  }

void float_array_fill_pixel_row(float_array_t *A, int32_t y, float v)
  { float *pxp = float_array_get_elem_address(A, 0, 0, y);
    int32_t NC = A->sz[0];
    int32_t NX = A->sz[1];
    int32_t c, x;
    for (x = 0; x < NX; x++)
      { float *sp = pxp;
        for (c = 0; c < NC; c++) 
          { (*sp) = v; sp += A->st[0]; }
        pxp += A->st[1];
      }
  }

void float_array_fill_channel(float_array_t *A, int32_t c, float v)
  { int32_t NC = A->sz[0];
    int32_t NX = A->sz[1]; 
    int32_t NY = A->sz[2];
    demand((c >= 0) && (c < NC), "invalid channel");
    int32_t x, y;
    for (y = 0; y < NY; y++)
      { for (x = 0; x < NX; x++)
          { float_array_set_elem(A, c, x, y, v); }
      }
  }

void float_array_fill(float_array_t *A, float v)
  { int32_t NC = A->sz[0];
    int32_t c;
    for (c = 0; c < NC; c++) { float_array_fill_channel(A, c, v); }
  }

/* IMAGE INTERPOLATION */

float float_array_interpolate_elems(float_array_t *A, int32_t c, double x, double y)
  { /* Get pixel addresses and weights: */
    double wx0, wx1;    /* Column weights. */
    double wy0, wy1;    /* Sample weights. */
    float *p00, *p10, *p01, *p11; /* Sample addresses. */
    float_array_get_interpolation_data
      (A, c, x, y, &p00, &p10, &p01, &p11, &wx0, &wx1, &wy0, &wy1);
    
    /* Apply interpolation formula: */
    return ((*p00)*wx0 + (*p10)*wx1)*wy0 + ((*p01)*wx0 + (*p11)*wx1)*wy1;
  }
  
void float_array_interpolate_pixels(float_array_t *A, double x, double y, float v[])
  { /* Get pixel addresses and weights: */
    double wx0, wx1;    /* Column weights. */
    double wy0, wy1;    /* Sample weights. */
    float *p00, *p10, *p01, *p11; /* Sample addresses. */
    float_array_get_interpolation_data
      (A, 0, x, y, &p00, &p10, &p01, &p11, &wx0, &wx1, &wy0, &wy1);
    
    /* Apply interpolation formula to all channels: */
    int32_t NC = A->sz[0];         /* Number of channels. */
    ix_step_t cst = A->st[0];  /* Position increment between channels. */
    int32_t c;
    for (c = 0; c < NC; c++)
      { /* Interpolate: */
        v[c] = ((*p00)*wx0 + (*p10)*wx1)*wy0 + ((*p01)*wx0 + (*p11)*wx1)*wy1;
        /* Advance sample pointers to next channel: */
        p00 += cst; p10 += cst;  p01 += cst;  p11 += cst;
      }
  }

void float_array_get_interpolation_data
  ( float_array_t *A,
    int32_t c,
    double x, 
    double y,
    float **p00,
    float **p10,
    float **p01,
    float **p11,
    double *wx0,
    double *wx1,
    double *wy0, 
    double *wy1
  )
  {
    int32_t ix0, ix1;       /* Column indices to interpolate. */
    float_array_get_interp_indices_and_weights(x, A->sz[1], &ix0, &ix1, wx0, wx1);
    int32_t dx = ix1 - ix0; /* Usually 1, may be 0. */
    
    int32_t iy0, iy1;       /* Row indices to interpolate. */
    float_array_get_interp_indices_and_weights(y, A->sz[2], &iy0, &iy1, wy0, wy1);
    int32_t dy = iy1 - iy0; /* Usually 1, may be 0. */
    
    /* Get sample addresses: */
    (*p00) = float_array_get_elem_address(A, c, ix0, iy0);
    (*p10) = (*p00) + dx*A->st[1];
    (*p01) = (*p00) + dy*A->st[2];
    (*p11) = (*p10) + dy*A->st[2];
  }
  
void float_array_get_interp_indices_and_weights
  ( double z, 
    int32_t N, 
    int32_t *i0, 
    int32_t *i1, 
    double *w0,
    double *w1
  )
  {
    /* Get the `low' index {i0}: */
    (*i0) = (int32_t)floor(z - 0.5);
    
    /* Compute the interpolation weights {w0} and {w1}: */
    (*w1) = z - (*i0) - 0.5;
    if ((*w1) < 0.0) { (*i0)--; (*w1) += 1.0; }
    if ((*w1) > 1.0) { (*i0)++; (*w1) -= 1.0; }
    (*w0) = 1.0 - (*w1);
   
    /* Compute the superior index {i1}: */
    (*i1) = (*i0) + 1;

    /* Clip {i0} and {i1} to {0..N-1}: */
    if ((*i0) < 0) 
      { (*i0) = (*i1) = 0; }
    else if ((*i1) >= N)
      { (*i0) = (*i1) = N-1; }
  }

void float_array_update_elem_range(float_array_t *A, int32_t c, float *vmin, float *vmax)
  { int32_t NC = A->sz[0];
    int32_t NX = A->sz[1];
    int32_t NY = A->sz[2];
    if ((c < 0) || (c >= NC))
      { /* Invalid channel, all samples are zero. (But there may be no samples!) */
        if ((NX > 0) && (NY > 0))
          { float v = 0.0;
            if (v > (*vmax)) { (*vmax) = v; }
            if (v < (*vmin)) { (*vmin) = v; }
          }
      }
    else
      { /* Valid channel. */
        int32_t x, y;
        for (y = 0; y < NY; y++)
          { for (x = 0; x < NX; x++)
              { float v = float_array_get_elem(A, c, x, y);
                if (v > (*vmax)) { (*vmax) = v; }
                if (v < (*vmin)) { (*vmin) = v; }
              }
          }
      }
  }

#define float_array_max_elems (1024u*1024u*1024u)
  /* Maximum total elements (2^30), for sanity checks. */

void float_image_get_interpolation_data
  ( float_image_t *fim,
    int32_t c,
    double x, 
    double y,
    float **p00,
    float **p10,
    float **p01,
    float **p11,
    double *wx0,
    double *wx1,
    double *wy0, 
    double *wy1
  );
  /* Computes the addresses {*p00,*p01,*p10,*p11} of the four samples
    of channel {c} that are involved in the interpolation of image {fim}
    at the point {(x,y)}, and the corresponding weights.

    The first digit in the pointer's name is the X (column) increment,
    the second index is the row increment. Therefor, {*p00} is the
    pixel with lowest X and Y indices among the four; {*p10} is its
    neighbor in the next column of the same row; and so on.

    The pixel addresses are clipped against the number of columns {NX}
    and the number of rows {NY}. Namely, if {(x,y)} lies outside the
    domain rectangle {[0_NX] × [0_NY]}, then the addresses
    {*p00,*p01,*p10,*p11} get adjusted so that they point to valid
    surrogate samples. These adjustment have the same effect as
    extending the image {fim} to infinity, by duplicating its border
    rows and columns. */

void float_image_get_interp_indices_and_weights
  ( double z, 
    int32_t N, 
    int32_t *i0, 
    int32_t *i1, 
    double *w0,
    double *w1
  );
  /* Computes the indices {i0,i1} of the pixels nedded to interpolate
    a row (or column) of samples at the fractional coordinate {z}, and
    their respective weights {w0,w1}. The indices are consecutive,
    except that they are clipped to the range {0..N-1}; therefore
    {i1 - i0} is either 0 or 1. */

