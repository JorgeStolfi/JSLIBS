/* Last edited on 2022-02-28 12:52:45 by stolfi */

double *hr2_pmap_aff_elem_address(hr2_pmap_t *M, int32_t s);
  /* Returns the address of element {s} of the affine projective map whose direct matrix is {A}.
    The elements {A.c[0..2][0]} must be {1,0,0}.
    The integer {s} must be in {0..5}.  For {s} is in {0..3},
    returns the addresses of the elements of the linear 2x2 subarray {A.c[1..2][1..2]}, in row-by roe order.
    For {s} in {4..5}, returns the addresses of the displacement coefficients {A.c[0][1..2]}. */

double *hr2_pmap_aff_elem_address(hr2_pmap_t *M, int32_t s)
  { r3x3_t *A = &(M->dir);
    hr2_pmap_aff_norm_check(A);
    demand((s >= 0) && (s < 6), "invalid elem index");
    if (s < 4) 
      { int32_t i = 1 + s / 2, j = 1 + s % 2;
        return &(A->c[i][j]);
      }
    else
      { int32_t j = 1 + (s - 4);
        return &(A->c[0][j]);
      }
  }
