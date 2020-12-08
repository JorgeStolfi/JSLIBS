/* Last edited on 2013-03-26 23:39:05 by stolfilocal */

void pst_adjust_weights
  (
    int n,
    int d,
    double x[],      /* IN: {x[i*d+k]} is coord {k} of observation {i}. */
    double V[],      /* IN/OUT: {V[i*d+k]} is the noise variance in {x[i*d+k]}. */
    double p_bad[],  /* IN/OUT: {p_bad[i]} is the prob of obs {i} being outlier. */
    double x_bad[],  /* IN: {x_bad[k]} is the mean value of {x[i][k]} if {x[i]} is outlier. */
    double V_bad[],  /* IN: {V_bad[k]} is the noise variance in {x[i][k]} if {x[i]} is outlier. */
    double x_gud[],  /* IN/OUT: {x_gud[k]} is the estimated coord {k} of the true vector. */
    double V_gud[]   /* IN/OUT: {V_gud[k]} is the estimated variance of noise in {x_gud[k]}. */
  )
  {
    return;
    /* !!! To write: !!! */
    // double *V_tmp = 
    // p_fac = 0; /* Trust factor in inlier probs. */
    // do
    //   { 
    //     /* Compute the variances from probs */
    //     /* Compute average */
    //     /* Compute probabilities */
    //   }
    // while (! converged);
    // 
    // /* Copy variances out */
  }
  
