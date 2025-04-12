
void tspo_find_closest_pair
  ( int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double spotDist,
    uint32_t *kq_P,
    uint32_t queue[],
    uint32_t *icMin_P, uint32_t *jcMin_P,
    double *drMin_P
  );
  /* Let {kq} be the value of {*kq_P} on entry. For all distinct spot
    index pairs indices {ic,jc} in {queue[kq..NS-1]} let {dr} be the
    ratio {d/(srad[ic]+srad[jc])}. Finds the minimum {drMin} of all
    those ratios, and the pair {icMin,jcMin} that realizes it. Assumes
    that all centers are reduced modulo {NX,NY}.
    
    Returns the results in {*icMin_P,*jcMin_P,*drMin_P}. Also moves
    those two entries in {queue} to {queue[kq]} and {queue[kq+1]} and
    increments {*kq_P} by 2. */
