
void tspo_find_closest_pair
  ( int32_t NX, int32_t NY,
    uint32_t NS,
    r2_t sctr[], double srad[], double spotDist,
    uint32_t *kq_P,
    uint32_t queue[],
    uint32_t *icMin_P, uint32_t *jcMin_P,
    double *drMin_P
  )
  {
    uint32_t kq = (*kq_P);
    uint32_t iqMin, jqMin;  /* Indices of pair in {queue}. */
    uint32_t icMin, jcMin;  /* Indices of pair in {sctr}. */
    double drMin = +INF;
    for (uint32_t jq = kq; jq < NS; jq++)
      { uint32_t jc = queue[jq];
        r2_t *cj = &(sctr[jc]);
        double xj = cj->c[0], yj = cj->c[1];
        assert((xj >= 0) && (xj < NX) && (yj >= 0) && (yj < NY));
        /* Compare with other vertices: */
        for (uint32_t iq = kq; iq < jq; iq++)
          { uint32_t ic = queue[iq];
            r2_t *sci = &(sctr[ic]);
            double xi = sci->c[0], yi = sci->c[1];
            /* Center {sci} must be in domain per previous assert. */
            /* Tries to take into account the toroidal geometry: */
            for (int32_t ux = -1; ux <= +1; ux++)
              { for (int32_t uy = -1; uy <= +1; uy++)
                  { double dij = hypot(xi - xj + ux*NX, yi - yj + uy*NY);
                    double drij = dij/(srad[ic]+srad[jc]);
                    if (drij < drMin)
                      { drMin = drij;
                        jqMin = jq; iqMin = iq;
                        jcMin = jc; icMin = ic;
                      }
                  }
              }
          }
      }
    /* Bring {quque[jqMin], queue[iqMin]} to front: */
    assert(iqMin < jqMin);
    assert(iqMin >= kq);
    /* Must fix in this order, {iqMin} before {jqMin}: */
    if (iqMin > kq) { queue[iqMin] = queue[kq]; queue[kq] = icMin; }
    kq++;
    assert(jqMin >= kq);
    if (jqMin > kq) { queue[jqMin] = queue[kq]; queue[kq] = jcMin; }
    kq++;
    
    (*icMin_P) = icMin;
    (*jcMin_P) = jcMin;
    (*drMin_P) = drMin;
    (*kq_P) = kq;
  }
