/* Last edited on 2025-01-11 14:55:47 by stolfi */

void pst_slope_map_build_system_2(float_image_t *IG, float_image_t *IW, imgsys_t* S)
  {

    /* Get/check the sizes of the slope maps: */
    int NC = IG->sz[0]; assert(NC == 2);
    int NX = IG->sz[1];
    int NY = IG->sz[2];
    
    if (IW != NULL)
      { assert(IW->sz[0] == 1);
        assert(IW->sz[1] == NX);
        assert(IW->sz[2] == NY);
      }
    
    /* Check the size of the system: */
    int NH = (NX+1)*(NY+1); /* Number of unknowns (heights). */
    assert(S->N == NH);
    
    /* Conceptually, we define two quadratic mismatch terms
      {qx[x,y](h)} and {qy[x,y](h)} for each pixel in the grid.
      
      The term {qx[x,y](h)} is the square of the difference between
      the mean slope in the pixel {[x,y]}, as
      obtained from the slope map {IG}, and the difference of the 
      heights at the midpoints of the two vertical sides of the pixel,
      obtained by interpolating the (unknown) corner heights.
      
        { qx[x,y](h) = (dx[x,y] - (- z00 + z10 - z01 + z11)/2)^2 }

      where {z00}, {z10}, {z01}, and {z10} are the {h} variables
      corresponding to heights {z[x+ix,y+iy]}, where {(ix,iy)} is
      {(0,0)}, {(1,0)}, {(0,1)}, and {(1,1)}, respectively; and
      {dx[x,y]} is the X-slope in the pixel, namely {dx[x,y] =
      IG[0,x,y]}.
      
      The term {qy[x,y](h)} is similarly defined, using the edges
      parallel to the Y axis:
      
        { qy[x,y](h) = (dy[x,y] - (- z00 - z01 + z10 + z11)/2)^2 }
       
      where {dy[x,y] = IG[1,x,y]}.

      The terms {qx[x,y](h)} and {qy[x,y](h)} exists for all {x} in
      {0..NX-1} and all {y} in {0..NY-1}.
      
      We want to minimize the sum {Q(h)} of all terms {qx[x,y](h)} and 
      {qy[x,y](h)}, each weighted by {W[x,y]}.

      Note that if the slope map {IG} is zero, we can anihilate all
      those terms (and therefore Q) by setting (1) all
      heights equal to 1, or (2) all weights alternating +1 or -1 in
      checkerboard fashion. In general, adding solutions (1) or (2) to
      any height map will not change {Q}. Therefore, if the system is
      to be solved iteratively, one must take care to exclude these
      two solutions.
      
      Minimizing {Q} is equivalent to finding the {h} vector that
      anihilates the gradient {G} of {Q}. The gradient of {Q} is an
      affine function of {h}, namely {G(h) = A h - b}.
      
      Each height {z[x,y]} appears into four {qx} terms and four {qy}
      terms, corresponding to the four pixels that share the corner
      {(x,y)}. Those terms are a quadratic function of nine height
      values, {z[x+ix,y+iy]} for all {ix,iy} in {-1,0,+1}. Therefore,
      the derivative of {Q} relative to {z[x,y]} (the "force" acting on
      {z[x,y]}) is an affine function
      {SUM{A[r,s]*h[s] : s} + b[r]} of {h} that depends on those nine
      elements of {h} only.   Actually it turns out that the effects of 
      four of those elements --- those on the same row or column as 
      {z[x,y]} --- cancel out. So, each row of the {A} matrix has
      only five nonzero elements, at most.
    */
    
    auto long int ind_h(int x, int y); 
      /* This function computes the index into the solution vector {h}
        corresponding to vertex {[x,y]} of the images. */
    
    long int ind_h(int x, int y) { return x + y*(NX+1); }
    
    /* Make sure that we can build the system: */
    assert(MAX_COEFFS >= 5);

    /* Enumerate the vertices of {IDX,IDY}, i.e. elems of the height map: */
    int x, y;
    for(y = 0; y <= NY; y++)
      { for(x = 0; x <= NX; x++)
          { /* Get the slopes and relative weights of the pixels incident to {(x,y)}: */
            double dxmm = 0.0, dxmp = 0.0, dxpm = 0.0, dxpp = 0.0;
            double dymm = 0.0, dymp = 0.0, dypm = 0.0, dypp = 0.0;
            double wmm = 0.0, wmp = 0.0, wpm = 0.0, wpp = 0.0;

            if ((x > 0) && (y > 0))
              { dxmm = get_sample(IG, 0, x-1, y-1);
                dymm = get_sample(IG, 1, x-1, y-1);
                wmm = get_sample(IW, 0, x-1, y-1);
                assert(wmm >= 0.0);
              }

            if ((x > 0) && (y < NY))
              { dxmp = get_sample(IG, 0, x-1, y);
                dymp = get_sample(IG, 1, x-1, y);
                wmp = get_sample(IW, 0, x-1, y);
                assert(wmp >= 0.0);
              }

            if ((x < NX) && (y > 0))
              { dxpm = get_sample(IG, 0, x, y-1);
                dypm = get_sample(IG, 1, x, y-1);
                wpm = get_sample(IW, 0, x, y-1);
                assert(wpm >= 0.0);
              }

            if ((x < NX) && (y < NY))
              { dxpp = get_sample(IG, 0, x, y);
                dypp = get_sample(IG, 1, x, y);
                wpp = get_sample(IW, 0, x, y);
                assert(wpp >= 0.0);
              }
           
            /* Compute the weighted diagonal increments {wdmm,wdpm,wdmp,wdpp}: */
            double wdmm = 0.0, wdpm = 0.0, wdmp = 0.0, wdpp = 0.0;
            
            /* Make sure that the diagonal coefficient is nonzero: */
            double wtot = wmm + wmp + wpm + wpp;
            if (wtot == 0.0)
              { /* Height {z[x,y]} is loose. Make it be the average of its cross neighbors: */
                wmm = ((x <= 0)  || (y <= 0)  ? 0.0 : 1.0);
                wmp = ((x <= 0)  || (y >= NY) ? 0.0 : 1.0);
                wpm = ((x >= NX) || (y <= 0)  ? 0.0 : 1.0);
                wpp = ((x >= NX) || (y >= NY) ? 0.0 : 1.0);
                wtot = wxm + wxp + wym + wyp;
                /* Leave {wdmm,wdpm,wdmp,wdpp} all zero. */
              }
            else
              { /* Compute {wdmm,wdpm,wdmp,wdpp} from given slopes: */
                wdmm = wmm*(- dxmm - dymm);
                wdpm = wpm*(+ dxpm - dypm);
                wdmp = wmp*(- dxmp + dymp);
                wdpp = wpp*(+ dxpp + dypp);
              }
            assert(wtot > 0.0);
                
            /* Get the index {k} into {h} of {IZ[x,y]} and its diagonal neighbors: */
            long int koo = ind_h(x, y); 
            long int kmm = ind_h(x-1, y-1); 
            long int kpm = ind_h(x+1, y-1); 
            long int kmp = ind_h(x-1, y+1); 
            long int kpp = ind_h(x+1, y+1); 
            
            /* Assemble the equation for height {z[x,y]}, divided by {-wtot}: */
            imgsys_equation_t *eqk = &(S->eq[k]);
            int nt = 0; /* Number of terms in equation. */
            eqk->rhs = 0.0;
            eqk->uid[nt] = k;  eqk->cf[nt] = 1.00; nt++;
            if (wmm != 0.0) 
              { eqk->uid[nt] = kmm; eqk->cf[nt] = -wmm/wtot; eqk->rhs += -wdmm/wtot; nt++; }
            if (wpm != 0.0) 
              { eqk->uid[nt] = kpm; eqk->cf[nt] = -wpm/wtot; eqk->rhs += -wdpm/wtot; nt++; }
            if (wmp != 0.0) 
              { eqk->uid[nt] = kmp; eqk->cf[nt] = -wmp/wtot; eqk->rhs += -wdmp/wtot; nt++; }
            if (wpp != 0.0) 
              { eqk->uid[nt] = kpp; eqk->cf[nt] = -wpp/wtot; eqk->rhs += -wdpp/wtot; nt++; }
            eqk->nt = nt;
          }
      }
  }
