/* Last edited on 2009-08-25 19:54:23 by stolfi */

  { 
    /* !!! Must do all drawing AFTER, otherwise the fill will overlap the draw at endpoint. !!! */
    if (draw)
      { if ((i0 > 0) && (i1 == 0)) 
          { /* Draw the top side of tile (bottom edge of quad): */
            pswr_segment(ps, f01[0], f01[1], f11[0], f11[1]);
          }
        if ((i0 == ns) && (i1 > 0))
          { /* Draw the right side of tile (right edge of quad): */
            pswr_segment(ps, f10[0], f10[1], f11[0], f11[1]);
          }
        if ((i1 == ns) && (i0 > 0))
          { /* Draw the top side of tile (top edge of quad): */
            pswr_segment(ps, f11[0], f11[1], f01[0], f01[1]);
          }
        if ((i0 == 0) && (i1 > 0))
          { /* Draw the right side of tile (left edge of quad): */
            pswr_segment(ps, f11[0], f11[1], f10[0], f10[1]);
          }
      }
  }   
