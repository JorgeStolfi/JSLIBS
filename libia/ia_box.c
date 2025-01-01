/* See box.h. */
/* Last edited on 2024-12-31 00:52:04 by stolfi */

#include <math.h>

#include <ia.h>
#include <flt.h>

#include <ia_box.h>

void ia_box_corner(ia_box_Dim d, ia_box_Interval *b, ia_box_Corner c, Float *x)
{
  for (int32_t i = 0; i < d; i ++)
    { x[i] = (c & 1 ? b[i].hi : b[i].lo);
      c  = (ia_box_Corner)(c << 1);
    }
}

void ia_box_center(ia_box_Dim d, ia_box_Interval *b, Float *x)
{
  for (int32_t i = 0; i < d; i ++) { x[i] = ia_mid(b[i]); }
}

ia_box_Axis ia_box_max_axis(ia_box_Dim d, ia_box_Interval *b)
{
  int32_t imax = 0;
  Float dmax = MinusInfinity;
  for (int32_t i = 0; i < d; i ++) 
    { Float ri = ia_rad(b[i]); 
      if (ri > dmax) { dmax = ri; imax = i; }
    }
  return (ia_box_Axis)imax;
}

Float ia_box_max_rad(ia_box_Dim d, ia_box_Interval *b)
{
  Float dmax = MinusInfinity;
  for (int32_t i = 0; i < d; i ++) 
    { Float ri = ia_rad(b[i]); 
      if (ri > dmax) { dmax = ri; }
    }
  return dmax;
}

double ia_box_radius(ia_box_Dim d, ia_box_Interval *b)
{
  double r, r2 = Zero;
  ROUND_UP;
  for (int32_t i = 0; i < d; i ++) 
    { double ri = (double)(ia_rad(b[i])); 
      r2 += ri*ri;
    }
  r = sqrt(r2);
  ROUND_NEAR;
  return r;
}

void ia_box_from_corners(ia_box_Dim d, Float *x, Float *y, ia_box_Interval *b)
{
  for (int32_t i = 0; i < d; i ++) 
    { Float xi = x[i];
      Float yi = y[i];
      if (xi <= yi) 
        { b[i].lo = xi; b[i].hi = yi; }
      else
        { b[i].lo = yi; b[i].hi = xi; }
    }
}

void ia_box_split(ia_box_Dim d, ia_box_Axis k, ia_box_Dir dir, ia_box_Interval *b, ia_box_Interval *h)
{
  for (int32_t i = 0; i < d; i ++) 
    { ia_box_Interval *hp = &(h[i]);
      *hp = b[i];
      if (i == k)
        { Float md = ia_mid(*hp);
          if (dir != HI) { hp->hi = md; }
          if (dir != LO) { hp->lo = md; }
        }
    }
}

 
