/* Test function for 2d plots */

Float foit_f (Float x, Float y)
  {
    ROUND_NEAR;
    return (x*x + y*y + x*y - (x*y)*(x*y)*Half - Quarter);
  }

Interval foit_fv (Interval x, Interval y)
  {
    Interval x2 = iv_sqr(x);
    Interval y2 = iv_sqr(y);
    Interval xy = iv_mul(x, y);
    Interval x2y2d2 = iv_affine(iv_sqr(xy), Half, Zero, Zero);
    Interval sum = iv_add(iv_add(iv_add(x2, y2), xy), iv_neg(x2y2d2));
    Interval res = iv_affine (sum, One, -Quarter, Zero);
    return (res);
  }

FOIP foit_ff (FOIP x, FOIP y)
  {
    MemP frame = foi_top();
    FOIP x2  = foi_sqr(x);
    FOIP y2  = foi_sqr(y);
    FOIP xy  = foi_mul(x, y);
    FOIP x2y2d2 = foi_affine(foi_sqr(xy), Half, Zero, Zero);
    FOIP sum = foi_add(foi_add(foi_add(x2, y2), xy), foi_neg(x2y2d2));
    FOIP res = foi_affine(sum, One, -Quarter, Zero);
    return (foi_return(frame, res));
  }

char *foit_fname = "f = x^2 + y^2 + x*y - (x*y)^2/2 - 1/4";

Interval foit_fxd = {-Two, Two};
Interval foit_fyd = {-Two, Two};

int foit_fn = 32;
