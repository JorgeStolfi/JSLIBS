/* Test function for 2d plots */

Float foit_g (Float x)
  {
    ROUND_NEAR;
    {
      Float y = -x;
      Float d = x - y;
      Float t = Quarter * d * (One - d*d);
      Float i = One / (One + t*t);
      Float h = t * i;
      Float res = i;
      return (res);
    }
  }

Interval foit_gv (Interval x)
  {
    Interval y = iv_neg(x);
    Interval d = iv_sub(x, y);
    Interval t = iv_mul(d, iv_affine(iv_sqr(d), -Quarter, Quarter, Zero));
    Interval i = iv_inv(iv_affine(iv_sqr(t), One, One, Zero));
    Interval h = iv_mul(t, i);
    Interval res = i;
    return (res);
  }

FOIP foit_gf (FOIP x)
  {
    MemP frame = foi_top();
    FOIP y = foi_from_interval(iv_neg(foi_range(x)));
    FOIP d = foi_sub(x, y);
    FOIP t = foi_mul(d, foi_affine(foi_sqr(d), -Quarter, Quarter, Zero));
    FOIP i = foi_inv(foi_affine(foi_sqr(t), One, One, Zero));
    FOIP h = foi_mul(t, i);
    FOIP res = i;
    return (foi_return(frame, res));
  }

char *foit_gname = 
  "y = -range(x); d = x-y; t = d(1-d^2); i = 1/(1 + t^2); h = t*i; g = i";

Interval foit_gxd = {-Two, Two};
Interval foit_gyd = {-Two, Two};

int foit_gn = 8;
