/* Test function for 2d plots */

Float foit_f (Float x, Float y)
  {
    ROUND_NEAR;
    {
      Float x2 = x*x;
      Float y2 = y*y;
      Float r2 = x2 + y2;
      Float m = x2 - Three * y2;
      Float m2 = m*m;
      Float res = x2*m2 - r2;
      return (res);
    }
  }

Interval foit_fv (Interval x, Interval y)
  {
    Interval x2 = iv_sqr(x);
    Interval y2 = iv_sqr(y);
    Interval r2 = iv_add(x2, y2);
    Interval m = iv_sub(x2, iv_affine(y2, Three, Zero, Zero));
    Interval m2 = iv_sqr(m);
    Interval res = iv_sub(iv_mul(x2, m2),  r2);
    return (res);
  }

FOIP foit_ff (FOIP x, FOIP y)
  {
    MemP frame = foi_top();
    FOIP x2 = foi_sqr(x);
    FOIP y2 = foi_sqr(y);
    FOIP r2 = foi_add(x2, y2);
    FOIP m = foi_sub(x2, foi_affine(y2, Three, Zero, Zero));
    FOIP m2 = foi_sqr(m);
    FOIP res = foi_sub(foi_mul(x2, m2),  r2);
    return (foi_return(frame, res));
  }

char *foit_fname = 
  "r2 = x^2 + y^2;\nm2 = (x^2 - 3*y^2)^2;\nf = x^2*m^2 - r^2";

Interval foit_fxd = {-Three, Three};
Interval foit_fyd = {-Three, Three};

int foit_fn = 32;
