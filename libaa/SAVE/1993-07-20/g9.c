/* Test function for single-variable graph */

Float foit_h (Float x);
Interval foit_hv (Interval x);
FOIP foit_hf (FOIP x);

Float foit_h (Float x)
  {
    ROUND_NEAR;
    return (sqrt(x*x - x + 0.50)/sqrt(x*x + 0.50));
  }

Float foit_g (Float x)
  {
    ROUND_NEAR;
    return foit_h(foit_h(x));
  }

Interval foit_hv (Interval x)
  {
    Interval x2  = iv_sqr(x);
    Interval hlf = {Half, Half};
    Interval dif = iv_sub(x2, x);
    Interval sum1 = iv_add(dif, hlf);
    Interval sum2 = iv_add(x2, hlf);
    Interval res = iv_mul(iv_sqrt(sum1), iv_inv(iv_sqrt(sum2)));
    return (res);
  }

Interval foit_gv (Interval x)
  {
    return (foit_hv(foit_hv(x)));
  }

FOIP foit_hf (FOIP x)
  {
    FOIP x2  = foi_sqr(x);
    FOIP hlf = foi_const(Half, Zero);
    FOIP dif = foi_sub(x2, x);
    FOIP sum1 = foi_add(dif, hlf);
    FOIP sum2 = foi_add(x2, hlf);
    FOIP res = foi_mul(foi_sqrt(sum1), foi_inv(foi_sqrt(sum2)));
    return (res);
  }

FOIP foit_gf (FOIP x)
  {
    MemP frame = foi_top();
    FOIP res = foit_hf(foit_hf(x));
    return (foi_return(frame, res));
  }

char *foit_gname = "h(x) = sqrt(x^2 - x + 1/2)/sqrt(x^2 + 1/2) g(x) = h(h(x))";

Interval foit_gxd = {-Two, Two};
Interval foit_gyd = {-Two, Two};

int foit_gn = 64;
