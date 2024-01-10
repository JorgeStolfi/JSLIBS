/* Test function for single-variable graph */

Float foit_g (Float x)
  {
    ROUND_NEAR;
    return (One/sqrt(x*x + Half));
  }

Interval foit_gv (Interval x)
  {
    Interval x2  = iv_sqr(x);
    Interval hlf = {Half, Half};
    Interval sum = iv_add(x2, hlf);
    Interval res = iv_inv(iv_sqrt(sum));
    return (res);
  }

FOIP foit_gf (FOIP x)
  {
    MemP frame = foi_top();
    FOIP x2  = foi_sqr(x);
    FOIP hlf = foi_const(Half, Zero);
    FOIP sum = foi_add(x2, hlf);
    FOIP res = foi_inv(foi_sqrt(sum));
    return (foi_return(frame, res));
  }

char *foit_gname = "g(x) = 1/sqrt(x^2 + 1/2)";

Interval foit_gxd = {-Three, Three};
Interval foit_gyd = {-Three, Three};

int foit_gn = 32;
