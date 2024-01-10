/* Test function for single-variable graph */

Float foit_g (Float x)
  {
    ROUND_NEAR;
    {
      Float yc = 27.0/16.0; /* Should be exact */
      Float y = yc;
      Float x2 = x*x;
      Float y2 = y*y;
      Float r2 = x2 + y2;
      Float m = x2 - Three * y2;
      Float m2 = m*m;
      Float res = x2*m2 - r2;
      return (res);
    }
  }

Interval foit_gv (Interval x)
  {
    ROUND_NEAR;
    {
      Float yc = 27.0/16.0; /* Should be exact */
      Interval y = {yc, yc};
      Interval x2 = iv_sqr(x);
      Interval y2 = iv_sqr(y);
      Interval r2 = iv_add(x2, y2);
      Interval m = iv_sub(x2, iv_affine(y2, Three, Zero, Zero));
      Interval m2 = iv_sqr(m);
      Interval res = iv_sub(iv_mul(x2, m2),  r2);
      return (res);
    }
  }

FOIP foit_gf (FOIP x)
  {
    MemP frame = foi_top();
    ROUND_NEAR;
    {
      Float yc = 27.0/16.0; /* Should be exact */
      FOIP y = foi_const(yc, Zero);
      FOIP x2 = foi_sqr(x);
      FOIP y2 = foi_sqr(y);
      FOIP r2 = foi_add(x2, y2);
      FOIP m = foi_sub(x2, foi_affine(y2, Three, Zero, Zero));
      FOIP m2 = foi_sqr(m);
      FOIP res = foi_sub(foi_mul(x2, m2),  r2);
      return (foi_return(frame, res));
    }
  }

char *foit_gname = "y = 27/16;\nr2 = x^2 + y^2;\nm2 = (x^2 - 3*y^2)^2;\ng = x^2*m2 - r2";

Interval foit_gxd = {-6.0, 6.0};
Interval foit_gyd = {-100.0, 100.0};

int foit_gn = 64;
