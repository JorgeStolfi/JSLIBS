/* Test function for 2d plots */

Float foit_f (Float x, Float y)
  {
    ROUND_NEAR;
    {
      Float ha = x*x + y*y + x*y - Quarter;
      Float fa = ha*ha*ha;
      Float hb = x*x + y*y - x*y - Quarter;
      Float fb = hb*hb*hb;
      Float res = (fa + fb - fa*fb);
      return (res);
    }
  }

Interval foit_fv (Interval x, Interval y)
  {
    Interval x2 = iv_sqr(x);
    Interval y2 = iv_sqr(y);
    Interval r2 = iv_add(x2, y2);
    Interval xy = iv_mul(x, y);
    Interval ha = iv_affine(iv_add(r2, xy), One, -Quarter, Zero);
    Interval fa = iv_mul(iv_sqr(ha), ha);
    Interval hb = iv_affine(iv_sub(r2, xy), One, -Quarter, Zero);
    Interval fb = iv_mul(iv_sqr(hb), hb);
    Interval fafb = iv_mul(fa, fb);
    Interval res = iv_sub (iv_add(fa, fb), fafb);
    return (res);
  }

FOIP foit_ff (FOIP x, FOIP y)
  {
    MemP frame = foi_top();
    FOIP x2 = foi_sqr(x);
    FOIP y2 = foi_sqr(y);
    FOIP r2 = foi_add(x2, y2);
    FOIP xy = foi_mul(x, y);
    FOIP ha = foi_affine(foi_add(r2, xy), One, -Quarter, Zero);
    FOIP fa = foi_mul(foi_sqr(ha), ha);
    FOIP hb = foi_affine(foi_sub(r2, xy), One, -Quarter, Zero);
    FOIP fb = foi_mul(foi_sqr(hb), hb);
    FOIP fafb = foi_mul(fa, fb);
    FOIP res = foi_sub (foi_add(fa, fb), fafb);
    return (foi_return(frame, res));
  }

char *foit_fname = 
  "r2 = x^2 + y^2;\nxy = x*y;\nfa = (r^2 + xy - 1/4)^3;\nfb = (r^2 - xy - 1/4)^3;\nf = fa + fb - fa*fb";

Interval foit_fxd = {-One, One};
Interval foit_fyd = {-One, One};

int foit_fn = 32;
