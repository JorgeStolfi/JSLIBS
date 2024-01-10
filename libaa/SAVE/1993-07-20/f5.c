/* Test function for 2d plots */

Float foit_f (Float x, Float y)
  {
    ROUND_NEAR;
    {
      Float x2 = x*x;
      Float y2 = y*y;
      Float r2 = x2 + y2;
      Float m2 = x2 - Three * y2;
      Float s = x*m2 + r2*r2;
      Float t = (1.0/8.0)/(2.0 + r2);
      Float res = s + t;
      return (res);
    }
  }

Interval foit_fv (Interval x, Interval y)
  {
    Float Eighth = Half*Quarter;
    Interval x2 = iv_sqr(x);
    Interval y2 = iv_sqr(y);
    Interval r2 = iv_add(x2, y2);
    Interval m2 = iv_sub(x2, iv_affine(y2, Three, Zero, Zero));
    Interval s = iv_add(iv_mul(x, m2), iv_sqr(r2));
    Interval t = iv_affine(iv_inv(iv_affine(r2, One, Two, Zero)), Eighth, Zero, Zero);
    Interval res = iv_add(s, t);
    return (res);
  }

FOIP foit_ff (FOIP x, FOIP y)
  {
    MemP frame = foi_top();
    Float Eighth = Half*Quarter;
    FOIP x2 = foi_sqr(x);
    FOIP y2 = foi_sqr(y);
    FOIP r2 = foi_add(x2, y2);
    FOIP m2 = foi_sub(x2, foi_affine(y2, Three, Zero, Zero));
    FOIP s = foi_add(foi_mul(x, m2), foi_sqr(r2));
    FOIP t = foi_affine(foi_inv(foi_affine(r2, One, Two, Zero)), Eighth, Zero, Zero);
    FOIP res = foi_add(s, t);
    return (foi_return(frame, res));
  }

char *foit_fname = 
  "r2 = x^2 + y^2;\nm2 = (x^2 - 3*y^2)^2;\ns = x*m2 - r2^2;\nt = (1/8)/(2 + r2);\nf = s + t";

Interval foit_fxd = {-Three/Two, Three/Two};
Interval foit_fyd = {-Three/Two, Three/Two};

int foit_fn = 32;
