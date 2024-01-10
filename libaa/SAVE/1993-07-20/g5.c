/* Test function for 2d plots */

Float foit_g (Float x)
  {
    ROUND_NEAR;
    {
      Float u = Three*Quarter * (x - One);
      Float v = Three*Half * x;
      Float u2 = u*u;
      Float v2 = v*v;
      Float r2 = u2 + v2;
      Float m2 = u2 - Three * v2;
      Float s = u*m2 + r2*r2;
      Float t = (1.0/8.0)/(2.0 + r2);
      Float res = s + t;
      return (res);
    }
  }

Interval foit_gv (Interval x)
  {
    Float Eighth = Half*Quarter;
    Interval u = iv_affine(x, Three*Quarter, -Three*Quarter, Zero);
    Interval v = iv_affine(x, Three*Half, Zero, Zero);
    Interval u2 = iv_sqr(u);
    Interval v2 = iv_sqr(v);
    Interval r2 = iv_add(u2, v2);
    Interval m2 = iv_sub(u2, iv_affine(v2, Three, Zero, Zero));
    Interval s = iv_add(iv_mul(u, m2), iv_sqr(r2));
    Interval t = iv_affine(iv_inv(iv_affine(r2, One, Two, Zero)), Eighth, Zero, Zero);
    Interval res = iv_add(s, t);
    return (res);
  }

FOIP foit_gf (FOIP x)
  {
    MemP frame = foi_top();
    Float Eighth = Half*Quarter;
    FOIP u = foi_affine(x, Three*Quarter, -Three*Quarter, Zero);
    FOIP v = foi_affine(x, Three*Half, Zero, Zero);
    FOIP u2 = foi_sqr(u);
    FOIP v2 = foi_sqr(v);
    FOIP r2 = foi_add(u2, v2);
    FOIP m2 = foi_sub(u2, foi_affine(v2, Three, Zero, Zero));
    FOIP s = foi_add(foi_mul(u, m2), foi_sqr(r2));
    FOIP t = foi_affine(foi_inv(foi_affine(r2, One, Two, Zero)), Eighth, Zero, Zero);
    FOIP res = foi_add(s, t);
    return (foi_return(frame, res));
  }

char *foit_gname = 
  "u = 3(x-1)/4; v = 3x/2; r2 = u^2 + v^2;\nm2 = (u^2 - 3*v^2)^2;\ns = u*m2 - r2^2;\nt = (1/8)/(2 + r2);\ng = s + t";

Interval foit_gxd = {-One, One};
Interval foit_gyd = {-One, One};

int foit_gn = 1;
