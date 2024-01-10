/* Test function for 2d plots */

Float foit_g (Float x)
  {
    ROUND_NEAR;
    {
      Float u = Three*Quarter * (x - One);
      Float v = Three*Half * x;
      Float u2 = u*u;
      Float v2 = v*v;
      Float m2 = u2 - Three * v2;
      Float s = u*m2;
      Float res = s;
      return (res);
    }
  }

Interval foit_gv (Interval x)
  {
    Interval u = iv_affine(x, Three*Quarter, -Three*Quarter, Zero);
    Interval v = iv_affine(x, Three*Half, Zero, Zero);
    Interval u2 = iv_sqr(u);
    Interval v2 = iv_sqr(v);
    Interval m2 = iv_sub(u2, iv_affine(v2, Three, Zero, Zero));
    Interval s = iv_mul(u, m2);
    Interval res = s;
    return (res);
  }

FOIP foit_gf (FOIP x)
  {
    MemP frame = foi_top();
    FOIP u = foi_affine(x, Three*Quarter, -Three*Quarter, Zero);
    FOIP v = foi_affine(x, Three*Half, Zero, Zero);
    FOIP u2 = foi_sqr(u);
    FOIP v2 = foi_sqr(v);
    FOIP m2 = foi_sub(u2, foi_affine(v2, Three, Zero, Zero));
    FOIP s = foi_mul(u, m2);
    FOIP res = s;
    return (foi_return(frame, res));
  }

char *foit_gname = 
  "u = 3(x-1)/4; v = 3x/2; r2 = u^2 + v^2;\nm2 = (u^2 - 3*v^2)^2;\ng = u*m2";

Interval foit_gxd = {-One, One};
Interval foit_gyd = {-One, One};

int foit_gn = 1;
