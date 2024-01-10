/* Test function for 2d plots */

Float foit_f (Float x, Float y)
  {
    ROUND_NEAR;
    {
      Float d = x - y;
      Float h = Quarter * d * (One - d*d);
      Float u = x + h;
      Float v = y + h;
      Float u2 = u*u;
      Float v2 = v*v;
      Float uv = u*v;
      Float m = u2 + v2 + uv + uv;
      Float res = m - Quarter;
      return (res);
    }
  }

Interval foit_fv (Interval x, Interval y)
  {
    Interval d = iv_sub(x, y);
    Interval h = iv_mul(d, iv_affine(iv_sqr(d), -Quarter, Quarter, Zero));
    Interval u = iv_add(x, h);
    Interval v = iv_add(y, h);
    Interval u2 = iv_sqr(u);
    Interval v2 = iv_sqr(v);
    Interval uv = iv_mul(u, v);
    Interval m = iv_add(iv_add(u2, v2), iv_add(uv, uv));
    Interval res = iv_affine(m, One, - Quarter, Zero);
    return (res);
  }

FOIP foit_ff (FOIP x, FOIP y)
  {
    MemP frame = foi_top();
    FOIP d = foi_sub(x, y);
    FOIP h = foi_mul(d, foi_affine(foi_sqr(d), -Quarter, Quarter, Zero));
    FOIP u = foi_add(x, h);
    FOIP v = foi_add(y, h);
    FOIP u2 = foi_sqr(u);
    FOIP v2 = foi_sqr(v);
    FOIP uv = foi_mul(u, v);
    FOIP m = foi_add(foi_add(u2, v2), foi_add(uv, uv));
    FOIP res = foi_affine(m, One, - Quarter, Zero);
    return (foi_return(frame, res));
  }

char *foit_fname = 
  "d = x-y; h = d(1-d^2); u = x+h; v = x+h; f = u^2 + v^2 + 2uv - 1/4";

Interval foit_fxd = {-Two, Two};
Interval foit_fyd = {-Two, Two};

int foit_fn = 32;
