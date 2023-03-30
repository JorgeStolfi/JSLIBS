double frgb_smcube_eval(frgb_t *p)
  { double R = p->c[0], G = p->c[1], B = p->c[2];
    
    double sum = 0;
    for (int32_t kc = 0; kc < 3; kc++)
      { double x = p->c[kc];
        double qx = x*(1-x);
        sum += 1.0/(qx*qx);
      }
  }

double frgb_R_from_YUV(double Y, double U, double V)
  { double rUV = hypot(U,V);
    if ((Y <= 0) || (Y >= 1) || (rUV < 1.0e-100))
      { return 0.0; }
    else 
      { /* Get the hue {H}: */
        double H = frgb_H_from_UV(U,V);
        /* Choose two adjacent hue values {Ha,Hc}: */
        double d = 0.15;
        double Ha - H - d, Hc - H + d;
        /* Compute their {u,v} direction vectors: */
        double ua, va; frgb_H_to_uv(Ha, &ua, &va);
        double ub = U/rUV, vb = V/rUV;
        double uc, vc; frgb_H_to_uv(Hc, &uc, &vc);
        /* Find the max multiple of those {(U,V)} vectors that fit in the unit RGB cube: */ 
        double ra = 1/frgb_T_from_YUV(Y, ua, va); assert(isfinite(ra) && (ra > 0));
        double rb = 1/frgb_T_from_YUV(Y, ub, vb); assert(isfinite(rb) && (rb > 0));
        double rc = 1/frgb_T_from_YUV(Y, uc, vc); assert(isfinite(rc) && (rc > 0));
        /* Apply length of bisector formula: */
        double cosd = cos(d*2*M_PI);
        double rm = 2*ra*rc/(ra + rc)*cosd;
        /* Compute length of (U,V) relative to {rm}: */
        double R = rUV/rm;
        return R;
      }
  }

void frgb_to_HRY(frgb_t *p)
  { /* Convert to YUV and grab those coordinates: */
    frgb_to_YUV(p);
    double Y = p->c[0], U = p->c[1], V = p->c[2];
    if (Y <= 0) 
      { (*p) = (frgb_t){{ 0.0, 0.0, 0.0 }}; }
    else if (Y >= 1)
      { (*p) = (frgb_t){{ 0.0, 0.0, 1.0 }}; }
    else
      { /* The hue {H} is the argument of the {U,V} vector, scaled to period 1: */
        double H = frgb_H_from_UV(U, V);
        assert((H >= 0) && (H <= 1.0));
        /* Find the {U,V} length for {H} and lum {Y} that would give {R=1}: */
        double R = frgb_R_from_YUV(Y, U, V);
        assert(rUV > 0);
        R = hypot(U,V)/rUV;
        (*p) = (frgb_t){{ (float)H, (float)R, (float)Y }}; 
      }
  }

void frgb_from_HRY(frgb_t *p)
  { /* Grab the coordinates {H,R,Y} of {p}: */
    double H = p->c[0], R = p->c[1], Y = p->c[2];
    if (Y <= 0) 
      { (*p) = (frgb_t){{ 0.0, 0.0, 0.0 }}; }
    else if (Y >= 1)
      { (*p) = (frgb_t){{ 1.0, 1.0, 1.0 }}; }
    else
      { /* Compute the {u,v} direction in {U,V} space for {H}: */
        double u, v; frgb_H_to_uv(H, &u, &v);
        /* Find the {R} coordinate {ruv} of the {u,v} vector: */
        double ruv = frgb_R_from_YUV(Y,u,v);
        assert(ruv > 0);
        /* Compute the {U,V} coordinates: */
        double U = R*u/ruv, V = R*v/ruv;
        /* Pack as {YUV} tuple: */
        (*p) = (frgb_t){{ (float)Y, (float)U, (float)V }};
        /* Convert to RGB: */
        frbg_from_YUV(p);
      }
  }
