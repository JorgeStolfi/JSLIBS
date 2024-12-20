/* Last edited on 2024-12-04 23:22:19 by stolfi */

r2_t tf_undistorted_sensor_coords_to_distorted_sensor_coords (camera_parameters_t cpar, r2_t pu)
{
  r2_t pd = tf_apply_kappa(&pu, cpar->kappa);
  return pd;
}

r2_t tf_distorted_sensor_coords_to_undistorted_sensor_coords (camera_parameters_t cpar, r2_t pd)
{
  r2_t pu = tf_apply_kappa(&pd, -cpar->kappa, &pd);
  return pu;
}

r2_t tf_apply_kappa(r2_t *p, double kappa)
{
  if (kappa == 0.0) {
    (*q) = (*p);
  } else { 
    /* Get the coordinates of {p}: */
    double pX = p->c[0];
    double pY = p->c[1];

    /* Compute distance squared {s} from optical axis: */
    double s = pX*pX + pY*pY;

    /* Compute the scaling factor {f}: */
    double m = 2*kappa*s;
    if (m >= 0.995) {
      /* Point is outside the maximum radius {R}: */
      m = 0.995;
    }
    double h = 1.0/(1.0 - m);
    double f = sqrt(h);

    /* Apply the radial scaling by {f}: */
    return (r2_t){{ pX*f, pY*f }};
  }
}

double tf_compute_maximum_safe_kappa
  ( int32_t Npx, 
    int32_t Npy, 
    double dpx,
    double dpy,
    double Cx,
    double Cy
  );
  /* Computes the maximum safe value of {kappa} from 
     the sensor size data. */
  
double tf_compute_maximum_safe_kappa
  ( int32_t Npx, 
    int32_t Npy, 
    double dpx,
    double dpy,
    double Cx,
    double Cy
  )
  {
    /* Find the max distance squared {R2} from any corner to optical axis: */
    double R2 = 0.0;
    int32_t sx,sy;
    for (sx = 0; sx <= 1; sx++) {
      for (sy = 0; sy <= 1; sy ++) {
        double duX = sx*Npx*dpx - Cx;
        double duY = sy*Npy*dpy - Cy;
        double du2 = duX*duX + duY*duY;
        if (du2 > R2) { R2 = du2; }
      }
    }
    
    /* Throw in a safety factor of 2.0 in the radius: */
    R2 = 2*2*R2;
    
    /* Compute the max {kappa} such that {2*kappa*R2 <= 0.995}: */
    double max_kappa = 0.995/R2/2;
    
    return max_kappa;
  }
