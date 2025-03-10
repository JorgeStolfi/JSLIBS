/* see r2_extra.h */
/* Last edited on 2024-12-01 10:20:44 by stolfi */

#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <jsmath.h>
#include <jsqroots.h>
#include <r2.h>
#include <r2x2.h>
#include <r3x3.h>
#include <affirm.h>

#include <r2_extra.h>

void r2_map_projective(r2_t *p, r3x3_t *M, r2_t *q, r2x2_t *J)
  { 
    bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "!! > %s\n", __FUNCTION__); }
    if (debug) { fprintf(stderr, "!! M ="); r3x3_print(stderr, M); fputc('\n',stderr); }
    if (debug) { fprintf(stderr, "!! p = "); r2_print(stderr, p); fputc('\n',stderr); }
    
    /* If {p} is already invalid, return it: */
    if (isnan(p->c[0]) || isnan(p->c[1])) { (*q) = (*p); return; }
    
    /* Convert {p} to homogeneous coordinates {hp}: */
    r3_t hp = (r3_t) {{ 1.0, p->c[0], p->c[1] }}; 
    if (debug) { fprintf(stderr, "!! hp = "); r3_print(stderr, &hp); fputc('\n',stderr); }

    /* Map {hp} to {hq}, get the homogeneous coordinates {qw,qx,qy}: */
    r3_t hq;
    r3x3_map_row(&hp, M, &hq);
    if (debug) { fprintf(stderr, "!! hq = "); r3_print(stderr, &hq); fputc('\n',stderr); }
    double qw = hq.c[0]; 
    double qx = hq.c[1]; 
    double qy = hq.c[2]; 
    if (debug) { fprintf(stderr, "!! hq = [ %24.16e  %24.16e  %24.16e ]\n", qw,qx,qy); }
    
    /* Check for validity of result: */
    if (qw <= 0.0)
      { /* Result is at infinity, or beyond: */
        q->c[0] = q->c[1] = NAN;; 
        return; 
      }
    else
      { /* Compute the Cartesian coordinates {qX,qY} of {q}: */
        double qX = qx/qw; 
        double qY = qy/qw;
        if (debug) { fprintf(stderr, "!! q = ( %24.16e  %24.16e )\n", qX, qY); }
        
        /* Update the given point {p}: */
        q->c[0] = qX; q->c[1] = qY;
        
        if (J != NULL)
          { /* Compute the Jacobian of {hq} (homogeneous) w.r.t {p} (Cartesian): */
            double dqw_dpX = M->c[1][0], dqw_dpY = M->c[2][0];
            double dqx_dpX = M->c[1][1], dqx_dpY = M->c[2][1];
            double dqy_dpX = M->c[1][2], dqy_dpY = M->c[2][2];
            if (debug) { fprintf(stderr, "!! dhq/dpX = [ %24.16e  %24.16e  %24.16e ]\n", dqw_dpX, dqx_dpX, dqy_dpX); }
            if (debug) { fprintf(stderr, "!! dhq/dpY = [ %24.16e  %24.16e  %24.16e ]\n", dqw_dpY, dqx_dpY, dqy_dpY); }

            /* Compute the Jacobian {K} of {q} (Cartesian) w.r.t {p} (Cartesian): */
            r2x2_t K;
            K.c[0][0]= (dqx_dpX - qX*dqw_dpX)/qw; /* dqX_dpX */ 
            K.c[0][1]= (dqy_dpX - qY*dqw_dpX)/qw; /* dqY_dpX */ 
            K.c[1][0]= (dqx_dpY - qX*dqw_dpY)/qw; /* dqX_dpY */ 
            K.c[1][1]= (dqy_dpY - qY*dqw_dpY)/qw; /* dqY_dpY */ 

            /* Update the Jacobian {J}: */
            r2x2_mul(J, &K, J);
          }
      }
    if (debug) { fprintf(stderr, "!! < %s\n", __FUNCTION__); }
  }

void r2_map_radial(r2_t *p, r2_t *h, double kappa, r2x2_t *J)
  {
    /* If {p} is already invalid, do nothing: */
    if (isnan(p->c[0]) || isnan(p->c[1])) { return; }
    
    /* If the radial map is not the identity, apply it: */
    if (kappa != 0)
      { /* Get the pixel sensor dimensions {hX,hY} in mm: */
        double hX = h->c[0];
        double hY = h->c[1];
      
        /* Get point {pX,pY} (projected coordinates, in mm): */
        double pX = p->c[0] * hX;
        double pY = p->c[1] * hY;
        
        /* Compute distance squared {s} from optical axis: */
        double s = pX*pX + pY*pY;
        
        /* Compute the scaling factor {f} and its derivative {df_ds}: */
        double m = 2*kappa*s;
        if (m >= 0.995) 
          { /* Point is outside the maximum radius {R}: */
            p->c[0] = p->c[1] = NAN;; 
            return; 
          }
        double h = 1.0/(1.0 - m);
        double f = sqrt(h);
        double df_ds = kappa*h*h/f;
        
        /* Compute the mapped point coordinates: */
        double qX = pX*f; 
        double qY = pY*f;
        
        /* Update the given point {p} (image coords, in pixels): */
        p->c[0] = qX / hX; p->c[1] = qY / hY;

        if (J != NULL)
          { /* Compute the Jacobian {K} of this map: */
            r2x2_t K;
            K.c[0][0]= f + pX*df_ds*2*pX;           /* dqX_dpX */ 
            K.c[0][1]= 0 + pY*df_ds*2*pX * hX / hY; /* dqY_dpX */ 
            K.c[1][0]= 0 + pX*df_ds*2*pY * hY / hX; /* dqX_dpY */ 
            K.c[1][1]= f + pY*df_ds*2*pY;           /* dqY_dpY */ 

            /* Update the Jacobian {J}: */
            r2x2_mul(J, &K, J);
          }
      }
  }


void r2_get_persp_rectangle_bbox
  ( interval_t tbox[],   /* Rectangle in true coordinates. */
    r3x3_t *T2I,         /* True-to-image projective map matrix. */
    interval_t ibox[]    /* (OUT) bonding box in image coordinates. */
  )
  {
    /* Initialize bounds: */
    for (uint32_t ax = 0;  ax < 2; ax++) { ibox[ax].end[0] = +INF; ibox[ax].end[1] = -INF; }
    /* Hack: map the corners of a slightly wider rectangle and get its bbox */
    double wx = HI(tbox[0]) - LO(tbox[0]);
    double wy = HI(tbox[1]) - LO(tbox[1]);
    for (uint32_t dx = 0;  dx <= 1; dx++)
      { for (uint32_t dy = 0;  dy <= 1; dy++)
          { /* Get a true corner {p} of the enlarged rectangle: */
            double txp = tbox[0].end[dx] + (2*dx - 1)*0.00001*wx;
            double typ = tbox[1].end[dy] + (2*dy - 1)*0.00001*wy;
            r2_t p = (r2_t){{ txp, typ }};
            /* Map {p} to the image coord system: */
            r2_map_projective(&p, T2I, &p, NULL);
            /* Expand box: */
            for (uint32_t ax = 0;  ax < 2; ax++)
              { if (p.c[ax] < ibox[ax].end[0]) { ibox[ax].end[0] = p.c[ax]; }
                if (p.c[ax] > ibox[ax].end[1]) { ibox[ax].end[1] = p.c[ax]; }
              }
	  }
      }
  }

void r2_get_persp_disk_bbox
  ( r2_t *ctr,             /* Disk center in true coordinates. */
    double rad,            /* Disk radius in true coordinates. */
    r3x3_t *T2I,           /* True-to-image projective map matrix. */
    interval_t ibox[]      /* (OUT) bonding box. */
  )
  {
    /* Initialize bounds: */
    for (uint32_t ax = 0;  ax < 2; ax++) { ibox[ax].end[0] = +INF; ibox[ax].end[1] = -INF; }
    /* Hack: map the corners of an enclosing 8-gon and get its bbox */
    double a = 1.00001*rad;             /* Max abs coord of an octagon corner. */
    double b = 1.00001*rad*tan(M_PI/8); /* Min abs coord of an octagon corner. */
    for (int32_t dx = -1; dx <= +1; dx += 2)
      { for (int32_t dy = -1; dy <= +1; dy += 2)
          { /* Get a corner {txp,typ} of the octagon: */
            double txp = ctr->c[0] + dx*a;
            double typ = ctr->c[1] + dy*b;
            /* Try the two transposed copies of the point: */
            for (uint32_t swap = 0;  swap < 2; swap++)
	      { /* Map it to the image coord system: */
                r2_t p = (r2_t){{ txp, typ }};
                r2_map_projective(&p, T2I, &p, NULL);
                /* Expand box: */
                for (uint32_t ax = 0;  ax < 2; ax++)
                  { if (p.c[ax] < ibox[ax].end[0]) { ibox[ax].end[0] = p.c[ax]; }
		    if (p.c[ax] > ibox[ax].end[1]) { ibox[ax].end[1] = p.c[ax]; }
		  }
                /* Swap true coords {txp,typ}: */
                double zz = txp; txp = typ; typ = zz;
	      }
	  }
      }
  }

bool_t r2_pixel_is_inside_persp_rectangle
  ( int32_t x,
    int32_t y,
    double mrg, 
    r3x3_t *I2T,
    interval_t tbox[]   /* Rectangle in true coordinates. */
  )
  {
    /* Check the four corners of the pixel: */
    for (int32_t dx = -1; dx <= +1; dx += 2)
      { for (int32_t dy = -1; dy <= +1; dy += 2)
          { /* Get a corner {p} of the pixel, expanded by {mrg}: */
            double ixp = x + 0.5 + dx*(0.5 + mrg);
            double iyp = y + 0.5 + dy*(0.5 + mrg);
	    r2_t p = (r2_t){{ ixp, iyp }};
            /* Map it to the true coord system: */
            r2_map_projective(&p, I2T, &p, NULL);
            /* Check if it is inside the rectangle: */
            if ((p.c[0] < LO(tbox[0])) || (p.c[0] > HI(tbox[0]))) { return FALSE; }
            if ((p.c[1] < LO(tbox[1])) || (p.c[1] > HI(tbox[1]))) { return FALSE; }
	  }
      }
    return TRUE;
  }

bool_t r2_pixel_is_inside_persp_disk
  ( int32_t x,
    int32_t y,
    double mrg, 
    r3x3_t *I2T,
    r2_t *ctr, 
    double rad
  )
  {
    /* Check the four corners of the pixel: */
    double r2 = rad*rad;
    for (int32_t dx = -1; dx <= +1; dx += 2)
      { for (int32_t dy = -1; dy <= +1; dy += 2)
          { /* Get a corner {p} of the pixel, expanded by {mrg}: */
            double ixp = x + 0.5 + dx*(0.5 + mrg);
            double iyp = y + 0.5 + dy*(0.5 + mrg);
	    r2_t p = (r2_t){{ ixp, iyp }};
            /* Map it to the true coord system: */
            r2_map_projective(&p, I2T, &p, NULL);
            /* Check if it is inside the disk: */
            double d2 = r2_dist_sqr(&p, ctr);
            if (d2 > r2) { return FALSE; }
	  }
      }
    return TRUE;
  }

void r2_map_twirl(r2_t *p,  r2_t *ctr, double rad, double ang, r2x2_t *J)
  { 
    r2_t u;
    r2_sub(p, ctr, &u);
    double R2 = rad*rad;
    double r2 = r2_norm_sqr(&u);
    double z = r2/R2;
    double w = 1/(1 + z*(1 + 0.5*z)); /* Intensity of whirl at {q}. */

    double a = -ang*w;      /* Actual rotation angle. */
    double ca = cos(a);           
    double sa = sin(a);           
    r2_t v = (r2_t){{ u.c[0]*ca - u.c[1]*sa, u.c[0]*sa + u.c[1]*ca }};
    r2_add(ctr, &v, p);
    if (J != NULL)
      { 
        double dvx_dux = +ca;
        double dvy_dux = +sa;
        double dvx_duy = -sa;
        double dvy_duy = +ca;

        double dvx_dw = -ang*(- u.c[0]*sa - u.c[1]*ca);
        double dvy_dw = -ang*(+ u.c[0]*ca - u.c[1]*sa);

        double dw_dz = -w*w*(1 + z);
        double dw_dr2 = dw_dz/R2;

        double dr2_dux = 2*u.c[0];
        double dr2_duy = 2*u.c[1];

        dvx_dux += dvx_dw*dw_dr2*dr2_dux;
        dvx_duy += dvx_dw*dw_dr2*dr2_duy;
        dvy_dux += dvy_dw*dw_dr2*dr2_dux;
        dvy_duy += dvy_dw*dw_dr2*dr2_duy;

        r2x2_t K;
        K.c[0][0] = dvx_dux;
        K.c[0][1] = dvy_dux;
        K.c[1][0] = dvx_duy;
        K.c[1][1] = dvy_duy;

        r2x2_mul(J, &K, J);
      }
  }

void r2_map_expand(r2_t *p, double xlo, double xhi, double ylo, double yhi, r2x2_t *J)
  {
    double dx, dy;
    expand_range(&(p->c[0]), xlo, xhi, &dx);
    expand_range(&(p->c[1]), ylo, yhi, &dy);
    if ((J != NULL) && isfinite(p->c[0]) && isfinite(p->c[1]))
      { r2x2_t K = (r2x2_t){{{ dx, 0 },{0, dy}}};
        r2x2_mul(J, &K, J);
      }
  }
                   
void r2_map_contract(r2_t *p, double xlo, double xhi, double ylo, double yhi, r2x2_t *J)
  {
    double dx, dy;
    contract_range(&(p->c[0]), xlo, xhi, &dx);
    contract_range(&(p->c[1]), ylo, yhi, &dy);
    if ((J != NULL) && (isfinite(p->c[0])) && (isfinite(p->c[1])))
      { r2x2_t K = (r2x2_t){{{ dx, 0 },{0, dy}}};
        r2x2_mul(J, &K, J);
      }
  }

void r2_clip_seg_to_unit_disk(r2_t *a, r2_t *b, double *ta, double *tb)
  { /* Set up the parametric equation {a + t*u} of the segment: */
    double ax = a->c[0];
    double ay = a->c[1];
    double ux = b->c[0] - ax;
    double uy = b->c[1] - ay;
    /* Insert in circle equation {x^2 + y^2 - 1}: */
    double A = ux*ux + uy*uy;
    double B = 2*(ux*ax + uy*ay);
    double C = ax*ax + ay*ay - 1.0;
    if (A == 0.0)
      { /* Points {a} and {b} coincide. Only check if inside disk. */
        if (C >= 0.0)
          { /* Inside or on border: */
            (*ta) = 0.0; (*tb) = 1.0;
          }
        else
          { /* Outside: */
            (*ta) = 0.9; (*tb) = 0.1;
          }
      }
    else
      { /* Points are distinct.  There must be two roots. */
        /* Find roots of equation: */
        double r1, r2, im; /* Real and imaginary parts. */
        int32_t sd; /* Sign of discriminant. */
        sd = roots_quadratic(A, B, C, &r1, &r2, &im);
        if (sd < 0)
          { /* Line {a--b} does not touch or cross disk, return empty interval: */
            assert(r1 == r2);
            assert(im > 0.0);
            (*ta) = 0.9; (*tb) = 0.1; 
          }
        else 
          { /* Line {a--b} touches or crosses disk: */
            assert(im == 0);
            if (sd == 0) { assert(r1 == r2); } else { assert(r1 < r2); }
            /* Parameter interval in closed disk is {[r1 _ r2]}. */
            if ((r2 < 0.0) || (r1 > 1.0))
              { /* Intersection is outside seg {a--b}, return empty interval: */
                (*ta) = 0.9; (*tb) = 0.1; 
              }
            else
              { (*ta) = (r1 < 0.0 ? 0.0 : r1);
                (*tb) = (r2 > 1.0 ? 1.0 : r2);
              }
          }
      }
  }

void r2_debug_point_jacobian(char *label, r2_t *p, r2x2_t *J, char *tail)
  {
    fprintf(stderr, "    %s", label);
    if (isnan(p->c[0]) || isnan(p->c[1]))
      { fprintf(stderr, " is invalid\n"); }
    else
      { fprintf(stderr, " = (%15.8f %15.8f)", p->c[0], p->c[1]);
        if (J != NULL)
          { fprintf
              ( stderr, 
                "  J = [[%6.1f %6.1f] [%6.1f %6.1f]]", 
                J->c[0][0], J->c[0][1], J->c[1][0], J->c[1][1]
              );
          }
      }
    fprintf(stderr, "%s", tail);
  }

void r2_map_compute_numeric_jacobian(r2_t *p, r2_map_jacobian_t *map, double step, r2x2_t *K, bool_t debug)
  {
    /* Check the partial derivatives numerically: */
    for (uint32_t i = 0; i < 2; i++)
      { /* Evaluate {map} at points {q[0..1]} , displaced from {p} along axis {i}: */
        r2_t q[2];
        for (uint32_t k = 0; k < 2; k++)
          { r2_t pk = (*p);
            pk.c[i] += (2*(int32_t)k - 1)*step;
            if (debug) { r2_debug_point_jacobian("    pk(1)", &pk, NULL, ""); }
            r2_t *qk = &(q[k]);
            (*qk) = pk; map(qk, NULL);
            if (debug) { r2_debug_point_jacobian("  qk", qk, NULL, "\n"); }
            if (isnan(qk->c[0]) || isnan(qk->c[1])) 
              { (*K) = (r2x2_t){{{ NAN, NAN}, {NAN, NAN}}};
                return;
              }
          }
        /* Compute the derivatives of : */
        for (uint32_t j = 0;  j < 2; j++)
          { /* Compute the derivative of coordinate {j} w.r.t coordinate {i}: */
            K->c[i][j] = (q[1].c[j] - q[0].c[j])/(2*step);
          }
      }
  }

void r2_map_check_jacobian(r2_t *p, r2_map_jacobian_t *map, char *mapname, double eps, bool_t debug)
  {
    if (debug) { fprintf(stderr, "    checking the Jacobian of %s ...\n", mapname); }
    if (debug) { fprintf(stderr, "!! p =  "); r2_print(stderr, p); fputc('\n', stderr); }
    /* Generate a `random' invertible matrix {M} to be the prior Jacobian: */
    r2x2_t M = (r2x2_t){{{ 3.14, 2.17 }, { -7.0, 1.41 }}};
    /* r2x2_t M = (r2x2_t){{{ 1, 0 }, { 0, 1 }}}; */
    /* Append to {M} the analytic Jacobian of {map}: */
    r2x2_t J = M;
    r2_t pt = (*p);
    if (debug) { fprintf(stderr, "!! J before "); r2x2_print(stderr, &J); fputc('\n', stderr); }
    map(&pt, &J);
    if (debug) { fprintf(stderr, "!! J after  "); r2x2_print(stderr, &J); fputc('\n', stderr); }
    /* Remove the prior Jacobian from {J}: */
    if (debug) { fprintf(stderr, "!! M dir    "); r2x2_print(stderr, &M); fputc('\n', stderr); }
    r2x2_inv(&M, &M);
    if (debug) { fprintf(stderr, "!! M inv    "); r2x2_print(stderr, &M); fputc('\n', stderr); }
    r2x2_mul(&M, &J, &J);
    if (debug) { fprintf(stderr, "!! J pure   "); r2x2_print(stderr, &J); fputc('\n', stderr); }
    /* Get the magnitude {mag_J} of the Jacobian: */
    double mag_J = fabs(J.c[0][0]) + fabs(J.c[0][1]) + fabs(J.c[1][0]) + fabs(J.c[1][1]); 
    if (mag_J < 1.0e-180) { return; }
    /* Get the magnitude {mag_p} of {p}: */
    double mag_p = fabs(p->c[0]) + fabs(p->c[1]);
    /* Choose a small enough step ... */
    double mag = mag_p/hypot(1.0, mag_J);
    double step = hypot(eps*mag, 1.0e-7);
    if (debug) { fprintf(stderr, "    step = %14.6e\n", step); }
    /* Compute the numeric Jacobian: */
    r2x2_t K;
    r2_map_compute_numeric_jacobian(p, map, step, &K, debug);
    /* Establish a limit for the relative error: */
    double tol = 1.0e-4;
    /* Check the partial derivatives numerically: */
    for (uint32_t i = 0;  i < 2; i++)
      { /* Check the derivatives of {ip} w.r.t. coordinate {i} of {op}: */
        for (uint32_t j = 0;  j < 2; j++)
          { /* Compute the derivative of coordinate {j} w.r.t coordinate {i}: */
            double dnum = K.c[i][j];
            /* Compare with the Jacobian returned by {map}: */
            double djac = J.c[i][j]; 
            double err = fabs(dnum-djac);
            if (err > tol*mag_J)
              { fprintf(stderr, "    ** Jacobian of %s is inconsistent (err = %+14.6e)\n", mapname, err); 
                fprintf(stderr, "    J[%d][%d] = %15.8e", i, j, djac);
                fprintf(stderr, "  numeric d_ip[%d]/d_op[%d] = %15.8e", j, i, dnum);
                fprintf(stderr, "\n");
                affirm(FALSE, "Aborted");
              }
          }
      }
  } 

