/* See voxm_obj.h */
/* Last edited on 2022-10-20 05:47:06 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <r3.h>
#include <affirm.h>

#include <voxm_obj.h>

double voxm_obj_ball(r3_t *p, double R, double fuzzR)
  {
    demand(R >= 0, "invalid radius");

    /* Get coordinates: */
    double X = p->c[0]; if (fabs(X) >= R + fuzzR) { return 0.0; }
    double Y = p->c[1]; if (fabs(Y) >= R + fuzzR) { return 0.0; }
    double Z = p->c[2]; if (fabs(Z) >= R + fuzzR) { return 0.0; }
    
    /* Compute the distance {d} from the object's nominal surface. */
    double d = sqrt(X*X + Y*Y + Z*Z) - R;

    /* Map distance to indicator: */
    return voxm_obj_classify(d, fuzzR);
  }

double voxm_obj_cube(r3_t *p, double R, double fillR, double fuzzR)
  {
    demand((R >= fillR) && (fillR >= 0), "invalid fillet radius");

    /* Get coordinates (absolute), reject easy cases: */
    double X = fabs(p->c[0]); if (X >= R + fuzzR) { return 0.0; }
    double Y = fabs(p->c[1]); if (Y >= R + fuzzR) { return 0.0; }
    double Z = fabs(p->c[2]); if (Z >= R + fuzzR) { return 0.0; }
    
    /* Sort the coordinates so that {X >= Y >= Z}: */
    if (X < Y) { double tmp = X; X = Y; Y = tmp; }
    if (X < Z) { double tmp = X; X = Z; Z = tmp; }
    if (Y < Z) { double tmp = Y; Y = Z; Z = tmp; }
    
    /* Compute the distance {d} from the object's nominal surface. */
    double d;
    double S = R - fillR; /* Min coordinate of fillet. */
    if (Y <= S)
      { /* Point {p'} is on the flat part of the face: */
        d = X - R;
      }
    else if (Z <= S)
      { /* Point {p'} is on a cylindrical edge fillet: */
        d = hypot(X - S, Y - S) - fillR;
      }
    else
      { /* Point {p'} is on a cylindrical edge fillet: */
        d = hypot(hypot(X - S, Y - S), Z - S) - fillR;
      }

    /* Map distance to indicator: */
    return voxm_obj_classify(d, fuzzR);
  }

double voxm_obj_box(r3_t *p, double RX, double RY, double RZ, double fillR, double fuzzR)
  {
    demand((fillR >= 0) && (RX >= fillR) && (RY >= fillR) && (RZ >= fillR), "invalid fillet radius");

    /* Get coordinates (absolute), reject easy cases: */
    double X = fabs(p->c[0]); if (X >= RX + fuzzR) { return 0.0; }
    double Y = fabs(p->c[1]); if (Y >= RY + fuzzR) { return 0.0; }
    double Z = fabs(p->c[2]); if (Z >= RZ + fuzzR) { return 0.0; }
    
    /* Compute the distance {d} from the object's nominal surface. */
    double d;
    double SX = RX - fillR;
    double SY = RY - fillR;
    double SZ = RZ - fillR;
    if ((Y <= SY) && (Z <= SZ))
      { /* Point {p'} is on the flat part of the face parallel to {Y,Z}: */
        d = X - RX;
      }
    else if ((X <= SX) && (Z <= SZ))
      { /* Point {p'} is on the flat part of the face parallel to {X,Z}: */
        d = Y - RY;
      }
    else if ((X <= SX) && (Y <= SY))
      { /* Point {p'} is on the flat part of the face parallel to {X,Y}: */
        d = Z - RZ;
      }
    else if (Z <= SZ)
      { /* Point {p'} is on the fillet of the cylindrical edge parallel to {Z}: */
        d = hypot(X - SX, Y - SY) - fillR;
      }
    else if (Y <= SY)
      { /* Point {p'} is on the fillet of the cylindrical edge parallel to {Y}: */
        d = hypot(X - SX, Z - SZ) - fillR;
      }
    else if (X <= SX)
      { /* Point {p'} is on the fillet of the cylindrical edge parallel to {X}: */
        d = hypot(Y - SY, Z - SZ) - fillR;
      }
    else
      { /* Point {p'} is on a spherical corner fillet: */
        d = hypot(hypot(X - SX, Y - SY), Z - SZ) - fillR;
      }

    /* Map distance to indicator: */
    return voxm_obj_classify(d, fuzzR);
  }

double voxm_obj_rounded_box(r3_t *p, double RX, double RY, double RZ, double roundR, double fillR, double fuzzR)
  {
    demand((fillR >= 0) && (RZ >= fillR), "invalid fillet radius");
    demand((roundR >= fillR) && (RX >= roundR) && (RY >= roundR), "invalid roundoff radius");

    /* Get coordinates (absolute), reject easy cases: */
    double X = fabs(p->c[0]); if (X >= RX + fuzzR) { return 0.0; }
    double Y = fabs(p->c[1]); if (Y >= RY + fuzzR) { return 0.0; }
    double Z = fabs(p->c[2]); if (Z >= RZ + fuzzR) { return 0.0; }
    
    /* Compute the distance {d} from the object's nominal surface. */
    double d;
    double TX = RX - roundR;
    double TY = RY - roundR;
    
    double SX = RX - fillR;
    double SY = RY - fillR;
    double SZ = RZ - fillR;
    if (Z <= SZ)
      { /* Point {p'} is on the flat or cylindrical side: */
        if (Y <= TY)
          { /* Point {p'} is on the flat part of the face parallel to {Y,Z}: */
            d = X - RX;
          }
        else if (X <= TX)
          { /* Point {p'} is on the flat part of the face parallel to {X,Z}: */
            d = Y - RY;
          }
        else 
          { /* Point {p'} is on the cylindrical side: */
            d = hypot(X - TX, Y - TY) - roundR;
          }
      }
    else
      { if (X <= TX)
          { if (Y <= SY)
              { /* Point {p'} is on the flat part of the face parallel to {X,Y}: */
                d = Z - RZ;
              }
            else
              { /* Point {p'} is on the fillet of the top edge parallel to {X}: */
                d = hypot(Y - SY, Z  - SZ) - fillR;
              }
          }
        else if (Y <= TY)
          { if (X <= SX)
              { /* Point {p'} is on the flat part of the face parallel to {X,Y}: */
                d = Z - RZ;
              }
            else
              { /* Point {p'} is on the fillet of the top edge parallel to {Y}: */
                d = hypot(X - SX, Z  - SZ) - fillR;
              }
          }
        else 
          { /* Point {p'} is on the top round corner or its fillet: */
            double rxy = hypot(X - TX, Y - TY) - (roundR - fillR);
            if (rxy <= 0)
              { /* Point {p'} is on the flat part of the round corner of the face parallel to {X,Y}: */
                d = Z - RZ;
              }
            else
              { /* Point {p'} is on the fillet of the round corner of the face parallel to {X,Y}: */
                d = hypot(rxy, Z - SZ) - fillR;
              }
          }
      }

    /* Map distance to indicator: */
    return voxm_obj_classify(d, fuzzR);
  }

double voxm_obj_donut(r3_t *p, double minR, double majR, int32_t ax, double fuzzR)
  {
    demand((majR >= minR) && (minR >= 0), "invalid radii");
    demand((ax >= 0) && (ax < 3), "invalid axis {ax}");

    /* Get {Z} and quick check: */
    double Z = fabs(p->c[ax]);
    if (Z >= minR + fuzzR) { return 0.0; }
    
    /* Get {X,Y} coordinates change to cylindrical, quick check: */
    double X = p->c[(ax+1)%3];
    double Y = p->c[(ax+2)%3];
    double XY = hypot(X, Y);
    if ((XY <= majR - minR - fuzzR) || (XY >= majR + minR + fuzzR)) { return 0.0; }

    /* Compute the distance {d} from the object's nominal surface. */
    double d = hypot(XY - majR, Z) - minR; /* Distance to donut midline. */

    /* Map distance to indicator: */
    return voxm_obj_classify(d, fuzzR);
  }

double voxm_obj_rod(r3_t *p, double H, double R, double fillR, double fuzzR)
  {
    demand((R >= fillR) && (H >= fillR) && (fillR >= 0), "invalid fillet radius");
    
    /* Get {Z} in absolute value and quick check: */
    double Z = fabs(p->c[2]);
    if (Z >= H + fuzzR) { return 0.0; }

    /* Get {X,Y} coordinates change to cylindrical, quick check: */
    double X = p->c[0];
    double Y = p->c[1];
    double XY = hypot(X, Y);
    if (XY >= R + fuzzR) { return 0.0; }

    /* Compute the distance {d} from the object's nominal surface. */
    double d;
    if (Z < H - fillR) 
      { /* Main cylindrical part: */
        d = XY - R;
      }
    else if (XY < R - fillR)
      { /* Flat part at top: */
        d = Z - H;
      }
    else
      { /* Toroidal fillet: */
        double midR = R - fillR; /* Radius of fillet midline. */
        double midZ = H - fillR; /* {Z}-coordinate of fillet midline. */
        d = hypot(XY - midR, Z - midZ) - fillR;
      }

    /* Map distance to indicator: */
    return voxm_obj_classify(d, fuzzR);
  }
  
double voxm_obj_tube(r3_t *p, double H, double Ri, double Ro, double fillR, double fuzzR)
  {
    demand((Ro >= Ri) && (Ri > 0), "invalid radii");
    double thk = Ro - Ri;  /* Wall thickness. */
    demand((2*fillR <= thk) && (fillR <= H), "invalid fillet radius");
    
    /* Get {Z}  and quick check: */
    double Z = p->c[2];
    double aZ = fabs(Z);
    if (aZ >= H + fuzzR) { return 0.0; }

    /* Get {X,Y} coordinates change to cylindrical, quick check: */
    double X = p->c[0];
    double Y = p->c[1];
    double XY = hypot(X, Y);
    if ((XY >= Ro + fuzzR) || (XY <= Ri - fuzzR)) { return 0.0; }
    if ((aZ <= H - fmax(fuzzR,fillR)) && (XY >= Ri + fuzzR) && (XY <= Ro - fuzzR)) { return 1.0; }

    /* Parameters of the fillets on top edge of tube: */
    double lipaxZ = H - fillR;   /* {Z}-coordinate of axial line of fillet. */
    double lipaxRi = Ri + fillR;   /* Radius of axial line of inner fillet. */
    double lipaxRo = Ro - fillR;   /* Radius of axial line of outer fillet. */

    /* Compute the distance {d} from the object's nominal surface. */
    double d;
    if (aZ > lipaxZ) 
      { /* Nearest point is on top base or upper fillets: */
        if (XY < lipaxRi)
          { /* Nearest point is in inner fillet: */
            d = hypot(lipaxRi - XY, aZ - lipaxZ) - fillR;
	  }
	else if (XY > lipaxRo)
	  { /* Nearest point is in outer fillet: */
            d = hypot(XY - lipaxRo, aZ - lipaxZ) - fillR;
	  }
	else
	  { /* Nearest point is in flat part of base: */
            d = aZ - H;
	  }
      }
    else 
      { /* Nearest point is on one of the walls: */
        double Rm = (Ri + Ro)/2; /* radius of mid-wall surface. */
      if (XY <= Rm)
          { /* Nearest point is on inner cylindrical wall: */
            d = Ri - XY;
          }
        else 
          { /* Nearest point is on outer cylindrical wall: */
            d = XY - Ro;
          }
      }

    /* Map distance to indicator: */
    return voxm_obj_classify(d, fuzzR);
  }

double voxm_obj_round_cup
  ( r3_t *p, 
    double H,
    double R, 
    double thk,
    double fillR,
    double fuzzR
  )
  {
    demand((R >= thk) && (thk >= 0), "invalid thickness");
    demand((fillR <= R) && (fillR >= thk) && (fillR + thk/2 <= 2*H), "invalid base fillet radius");
    
    /* Get {Z}  and quick check: */
    double Z = p->c[2];
    if (fabs(Z) >= H + fuzzR) { return 0.0; }

    /* Get {X,Y} coordinates change to cylindrical, quick check: */
    double X = p->c[0];
    double Y = p->c[1];
    double XY = hypot(X, Y);
    if (XY >= R + fuzzR) { return 0.0; }

    double wallR = thk/2; /* Half of wall thickness. */

    /* Parameters of the full fillet on top edge of cup: */
    double lipaxZ = H - wallR;   /* {Z}-coordinate of axial line of fillet. */

    /* Parameters of the fillet around the base: */
    double bfilaxZ = - H + fillR; /* {Z} coordinate of axial line of base fillet. */
    double bfilaxR = R - fillR;   /* Radius of axial line (major radius) of base fillet. */

    /* Compute the distance {d} from the object's nominal surface. */
    double d;
    if (Z > lipaxZ) 
      { /* Nearest point is on upper fillet: */
        double lipaxR = R - wallR;   /* Radius of axial line (major radius) of fillet on top edge.  */
        double lipmiR = wallR;       /* Minor radius of fillet on top edge. */
        d = hypot(XY - lipaxR, Z - lipaxZ) - lipmiR;
      }
    else if (Z > bfilaxZ) 
      { /* Nearest point is on cylindrical wall: */
        double mwR = R - wallR;   /* Radius of cylindrical midsurf. */
        if (XY > mwR)
          { double otR = R; /* Radius of outer cylindrical surface. */
            d = XY - otR; }
        else
          { double inR = R - thk; /* Radius of inner cylindrical surface. */
            d = inR - XY;
          }
      }
    else if (XY < bfilaxR)
      { /* Nearest point is in the flat bottom wall: */
        double mwZ = - H + wallR; /* {Z} coordinate of bottom midsurf. */
        if (Z < mwZ)
          { double otZ = -H; /*{Z} of flat bottom's outer surface. */
            d = otZ - Z; 
          }
        else 
          { double inZ = - H + thk; /* {Z} of flat bottom's inner surface. */
            d = Z - inZ;
          }
      }
    else
      { /* Nearest point is on base fillet: */
        d = hypot(XY - bfilaxR, Z - bfilaxZ);
        double mwR = fillR - wallR; /* Minor radius of fillet's midsurf. */
        if (d > mwR)
          { double otR = fillR;    /* Minor radius of fillet's outer surface. */
            d = d - otR;
          }
        else
          { double inR = fillR - thk;  /* Minor radius of fillet's inner surface. */
            d = inR - d;
          }
      }

    /* Map distance to indicator: */
    return voxm_obj_classify(d, fuzzR);
  }

double voxm_obj_classify(double d, double fuzzR)
  {
    if (d <= -fuzzR) 
      { return 1.0; }
    else if (d >= +fuzzR)
      { return 0.0; }
    else
      { double s = (fuzzR - d)/(2*fuzzR);
        return s;
      }
  }
