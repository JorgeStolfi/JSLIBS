/* See voxb_obj.h */
/* Last edited on 2022-10-20 05:48:29 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <r3.h>
#include <affirm.h>

#include <voxb_obj.h>

bool_t voxb_obj_ball(r3_t *p, double R)
  {
    demand(R >= 0, "invalid radius");

    /* Get coordinates: */
    double X = p->c[0];
    double Y = p->c[1];
    double Z = p->c[2];
    
    /* Check distance from the object's nominal surface. */
    return X*X + Y*Y + Z*Z < R*R;
  }

bool_t voxb_obj_cube(r3_t *p, double R, double fillR)
  {
    demand((R >= fillR) && (fillR >= 0), "invalid fillet radius");

    /* Get coordinates (absolute): */
    double X = fabs(p->c[0]); if (X >= R) { return FALSE; }
    double Y = fabs(p->c[1]); if (Y >= R) { return FALSE; }
    double Z = fabs(p->c[2]); if (Z >= R) { return FALSE; }
    
    /* Sort the coordinates so that {X >= Y >= Z}: */
    if (X < Y) { double tmp = X; X = Y; Y = tmp; }
    if (X < Z) { double tmp = X; X = Z; Z = tmp; }
    if (Y < Z) { double tmp = Y; Y = Z; Z = tmp; }
    
    /* Check distance from the object's nominal surface. */
    double S = R - fillR; /* Min coordinate of fillet. */
    if (Y <= S)
      { /* Point {p'} is on the flat part of the face: */
        return TRUE;
      }
    else if (Z <= S)
      { /* Point {p'} is on a cylindrical edge fillet: */
        return (X-S)*(X-S) + (Y-S)*(Y-S) < fillR*fillR;
      }
    else
      { /* Point {p'} is on a spherical corner fillet: */
        return (X-S)*(X-S) + (Y-S)*(Y-S) + (Z-S)*(Z-S) < fillR*fillR;
      }
  }

bool_t voxb_obj_box(r3_t *p, double RX, double RY, double RZ, double fillR)
  {
    demand((fillR >= 0) && (RX >= fillR) && (RY >= fillR) && (RZ >= fillR), "invalid fillet radius");

    /* Get coordinates (absolute): */
    double X = fabs(p->c[0]); if (X >= RX) { return FALSE; }
    double Y = fabs(p->c[1]); if (Y >= RY) { return FALSE; }
    double Z = fabs(p->c[2]); if (Z >= RZ) { return FALSE; }
    
    /* Check distance from the object's nominal surface. */
    double SX = RX - fillR;
    double SY = RY - fillR;
    double SZ = RZ - fillR;
    if ((Y <= SY) && (Z <= SZ))
      { /* Point {p'} is on the flat part of the face parallel to {Y,Z}: */
        return TRUE;
      }
    else if ((X <= SX) && (Z <= SZ))
      { /* Point {p'} is on the flat part of the face parallel to {X,Z}: */
        return TRUE;
      }
    else if ((X <= SX) && (Y <= SY))
      { /* Point {p'} is on the flat part of the face parallel to {X,Y}: */
        return TRUE;
      }
    else if (Z <= SZ)
      { /* Point {p'} is on the fillet of the cylindrical edge parallel to {Z}: */
        return (X-SX)*(X-SX) + (Y-SY)*(Y-SY) < fillR*fillR;
      }
    else if (Y <= SY)
      { /* Point {p'} is on the fillet of the cylindrical edge parallel to {Y}: */
        return (X-SX)*(X-SX) + (Z-SZ)*(Z-SZ) < fillR*fillR;
      }
    else if (X <= SX)
      { /* Point {p'} is on the fillet of the cylindrical edge parallel to {X}: */
        return (Y-SY)*(Y-SY) +(Z-SZ)*(Z-SZ) < fillR*fillR;
      }
    else
      { /* Point {p'} is on a spherical corner fillet: */
        return (X-SX)*(X-SX) + (Y-SY)*(Y-SY) + (Z-SZ)*(Z-SZ) < fillR*fillR;
      }
  }

bool_t voxb_obj_rounded_box(r3_t *p, double RX, double RY, double RZ, double roundR, double fillR)
  {
    demand((fillR >= 0) && (RZ >= fillR), "invalid fillet radius");
    demand((roundR >= fillR) && (RX >= roundR) && (RY >= roundR), "invalid roundoff radius");

    /* Get coordinates (absolute), reject easy cases: */
    double X = fabs(p->c[0]); if (X >= RX) { return FALSE; }
    double Y = fabs(p->c[1]); if (Y >= RY) { return FALSE; }
    double Z = fabs(p->c[2]); if (Z >= RZ) { return FALSE; }
    
    /* Check distance from the object's nominal surface. */
    double TX = RX - roundR;
    double TY = RY - roundR;
    
    double SX = RX - fillR;
    double SY = RY - fillR;
    double SZ = RZ - fillR;
    if (Z <= SZ)
      { /* Point {p'} is on the flat or cylindrical side: */
        if (Y <= TY)
          { /* Point {p'} is on the flat part of the face parallel to {Y,Z}: */
            return TRUE;
          }
        else if (X <= TX)
          { /* Point {p'} is on the flat part of the face parallel to {X,Z}: */
            return TRUE;
          }
        else 
          { /* Point {p'} is on the cylindrical side: */
            return (X-TX)*(X-TX) + (Y-TY)*(Y-TY) < roundR*roundR;
          }
      }
    else
      { if (X <= TX)
          { if (Y <= SY)
              { /* Point {p'} is on the flat part of the face parallel to {X,Y}: */
                return TRUE;
              }
            else
              { /* Point {p'} is on the fillet of the top edge parallel to {X}: */
                 return (Y-SY)*(Y-SY) +(Z-SZ)*(Z-SZ) < fillR*fillR;
              }
          }
        else if (Y <= TY)
          { if (X <= SX)
              { /* Point {p'} is on the flat part of the face parallel to {X,Y}: */
                return TRUE;
              }
            else
              { /* Point {p'} is on the fillet of the top edge parallel to {Y}: */
                return (Y-SY)*(Y-SY) +(Z-SZ)*(Z-SZ) < fillR*fillR;
              }
          }
        else 
          { /* Point {p'} is on the top round corner or its fillet: */
            double rxy = hypot(X - TX, Y - TY) - (roundR - fillR);
            if (rxy <= 0)
              { /* Point {p'} is on the flat part of the round corner of the face parallel to {X,Y}: */
                return TRUE;
              }
            else
              { /* Point {p'} is on the fillet of the round corner of the face parallel to {X,Y}: */
                return rxy*rxy + (Z-SZ)*(Z-SZ) < fillR*fillR;
              }
          }
      }
  }

bool_t voxb_obj_donut(r3_t *p, double minR, double majR, int32_t ax)
  {
    demand((majR >= minR) && (minR >= 0), "invalid radii");
    demand((ax >= 0) && (ax < 3), "invalid axis {ax}");

    /* Get {Z} and quick check: */
    double Z = fabs(p->c[ax]);
    if (Z >= minR) { return FALSE; }
    
    /* Get {X,Y} coordinates change to cylindrical, quick check: */
    double X = p->c[(ax+1)%3];
    double Y = p->c[(ax+2)%3];
    double XY = hypot(X, Y);
    if ((XY <= majR - minR) || (XY >= majR + minR)) { return FALSE; }

    /* Check distance from the object's nominal surface. */
    return (XY - majR)*(XY-majR) + Z*Z < minR;
  }

bool_t voxb_obj_rod(r3_t *p, double H, double R, double fillR)
  {
    demand((R >= fillR) && (H >= fillR) && (fillR >= 0), "invalid fillet radius");
    
    /* Get {Z} in absolute value and quick check: */
    double Z = fabs(p->c[2]);
    if (Z >= H) { return FALSE; }

    /* Get {X,Y} coordinates change to cylindrical, quick check: */
    double X = p->c[0];
    double Y = p->c[1];
    double XY = hypot(X, Y);
    if (XY >= R) { return FALSE; }

    /* Check distance from the object's nominal surface. */
    if (Z < H - fillR) 
      { /* Main cylindrical part: */
        return TRUE;
      }
    else if (XY < R - fillR)
      { /* Flat part at top: */
        return TRUE;
      }
    else
      { /* Toroidal fillet: */
        double midR = R - fillR; /* Radius of fillet midline. */
        double midZ = H - fillR; /* {Z}-coordinate of fillet midline. */
        return (XY-midR)*(XY-midR) + (Z-midZ)*(Z-midZ) < fillR*fillR;
      }
  }
  
bool_t voxb_obj_tube(r3_t *p, double H, double Ri, double Ro, double fillR)
  {
    demand((Ro >= Ri) && (Ri > 0), "invalid radii");
    double thk = Ro - Ri;  /* Wall thickness. */
    demand((2*fillR <= thk) && (fillR <= H), "invalid fillet radius");
    
    /* Get {Z}  and quick check: */
    double Z = p->c[2];
    double aZ = fabs(Z);
    if (aZ >= H) { return FALSE; }

    /* Get {X,Y} coordinates change to cylindrical, quick check: */
    double X = p->c[0];
    double Y = p->c[1];
    double XY = hypot(X, Y);
    if ((XY >= Ro) || (XY <= Ri)) { return FALSE; }
    if ((aZ <= H - fillR) && (XY >= Ri) && (XY <= Ro)) { return TRUE; }

    /* Parameters of the fillets on top edge of tube: */
    double fZ = H - fillR;   /* {Z}-coordinate of axial line of fillet. */
    double fRi = Ri + fillR;   /* Radius of axial line of inner fillet. */
    double fRo = Ro - fillR;   /* Radius of axial line of outer fillet. */

    if (aZ > fZ) 
      { /* Nearest point is on top base or upper fillets: */
        if (XY < fRi)
          { /* Nearest point is in inner fillet: */
            return (XY-fRi)*(XY-fRi) + (aZ-fZ)*(aZ-fZ) < fillR*fillR;
	  }
	else if (XY > fRo)
	  { /* Nearest point is in outer fillet: */
            return (XY-fRo)*(XY-fRo) + (aZ-fZ)*(aZ-fZ) < fillR*fillR;
	  }
	else
	  { /* Nearest point is in flat part of base: */
            return TRUE;
	  }
      }
    else 
      { /* Nearest point is on one of the walls: */
        return (XY > Ri) && (XY < Ro);
      }
  }

bool_t voxb_obj_round_cup(r3_t *p, double H, double R, double thk, double fillR)
  {
    demand((R >= thk) && (thk >= 0), "invalid thickness");
    demand((fillR <= R) && (fillR >= thk) && (fillR + thk/2 <= 2*H), "invalid base fillet radius");
    
    /* Get {Z}  and quick check: */
    double Z = p->c[2];
    if (fabs(Z) >= H) { return FALSE; }

    /* Get {X,Y} coordinates change to cylindrical, quick check: */
    double X = p->c[0];
    double Y = p->c[1];
    double XY = hypot(X, Y);
    if (XY >= R) { return FALSE; }

    /* The upper edge of the cup has no flat part, just a half-torus. */
    /* Parameters of the full fillet on top edge of cup: */
    double tRmin = thk/2;    /* Minor radius of top fillet (half of wall thickness). */
    double tZ = H - tRmin;   /* {Z}-coordinate of axial line of top  edge fillet. */

    /* Parameters of the fillet around the base: */
    double bZ = - H + fillR; /* {Z} coordinate of axial line of base fillet. */
    double bRmaj = R - fillR;   /* Radius of axial line (major radius) of base fillet. */

    /* Check distance from the object's nominal surface. */
    if (Z > tZ) 
      { /* Point may be in upper lip fillet: */
        double tRmaj = R - tRmin;   /* Radius of axial line (major radius) of fillet on top edge.  */
        return (XY-tRmaj)*(XY-tRmaj) + (Z-tZ)*(Z-tZ) < tRmin*tRmin;
      }
    else if (Z >= bZ) 
      { /* Point may be inside wall: */
        double Ri = R - thk;
        return (XY > Ri);
      }
    else if (XY <= bRmaj)
      { /* Point may be inside bottom wall: */
        return (Z > -H) && (Z < -H + thk);
      }
    else
      { /* Point may be inside base fillet: */
        double d = hypot(XY - bRmaj, Z - bZ);
        return (d < fillR) && (d > fillR - thk);
      }
  }
