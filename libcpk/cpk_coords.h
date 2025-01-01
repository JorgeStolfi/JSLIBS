#ifndef cpk_coords_H
#define cpk_coords_H

/* Conversion of latitude-longitude coords to/from plane coords. */
/* Last edited on 2024-12-31 10:42:18 by stolfi */

/* These procedures use the UTM projection, except that the reference
  meridian {refLon} is selected independently for each dataset,
  usually some sort of mean longitude. See USGS Bulletin 1532. The
  numerical model for the shape of the Earth is the South American
  1969 ellipsoid.
  
  Latitude and longitude are expressed in signed degrees, respectively North
  of the equator (in the range -90 to +90) and East of the Greenwhich meridian
  (in the range (-180 to +180).
  
  Planar coordinates are called Easthing {dx} and Northing {dy}, both in meters.
  Easting is measured orthogonally to the reference meridian, East from it.
  Northing is measured along the reference meridian, North from the equator. */

void LL_to_EUTM(double Lat, double Lon, double *dx, double *dy, double refLon);
  /* Converts the latitude-longitude pair {Lat,Lon} to EUTM coordinates {dx,dy}
   in meters. */

void EUTM_to_LL(double dx, double dy, double *Lat, double *Lon,  double refLon);
  /* Converts EUTM coordinates {dx,dy} in meters to latitude-longitude
   pair {Lat,Lon}. */

#endif
