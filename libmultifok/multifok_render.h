/* Tools for computing apparent colors of surfaces. */
/* Last edited on 2025-02-06 22:14:33 by stolfi */

#ifndef multifok_render_H
#define multifok_render_H

#include <stdint.h>

#include <r3.h>
#include <frgb.h>

frgb_t multifok_render_compute_color
  ( r3_t *uVis,
    r3_t *sNrm,
    frgb_t *sGlo,
    frgb_t *sLam,
    r3_t *uLit,
    frgb_t *cLit,
    frgb_t *cIso,
    bool_t debug
  );
  /* Computes the color {clr} of an infintesimal smooth surface patch
    with a mix of glossy and Lambertian finishes, under a mix of
    isotropic ("ambient") and unidirectional illumination fields.
    
    Let {p} be the center of the surface patch. The parameters are
    
      * {uvis} a unit vector pointing from the observer towards {p}.
      
      * {sNrm} the outwards surface normal at {p}, a unit vector.
      
      * {sGlo} the fraction of BRDF at {p} that is glossy. 
      
      * {sLam} the albedo of the Lambertian component of the BRDF at {p}. 
      
      * {uLit} the uniy directin vector from {p} towards the light source.
      
      * {cLit} the intensity at {p} of the light from the point source.
      
      * {cIso} the intensity of the isotropic light field at {p}. 
      
    The quantities {sGlo,sLam,cIso,cLit} are specific to each channel {c}.
    
    If {uLit} is not {NULL},  should be a unit vector
    pointing towards the simulated point light source. For each color
    channel {c}, the procedure assumes that the surface is illuminated
    by a light field that is a mix of uniform an isotropic
    (pan-directional, "ambient") light and the light from a white point
    source at infinte distance, with intensities {cIso[c]} with weights
    and {(1-cIso[c])*cLit[c]}, respectively.
    
    Specifically, {cIso[c]} is taken to be be the apparent brightness
    that the isotropic light produces on a Lambertian white surface,
    with any viewing directions on the outwards side of the surface, if
    the directional light intensity {cLit[c]} were zero. If {cIso} is
    {NULL}, it is assumed to be {(0,0,0)} (no isotropic light).
    
    The value {cLit[c]} is taken to be apparent brightness that the
    point source produces when falling on a Lambertian white surface
    perpendicular to {uLit} and viewed from any direction on the same
    side as {uLit}, if the isotropic light intensity {cIso[c]} were
    zero. If {uLit} and/or {cLit} is {NULL}, {cLit} is assumed to be
    {(0,0,0)} (no unidirectional light).
    
    For each color channel {c}, the surface is assumed to have mix of a
    white glossy BRDF and a white Lambertian BRDF, with weights
    {sGlo[c]} and {(1-sGlo[c])*sLam[c]}, respectively. If {sGlo} is
    {NULL}, it is assumed to be {(0,0,0)}. If {sLam} is {NULL}, it is
    assumed to be {(0,0,0)}.
    
    In particular, if {cLit} is {(0,0,0)}, {cIso} is {(1,1,1)}, and
    {sGlo} is {(0,0,0)}, the rendered color surface color will be just
    the surface's lambedo {sLam}, independent of the surface's normal
    direction. */

#endif
