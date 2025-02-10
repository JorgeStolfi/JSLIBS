/* See {multifok_render.h}. */
/* Last edited on 2025-02-06 23:48:48 by stolfi */

#include <stdint.h>
#include <math.h>

#include <bool.h>
#include <r3.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <affirm.h>

#include <multifok_render.h>

/* INTERNAL PROTOTYPES */

frgb_t multifok_render_compute_color_unidir(r3_t *sNrm, r3_t *uVis, r3_t *uLit, frgb_t *eLit, frgb_t *sGlo, frgb_t *sLam, bool_t debug);
  /* Computes the apparent color {clr} of a surface due to a distant
    point light source, given unit normal vector {sNrm} pointing out of
    the surface, the unit direction vector {uVis} pointing towards the
    observer, the unit direction vector {uLit} pointing towards the
    light source, the intensity {eLit} of the light at the surface, and
    the mixing coefficients {sLam} (albedo) and {sGlo} of the Lambertian
    and Glossy BRDF, respectively.
    
    In each color channel {c}, the result {clr[c]} will be
    {eLit[c]*sLam[c]} if {sNrm=uLit} and {sGlo[c]} is zero. */ 
  
frgb_t multifok_render_compute_color_isotropic(r3_t *sNrm, r3_t *uVis, frgb_t *eIso, frgb_t *sGlo, frgb_t *sLam, bool_t debug);
  /* Computes the apparent color {clr} of a surface due to a isotropic (pan-directional, "ambient")
    light field, given unit normal vector {sNrm} pointing out of
    the surface, the unit direction vector {uVis} pointing towards the
    observer, the intensity {eIso} of the light at the surface, and
    the mixing coefficients {sLam} (albedo) and {sGlo} of the Lambertian
    and Glossy BRDF, respectively.
    
    In each color channel {c}, the result {clr[c]} will be
    {eIso[c]*sLam[c]} if {sNrm=uLit} and {sGlo[c]} is zero. */ 

/* IMPLEMENTATIONS */

frgb_t multifok_render_compute_color
  ( r3_t *uVis,
    r3_t *sNrm, 
    frgb_t *sGlo,     
    frgb_t *sLam, 
    r3_t *uLit,
    frgb_t *cLit, 
    frgb_t *cIso,
    bool_t debug
  )
  { demand(uVis != NULL, "null view direction {uVis}");
    double mNrm = (sNrm == NULL ? 0.0 : r3_norm(sNrm));
    demand(isfinite(mNrm), "invalid normal direction");
    if (mNrm < 1.0e-20) 
      { /* Undefined norma; return generic gray. */
        if (debug) { fprintf(stderr, "        -- null normal\n"); }
        return (frgb_t){{ 0.500, 0.500, 0.500 }};
      }
    double dotNrmVis = - r3_dot(sNrm, uVis);
    demand(dotNrmVis >= -1.0e-8, "hit point on back side of object");
    if (dotNrmVis <= 0)
      { /* Should not happen, but in any case: */
        if (debug) { fprintf(stderr, "        -- glancing hit\n"); }
        return (frgb_t){{ 0.000, 0.000, 0.000 }};
      }
    double dotNrmLit = (uLit == NULL ? 0.0 : r3_dot(sNrm, uLit));
    if (debug) { fprintf(stderr, "        -- dor(sNrm,uLit) = %+8.5f\n", dotNrmLit); }
    frgb_t eIso, eLit; /* Local illumination energy per area from isotropic and point source. */
    for (int32_t c = 0; c < 3; c++)
      { demand((cIso->c[c] >= 0.0) && (cIso->c[c] <= 1.0), "invalid {cIso} fraction");
        eIso.c[c] = cIso->c[c];
        eLit.c[c] = (float)((1 - cIso->c[c])*(dotNrmLit <= 0 ? 0.0 : cLit->c[c]*dotNrmLit));
      }
    if (debug) 
      { frgb_print(stderr, "        -- eIso = ( ", &eIso, 3, "%6.4f", " )\n"); 
        frgb_print(stderr, "        -- eLit = ( ", &eLit, 3, "%6.4f", " )\n"); 
      }
    /* Compute brightness seen from {uVis} from isotropic and point source: */
    frgb_t bLit = multifok_render_compute_color_unidir(sNrm, uVis, uLit, &eLit, sGlo, sLam, debug);
    frgb_t bIso = multifok_render_compute_color_isotropic(sNrm, uVis, &eIso, sGlo, sLam, debug); 
    if (debug) 
      { frgb_print(stderr, "        -- bIso = ( ", &bIso, 3, "%6.4f", " )\n"); 
        frgb_print(stderr, "        -- bLit = ( ", &bLit, 3, "%6.4f", " )\n"); 
      }
    /* Add the two: */
    frgb_t clr;
    for (int32_t c = 0; c < 3; c++) { clr.c[c] = bIso.c[c] + bLit.c[c]; }
    return clr;
  }
  
#define GLO_A 0.9
  /* Transmission of glossy layer for normal light. */
  
#define GLO_C 4.5
  /* Scaling factor for unidir light and glossy BRDF. !!! Should be determined by integral !!! */
  
#define GLO_D 1.5
  /* Scaling factor for ambient light and glossy BRDF. !!! Should be determined by integral !!! */
  
frgb_t multifok_render_compute_color_unidir(r3_t *sNrm, r3_t *uVis, r3_t *uLit, frgb_t *eLit, frgb_t *sGlo, frgb_t *sLam, bool_t debug)
  { 
    demand(sNrm != NULL, "null surface normal {sNrm}");
    demand(uVis != NULL, "null view direction {uVis}");
    demand(eLit != NULL, "null unidirectional light intensty {eLit}");

    if (debug) { fprintf(stderr, "        entering {multifok_render_compute_color_unidir}\n"); }

    if ((uLit == NULL) || ((sLam == NULL) && (sGlo == NULL)) || frgb_is_all_zeros(eLit)) 
      { if (debug) { fprintf(stderr, "          ~~ no brightness from unidir source\n"); }
        return (frgb_t){{ 0,0,0 }};
      }

    double dotNrmVis = - r3_dot(uVis, sNrm);
    demand(dotNrmVis > 0.0, "hit on back of object");
    
    double fGlo; /* Fraction of {eLit} that is scattered by the glossy layer, apart from {sGlo}: */
    double aGlo; /* Glossy BRDF, apart from the {sGlo} factor. */
    if ((sGlo == NULL) || frgb_is_all_zeros(sGlo))
      { fGlo = 0.0; aGlo = 0.0; }
    else
      { /* Compute the fraction of the light scattered by the glossy layer, apart from {sGlo}: */
        double dotNrmLit = r3_dot(sNrm, uLit);
        if (dotNrmLit <= 0.0)
          { fGlo = 0.0; aGlo = 0.0; }
        else
          { /* Compute the mirrored direction: */
            r3_t uMir; r3_mix(2*dotNrmLit, sNrm, -1.0, uLit, &uMir);
            double dotMirVis = - r3_dot(uVis, &uMir);
            if (dotMirVis <= 0.0)
              { fGlo = 0.0; aGlo = 0.0; }
            else
              { /* Evaluate the glossy BRDF, apart from the {sGlo} factor: */
                fGlo = 1.0 - GLO_A*dotNrmLit;
                aGlo = GLO_C*dotMirVis*dotMirVis*dotNrmVis;
              }
          }
      }
    if (debug) { fprintf(stderr, "          ~~ fGlo = %8.6f  aGlo = %8.6f\n", fGlo, aGlo); }
      
    frgb_t eGlo, eLam, bGlo, bLam;
    for (int32_t c = 0; c < 3; c++)
      { /* Compute the amount of unidir light scattered by the glossy layer: */
        eGlo.c[c] = (float)(sGlo == NULL ? 0.0 : fGlo*sGlo->c[c]*eLit->c[c]);
        /* Compute the amount of unidir light that goes through the glossy layer: */
        eLam.c[c] = eLit->c[c] - eGlo.c[c];
        /* Compute the apparent brightness of the glossy layer: */
        bGlo.c[c] = (float)(aGlo*eGlo.c[c]);
        /* Compute the apparent brightness of the Lambertian layer: */
        bLam.c[c] = (float)(sLam == NULL ? 0.0 : eLam.c[c]*sLam->c[c]);
      }
    /* Combine {bGlo,bLam}: */
    frgb_t clr;
    for (int32_t c = 0; c < 3; c++) { clr.c[c] = (float)(bGlo.c[c] + bLam.c[c]); }

    if (debug) 
      { frgb_print(stderr, "          ~~ eGlo = ( ", &eGlo, 3, "%8.6f", " )");
        frgb_print(stderr, " eLam = ( ", &eLam, 3, "%8.6f", " )\n");
        frgb_print(stderr, "          ~~ bGlo = ( ", &bGlo, 3, "%8.6f", " )");
        frgb_print(stderr, " bLam = ( ", &bLam, 3, "%8.6f", " )\n");
        frgb_print(stderr, "          ~~ clr =  ( ", &clr, 3, "%8.6f", " )\n");
      }

    return clr;        
  }
   
frgb_t multifok_render_compute_color_isotropic(r3_t *sNrm, r3_t *uVis, frgb_t *eIso, frgb_t *sGlo, frgb_t *sLam, bool_t debug)
  { 
    demand(sNrm != NULL, "null surface normal {sNrm}");
    demand(uVis != NULL, "null view direction {uVis}");
    demand(eIso != NULL, "null isotropic light intensty {eLit}");

    if (debug) { fprintf(stderr, "        entering {multifok_render_compute_color_isotropic}\n"); }

    if (((sLam == NULL) && (sGlo == NULL)) || frgb_is_all_zeros(eIso)) 
      { if (debug) { fprintf(stderr, "          ~~ no brightness from ambient light\n"); }
        return (frgb_t){{ 0,0,0 }};
      }

    double dotNrmVis = - r3_dot(uVis, sNrm);
    if (dotNrmVis <= 0.0) { return (frgb_t){{ 0,0,0 }}; }
    
    double fGlo; /* Fraction of {eIso} that is scattered by the glossy layer, apart from {sGlo} factor: */
    double aGlo; /* Glossy BRDF, apart from the {sGlo} factor. */
    if ((sGlo == NULL) || frgb_is_all_zeros(sGlo))
      { fGlo = 0.0; aGlo = 0.0; }
    else
      { /* Compute the fraction of the light scattered by the glossy layer, apart from {sGlo}: */
        /* Evaluate the glossy BRDF, apart from the {sGlo} factor: */
        fGlo = 1.0 - 0.5*GLO_A; /* Should be consistent with unidir formula. */
        aGlo = GLO_D*dotNrmVis;
      }
    if (debug) { fprintf(stderr, "          ~~ fGlo = %8.6f  aGlo = %8.6f\n", fGlo, aGlo); }
      
    frgb_t eGlo, eLam, bGlo, bLam;
    for (int32_t c = 0; c < 3; c++)
      { /* Compute the amount of isotropic light scattered by the glossy layer: */
        eGlo.c[c] = (float)(sGlo == NULL ? 0.0 : fGlo*sGlo->c[c]*eIso->c[c]);
        /* Compute the amount of isotropic light that goes through the glossy layer: */
        eLam.c[c] = eIso->c[c] - eGlo.c[c];
        /* Compute the apparent brightness of the glossy layer: */
        bGlo.c[c] = (float)(aGlo*eGlo.c[c]);
        /* Compute the apparent brightness of the Lambertian layer: */
        bLam.c[c] = (float)(sLam == NULL ? 0.0 : eLam.c[c]*sLam->c[c]);
      }
    /* Combine {bGlo,bLam}: */
    frgb_t clr;
    for (int32_t c = 0; c < 3; c++) { clr.c[c] = (float)(bGlo.c[c] + bLam.c[c]); }

    if (debug) 
      { frgb_print(stderr, "          ~~ eGlo = ( ", &eGlo, 3, "%8.6f", " )");
        frgb_print(stderr, " eLam = ( ", &eLam, 3, "%8.6f", " )\n");
        frgb_print(stderr, "          ~~ bGlo = ( ", &bGlo, 3, "%8.6f", " )");
        frgb_print(stderr, " bLam = ( ", &bLam, 3, "%8.6f", " )\n");
        frgb_print(stderr, "          ~~ clr =  ( ", &clr, 3, "%8.6f", " )\n");
      }
      
    return clr;        
  }
