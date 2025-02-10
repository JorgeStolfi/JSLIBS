#define PROG_NAME "test_mfok_object_raytrace_normal"
#define PROG_DESC "test of {multifok_scene_object_raytrace.h},  {multifok_scene_object_normal.h}"
#define PROG_VERS "1.0"

/* Last edited on 2025-02-10 07:52:17 by stolfi */ 
/* Created on 2023-01-05 by J. Stolfi, UNICAMP */

#define test_mfok_object_normal_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <interval.h>
#include <box.h>
#include <r3.h>

#include <multifok_scene_object.h>
#include <multifok_scene_object_raytrace.h>
#include <multifok_scene_object_normal.h>

#define PROG_HELP \
  "  Please help yourself!"

#define PROG_INFO \
  "  There are five fingers in one hand."
  
#define ot_FLAT multifok_scene_object_type_FLAT
#define ot_RAMP multifok_scene_object_type_RAMP
#define ot_DISK multifok_scene_object_type_DISK
#define ot_BALL multifok_scene_object_type_BALL
#define ot_CONE multifok_scene_object_type_CONE
#define ot_PYRA multifok_scene_object_type_PYRA

#define ot_FRST multifok_scene_object_type_FRST
#define ot_LAST multifok_scene_object_type_LAST
   
/* !!! Should check also oblique rays !!! */

void morn_do_test(double sObj, double rRay);
  /* Performs the tests with a vertical ray that passes at distance {rRay}
    from the object's center, with nonzero {X} and {y} displacements.
    The object will have nominal radius {sObj}. */
 
void morn_set_bbox(multifok_scene_object_t *obj, r3_t *cObj, double sObj);
  /* Sets {obj.bbox} as appropriate for the given type, assuming the
    center {cObj} and nominal radius {sObj}.
    
    For foreground objects, {sObj} is the {X} and {Y} half-width. 
    Background objects will have {X} and {Y} size many times {sObj}. */

bool_t morn_should_hit(multifok_scene_object_t *obj, r3_t *pRay, r3_t *dRay);
  /* Returns true if and only if the ray with origin {pRay}
    and direction {dRay} is expected to hit the object {obj}.
    Currently works only for vertical downward rays with {pRay} above the object. */
    
void morn_check_hit(bool_t should_hit, multifok_scene_object_t *obj, r3_t *pHit, r3_t *nHit_cmp, r3_t *nHit_num);
  /* To be called if the test ray hit the object {obj} at point {pHit}
    the officially computed normal there is {nHit_cmp} and the numerically computed normal
    there is {nHit_num}. */

r3_t morn_compute_normal_numerically(multifok_scene_object_t *obj, r3_t *pRay, r3_t *dRay, double tMin, double tMax, r3_t *pHit);
  /* Assumes that {pHit} is the hit point on {obj} of ray {pRay,dRay}.
    Traces two additional rays with direction {dRay}, from points slightly perturbed from {pRay}, and returns the
    normal of the triangle whose corners are {pHit} and the hit points if those those two rays.
    If the attempt fails, retruns {(0,0,0)}. */
    
int32_t main(int32_t argn, char **argv);

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  { morn_do_test(2.0, 0.95);
    morn_do_test(2.0, 1.01);
    return 0;
  }

void morn_do_test(double sObj, double rRay)
  { 
    fprintf(stderr, "checking with radius at distance %12.10f ...\n", rRay);

    multifok_scene_object_t obj;
    obj.ID = 1;

    r3_t cObj; /* Logical center of object. */
    for (int32_t j = 0; j < 3; j++) { cObj.c[j] = 3*j; }

    obj.fgGlo = (frgb_t){{ 0.200f, 0.190f, 0.020f }};
    obj.bgGlo = (frgb_t){{ 0.000f, 0.000f, 0.000f }};

    obj.fgLam = (frgb_t){{ 0.950f, 0.750f, 0.200f }};
    obj.bgLam = (frgb_t){{ 0.900f, 0.300f, 0.000f }};
    
    double uAng = M_PI/6;
    double uRay_x = rRay*cos(uAng); 
    double uRay_y = rRay*sin(uAng);

    for (int32_t type = ot_FRST; type <= ot_LAST; type++)
      { fprintf(stderr, "\n");
        
        obj.type = (multifok_scene_object_type_t)type;
        char *xtype = multifok_scene_object_type_to_string(obj.type);
        fprintf(stderr, "  checking type %s ...\n", xtype);

        morn_set_bbox(&obj, &cObj, sObj);
      
        double zBot = obj.bbox[2].end[0]; 
        double zTop = obj.bbox[2].end[1]; 
      
        r3_t pRay = (r3_t){{ cObj.c[0] + uRay_x, cObj.c[1] + uRay_y, zTop + 0.1 }}; /* Ray origin. */
        r3_t dRay = (r3_t){{ 0.0, 0.0, -1.0 }}; /* Ray direction. */

        bool_t should_hit =  morn_should_hit(&obj, &pRay, &dRay);
        double tMin = (zTop - pRay.c[2])/dRay.c[2] - 1.0e-6;
        double tMax = (zBot - pRay.c[2])/dRay.c[2] + 1.0e-6;
        assert(tMax > tMin);
        double tHit = multifok_scene_object_raytrace(&obj, &pRay, &dRay, tMin, tMax, TRUE);
        if (isfinite(tHit))
          { r3_t pHit; r3_mix(1.0, &pRay, tHit, &dRay, &pHit);
            r3_t nHit_cmp = multifok_scene_object_normal(&obj, &pHit, TRUE);
            r3_t nHit_num = morn_compute_normal_numerically(&obj, &pRay, &dRay, tMin, tMax, &pHit);
            morn_check_hit(should_hit, &obj, &pHit, &nHit_cmp, &nHit_num);
            demand(should_hit, "ray hit but should not have");
          }
        else
          { demand((! should_hit), "ray should have hit but did not"); }
          
        fprintf(stderr, "\n");
      }
  }
  
void morn_set_bbox(multifok_scene_object_t *obj, r3_t *cObj, double sObj)
  {
    for (int32_t j = 0; j < 3; j++)
      { interval_t *bbj = &(obj->bbox[j]);
        switch(obj->type)
          { case ot_FLAT:
              bbj->end[0] = cObj->c[j] - (j == 2 ? 0.0 : 5*sObj);
              bbj->end[1] = cObj->c[j] + (j == 2 ? 0.0 : 5*sObj);
              break;
            
            case ot_RAMP:
              bbj->end[0] = cObj->c[j] - (j == 2 ? sObj : 5*sObj);
              bbj->end[1] = cObj->c[j] + (j == 2 ? sObj : 5*sObj);
              break;
            
            case ot_DISK:
              bbj->end[0] = cObj->c[j] - (j == 2 ? 0.0 : sObj);
              bbj->end[1] = cObj->c[j] + (j == 2 ? 0.0 : sObj);
              break;
            
            case ot_BALL:
            case ot_CONE:
            case ot_PYRA:
              bbj->end[0] = cObj->c[j] - sObj;
              bbj->end[1] = cObj->c[j] + sObj;
              break;
              
            default:
              assert(FALSE);
          }
      }
  }
            
bool_t morn_should_hit(multifok_scene_object_t *obj, r3_t *pRay, r3_t *dRay)
  { demand((dRay->c[0] == 0) && (dRay->c[1] == 0), "oblique rays not tested yet");

    r3_t cObj, rObj;
    for (int32_t j = 0; j < 3; j++)
      { interval_mid_rad(&(obj->bbox[j]), &(cObj.c[j]), &(rObj.c[j])); }
      
    r3_t uRay; r3_sub(pRay, &cObj, &uRay);

    switch(obj->type)
      { case ot_FLAT:
        case ot_RAMP:
          return TRUE;
          
        case ot_DISK:
        case ot_BALL:
        case ot_CONE:
          demand(fabs(rObj.c[0] - rObj.c[1]) < 1.0e-10, "DISK/BALL/CONE with unequal X and Y radii");
          return (hypot(uRay.c[0], uRay.c[1]) < rObj.c[0]);
        
        case ot_PYRA:
          demand(fabs(rObj.c[0] - rObj.c[1]) < 1.0e-10, "PYRA with unequal X and Y radii");
          return (fmax(fabs(uRay.c[0]), fabs(uRay.c[1])) < rObj.c[0]);
          
        default:
          assert(FALSE);
      }
  }

r3_t morn_compute_normal_numerically(multifok_scene_object_t *obj, r3_t *pRay, r3_t *dRay, double tMin, double tMax, r3_t *pHit)
  { double eps = 1.0e-5;
    r3_t u_pert[2];
    for (int32_t j = 0; j < 2; j++)
      { /* Perturb the ray slightly along axis {j}: */
        bool_t hit_pert = FALSE;
        for (int32_t sgn = -1; (sgn <= +1) && (! hit_pert); sgn += 2)
          { r3_t pRay_pert = (*pRay);
            pRay_pert.c[j] += sgn*eps;
            double tHit_pert = multifok_scene_object_raytrace(obj, &pRay_pert, dRay, tMin, tMax, TRUE);
            if (isfinite(tHit_pert))
              { r3_t pHit_pert; r3_mix(1.0, &pRay_pert, tHit_pert, dRay, &(pHit_pert));
                r3_sub(&pHit_pert, pHit, &(u_pert[j]));
                hit_pert = TRUE;
              }
          }
        demand(hit_pert, "cannot hit with a perturbed ray");
      }
    r3_t nHit_num;
    r3_cross(&(u_pert[0]), &(u_pert[1]), &nHit_num);
    (void)r3_dir(&nHit_num, &nHit_num);
    if (r3_dot(&nHit_num, dRay) > 0) { r3_neg(&nHit_num, &nHit_num); }
    return nHit_num;
  }

void morn_check_hit(bool_t should_hit, multifok_scene_object_t *obj, r3_t *pHit, r3_t *nHit_cmp, r3_t *nHit_num)  
  {
    r3_t cObj, rObj;
    for (int32_t j = 0; j < 3; j++)
      { interval_mid_rad(&(obj->bbox[j]), &(cObj.c[j]), &(rObj.c[j])); }
    double xyrad = hypot(rObj.c[0], rObj.c[1]);

    r3_gen_print(stderr, pHit, "%24.16e", "    hit point (abs) = ( ", " ", " )\n");

    r3_t uHit; r3_sub(pHit, &cObj, &uHit);
    r3_gen_print(stderr, &uHit, "%24.16e", "    hit point (rel) = ( ", " ", " )");
    r3_t vHit; r3_scale(1/xyrad, &uHit, &vHit);
    r3_gen_print(stderr, &vHit, "%24.16e", " = rXY * ( ", " ", " )\n");

    demand(should_hit, "should not have hit but did");

    r3_gen_print(stderr, nHit_cmp, "%24.16e", "    normal (computed) =  ( ", " ", " )");
    double nMag_cmp = r3_norm(nHit_cmp);
    fprintf(stderr, " length = %24.16e\n", nMag_cmp);
    demand(fabs(nMag_cmp - 1.0) < 1.0e-10, "normal is not normalized");

    r3_gen_print(stderr, nHit_num, "%24.16e", "    normal (numerical) = ( ", " ", " )");
    double nMag_num = r3_norm(nHit_num);
    fprintf(stderr, " length = %24.16e\n", nMag_num);
    
    if (nMag_num == 0)
      { fprintf (stderr, "    !! numeric normal failed\n"); }
    else
      { double nErr = r3_dist(nHit_cmp, nHit_num);
        fprintf(stderr, "    normal error = %24.16e\n", nErr);
        demand(nErr <= 1.0e-4, "normals do not match");
      }
  }
