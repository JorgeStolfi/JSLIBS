/* See {neuromat_{eeg_|}image.h}. */
/* Last edited on 2021-08-31 12:11:20 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <r2.h>
#include <r3.h>
#include <affirm.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_geom.h>

void neuromat_eeg_geom_get_20_schematic_2D_points(int32_t ne, r2_t pos2D[]);
  /* Puts into {posdD[0..ne-1]} the  two-dimensional
    schematic positions of the electrodes in 20-electrode
    experiments.  The value of {ne} must be 20. */

void neuromat_eeg_geom_get_128_schematic_2D_points(int32_t ne, r2_t pos2D[]);
  /* Puts into {posdD[0..ne-1]} the two-dimensional schematic positions
    of the electrodes in 128-electrode experiments. The value of {ne}
    must be 128 or 129; if the latter, includes the voltage reference
    electrode ("CZ") as the last element. */
    
r2_t *neuromat_eeg_geom_get_schematic_2D_points(int32_t ne)
  {
    r2_t *pos2D = (r2_t *)notnull(malloc(ne*sizeof(r2_t)), "no mem");
    if (ne == 20)
      { neuromat_eeg_geom_get_20_schematic_2D_points(ne, pos2D); }
    else if ((ne == 128) || (ne == 129))
      { neuromat_eeg_geom_get_128_schematic_2D_points(ne, pos2D); }
    else
      { demand(FALSE, "invalid electrode count"); }
    return pos2D;
  }

    
r2_t *neuromat_eeg_geom_get_schematic_2D_points_by_name(int32_t ne, char *chnames[], int32_t ne_full)
  { r2_t *pos2D_full = neuromat_eeg_geom_get_schematic_2D_points(ne_full);
    char **chnames_full = neuromat_eeg_get_channel_names(ne_full, 0, NULL); /* Names of full electrode set. */
    r2_t *pos2D = (r2_t *)notnull(malloc(ne*sizeof(r2_t)), "no mem");
    for (int32_t ie = 0; ie < ne; ie++)
      { char *name = chnames[ie];
        int32_t ie_full = neuromat_eeg_find_channel_by_name(name, 0, ne_full-1, chnames_full, TRUE);
        pos2D[ie] = pos2D_full[ie_full];
      }
    for (int32_t ie_full = 0; ie_full < ne_full; ie_full++) { free(chnames_full[ie_full]); }
    free(chnames_full);
    free(pos2D_full);
    return pos2D;
  }


void neuromat_eeg_geom_get_20_schematic_2D_points(int32_t ne, r2_t pos2D[])
  {
    demand(ne == 20, "invalid number of electrodes");
    /* Positions in arbitrary coordinates:  */
    r2_t ctr_ref = (r2_t){{ 334, 304 }};
    r2_t rad_ref = (r2_t){{ 250, 250 }};
    pos2D[ 0] = (r2_t){{ 156, 193 }}; /* F7 */ 
    pos2D[ 1] = (r2_t){{ 121, 304 }}; /* T3 */ 
    pos2D[ 2] = (r2_t){{ 162, 429 }}; /* T5 */ 
    pos2D[ 3] = (r2_t){{ 265, 102 }}; /* Fp1 */
    pos2D[ 4] = (r2_t){{ 244, 204 }}; /* F3 */ 
    pos2D[ 5] = (r2_t){{ 224, 304 }}; /* C3 */ 
    pos2D[ 6] = (r2_t){{ 248, 415 }}; /* P3 */ 
    pos2D[ 7] = (r2_t){{ 263, 502 }}; /* O1 */ 
    pos2D[ 8] = (r2_t){{ 512, 193 }}; /* F8 */ 
    pos2D[ 9] = (r2_t){{ 537, 304 }}; /* T4 */ 
    pos2D[10] = (r2_t){{ 506, 429 }}; /* T6 */ 
    pos2D[11] = (r2_t){{ 403, 102 }}; /* Fp2 */
    pos2D[12] = (r2_t){{ 424, 204 }}; /* F4 */ 
    pos2D[13] = (r2_t){{ 444, 304 }}; /* C4 */ 
    pos2D[14] = (r2_t){{ 418, 415 }}; /* P4 */ 
    pos2D[15] = (r2_t){{ 405, 502 }}; /* O2 */ 
    pos2D[16] = (r2_t){{ 334, 207 }}; /* Fz */ 
    pos2D[17] = (r2_t){{ 334, 304 }}; /* Cz */ 
    pos2D[18] = (r2_t){{ 334, 412 }}; /* Pz */ 
    pos2D[19] = (r2_t){{ 334, 509 }}; /* Oz */ 
    
    /* Rescale to unit disk: */
    for (int32_t ie = 0; ie < ne; ie++) 
      { for (int32_t j = 0; j < 2; j++) 
          { double z_raw = pos2D[ie].c[j];
            double z_rel = (z_raw - ctr_ref.c[j])/rad_ref.c[j];
            if (j == 1) { z_rel = -z_rel; }
            pos2D[ie].c[j] = z_rel;
          }
      }
  }

void neuromat_eeg_geom_get_128_schematic_2D_points(int32_t ne, r2_t pos2D[])
  {
    demand((ne = 128) || (ne == 129), "invalid {ne}");

    /* Positions in arbitrary coordinates:  */
    r2_t ctr_ref = (r2_t){{ 320, 320 }};
    r2_t rad_ref = (r2_t){{ 320, 320 }};
    
    auto r2_t sym(r2_t p);
      /* Mirrors {p} across the midline. */
    
    r2_t sym(r2_t p) 
      { return (r2_t){{ 640 - p.c[0], p.c[1] }}; }
    
    /* Midline: */
    pos2D[ 16] = (r2_t){{ 320,  62 }}; /* "C17" medial */
    pos2D[ 14] = (r2_t){{ 320, 121 }}; /* "C15" medial */
    pos2D[ 15] = (r2_t){{ 320, 160 }}; /* "C16" medial */
    pos2D[ 10] = (r2_t){{ 320, 189 }}; /* "C11" medial */
    pos2D[  5] = (r2_t){{ 320, 256 }}; /* "C6" medial */
    pos2D[ 54] = (r2_t){{ 320, 352 }}; /* "C55" medial */
    pos2D[ 61] = (r2_t){{ 320, 420 }}; /* "C62" medial */
    pos2D[ 71] = (r2_t){{ 320, 450 }}; /* "C72" medial */
    pos2D[ 74] = (r2_t){{ 320, 513 }}; /* "C75" medial */
    pos2D[ 80] = (r2_t){{ 320, 577 }}; /* "C81" medial */
    /* Right half: */
    pos2D[  0] = (r2_t){{ 477, 162 }}; /* "C1" -> "C32" */
    pos2D[  1] = (r2_t){{ 432, 165 }}; /* "C2" -> "C26" */
    pos2D[  2] = (r2_t){{ 396, 176 }}; /* "C3" -> "C23" */
    pos2D[  3] = (r2_t){{ 364, 198 }}; /* "C4" -> "C19" */
    pos2D[  4] = (r2_t){{ 340, 226 }}; /* "C5" -> "C12" */
    pos2D[  7] = (r2_t){{ 430, 127 }}; /* "C8" -> "C25" */
    pos2D[  8] = (r2_t){{ 393, 141 }}; /* "C9" -> "C22" */
    pos2D[  9] = (r2_t){{ 355, 164 }}; /* "C10" -> "C18" */
    pos2D[ 13] = (r2_t){{ 369, 100 }}; /* "C14" -> "C21" */
    pos2D[ 75] = (r2_t){{ 348, 479 }}; /* "C76" -> "C71" */
    pos2D[ 76] = (r2_t){{ 358, 442 }}; /* "C77" -> "C67" */
    pos2D[ 77] = (r2_t){{ 360, 409 }}; /* "C78" -> "C61" */
    pos2D[ 78] = (r2_t){{ 359, 371 }}; /* "C79" -> "C54" */
    pos2D[ 79] = (r2_t){{ 348, 330 }}; /* "C80" -> "C31" */
    pos2D[ 81] = (r2_t){{ 358, 541 }}; /* "C82" -> "C74" */
    pos2D[ 82] = (r2_t){{ 378, 503 }}; /* "C83" -> "C70" */
    pos2D[ 83] = (r2_t){{ 389, 464 }}; /* "C84" -> "C66" */
    pos2D[ 84] = (r2_t){{ 392, 424 }}; /* "C85" -> "C60" */
    pos2D[ 85] = (r2_t){{ 388, 384 }}; /* "C86" -> "C53" */
    pos2D[ 86] = (r2_t){{ 378, 340 }}; /* "C87" -> "C37" */
    pos2D[ 87] = (r2_t){{ 397, 565 }}; /* "C88" -> "C73" */
    pos2D[ 88] = (r2_t){{ 421, 519 }}; /* "C89" -> "C69" */
    pos2D[ 89] = (r2_t){{ 433, 474 }}; /* "C90" -> "C65" */
    pos2D[ 90] = (r2_t){{ 432, 431 }}; /* "C91" -> "C59" */
    pos2D[ 91] = (r2_t){{ 425, 395 }}; /* "C92" -> "C52" */
    pos2D[ 92] = (r2_t){{ 409, 350 }}; /* "C93" -> "C42" */
    pos2D[ 93] = (r2_t){{ 468, 527 }}; /* "C94" -> "C68" */
    pos2D[ 94] = (r2_t){{ 474, 476 }}; /* "C95" -> "C64" */
    pos2D[ 95] = (r2_t){{ 475, 425 }}; /* "C96" -> "C58" */
    pos2D[ 96] = (r2_t){{ 460, 391 }}; /* "C97" -> "C51" */
    pos2D[ 97] = (r2_t){{ 436, 361 }}; /* "C98" -> "C47" */
    pos2D[ 98] = (r2_t){{ 523, 469 }}; /* "C99" -> "C63" */
    pos2D[ 99] = (r2_t){{ 518, 415 }}; /* "C100" -> "C57" */
    pos2D[100] = (r2_t){{ 498, 375 }}; /* "C101" -> "C50" */
    pos2D[101] = (r2_t){{ 474, 346 }}; /* "C102" -> "C46" */
    pos2D[102] = (r2_t){{ 442, 324 }}; /* "C103" -> "C41" */
    pos2D[103] = (r2_t){{ 412, 304 }}; /* "C104" -> "C36" */
    pos2D[104] = (r2_t){{ 378, 297 }}; /* "C105" -> "C30" */
    pos2D[105] = (r2_t){{ 339, 293 }}; /* "C106" -> "C7" */
    pos2D[106] = (r2_t){{ 558, 393 }}; /* "C107" -> "C56" */
    pos2D[107] = (r2_t){{ 506, 326 }}; /* "C108" -> "C45" */
    pos2D[108] = (r2_t){{ 471, 297 }}; /* "C109" -> "C40" ("T7"?) */
    pos2D[109] = (r2_t){{ 437, 279 }}; /* "C110" -> "C35" */
    pos2D[110] = (r2_t){{ 398, 269 }}; /* "C111" -> "C29" */
    pos2D[111] = (r2_t){{ 352, 267 }}; /* "C112" -> "C13" */
    pos2D[112] = (r2_t){{ 567, 328 }}; /* "C113" -> "C49" */
    pos2D[113] = (r2_t){{ 534, 287 }}; /* "C114" -> "C44" */
    pos2D[114] = (r2_t){{ 496, 261 }}; /* "C115" -> "C39" */
    pos2D[115] = (r2_t){{ 459, 250 }}; /* "C116" -> "C34" */
    pos2D[116] = (r2_t){{ 416, 241 }}; /* "C117" -> "C28" */
    pos2D[117] = (r2_t){{ 374, 242 }}; /* "C118" -> "C20" */
    pos2D[118] = (r2_t){{ 592, 234 }}; /* "C119" -> "C48" */
    pos2D[119] = (r2_t){{ 554, 237 }}; /* "C120" -> "C43" */
    pos2D[120] = (r2_t){{ 512, 215 }}; /* "C121" -> "C38" */
    pos2D[121] = (r2_t){{ 471, 207 }}; /* "C122" -> "C33" */
    pos2D[122] = (r2_t){{ 430, 208 }}; /* "C123" -> "C27" */
    pos2D[123] = (r2_t){{ 394, 218 }}; /* "C124" -> "C24" */
    pos2D[124] = (r2_t){{ 511, 155 }}; /* "C125" -> "C128" */
    pos2D[125] = (r2_t){{ 357,  40 }}; /* "C126" -> "C127" */
    /* Left half: */
    pos2D[  6] = sym(pos2D[106-1]); /* "C7" <- "C106" */
    pos2D[ 11] = sym(pos2D[  5-1]); /* "C12" <- "C5" */
    pos2D[ 12] = sym(pos2D[112-1]); /* "C13" <- "C112" */
    pos2D[ 17] = sym(pos2D[ 10-1]); /* "C18" <- "C10" */
    pos2D[ 18] = sym(pos2D[  4-1]); /* "C19" <- "C4" */
    pos2D[ 19] = sym(pos2D[118-1]); /* "C20" <- "C118" */
    pos2D[ 20] = sym(pos2D[ 14-1]); /* "C21" <- "C14" */
    pos2D[ 21] = sym(pos2D[  9-1]); /* "C22" <- "C9" */
    pos2D[ 22] = sym(pos2D[  3-1]); /* "C23" <- "C3" */
    pos2D[ 23] = sym(pos2D[124-1]); /* "C24" <- "C124" */
    pos2D[ 24] = sym(pos2D[  8-1]); /* "C25" <- "C8" */
    pos2D[ 25] = sym(pos2D[  2-1]); /* "C26" <- "C2" */
    pos2D[ 26] = sym(pos2D[123-1]); /* "C27" <- "C123" */
    pos2D[ 27] = sym(pos2D[117-1]); /* "C28" <- "C117" */
    pos2D[ 28] = sym(pos2D[111-1]); /* "C29" <- "C111" */
    pos2D[ 29] = sym(pos2D[105-1]); /* "C30" <- "C105" */
    pos2D[ 30] = sym(pos2D[ 80-1]); /* "C31" <- "C80" */
    pos2D[ 31] = sym(pos2D[  1-1]); /* "C32" <- "C1" */
    pos2D[ 32] = sym(pos2D[122-1]); /* "C33" <- "C122" */
    pos2D[ 33] = sym(pos2D[116-1]); /* "C34" <- "C116" */
    pos2D[ 34] = sym(pos2D[110-1]); /* "C35" <- "C110" */
    pos2D[ 35] = sym(pos2D[104-1]); /* "C36" <- "C104" */
    pos2D[ 36] = sym(pos2D[ 87-1]); /* "C37" <- "C87" */
    pos2D[ 37] = sym(pos2D[121-1]); /* "C38" <- "C121" */
    pos2D[ 38] = sym(pos2D[115-1]); /* "C39" <- "C115" */
    pos2D[ 39] = sym(pos2D[109-1]); /* "C40" <- "C109" ("T7"?) */
    pos2D[ 40] = sym(pos2D[103-1]); /* "C41" <- "C103" */
    pos2D[ 41] = sym(pos2D[ 93-1]); /* "C42" <- "C93" */
    pos2D[ 42] = sym(pos2D[120-1]); /* "C43" <- "C120" */
    pos2D[ 43] = sym(pos2D[114-1]); /* "C44" <- "C114" */
    pos2D[ 44] = sym(pos2D[108-1]); /* "C45" <- "C108" */
    pos2D[ 45] = sym(pos2D[102-1]); /* "C46" <- "C102" */
    pos2D[ 46] = sym(pos2D[ 98-1]); /* "C47" <- "C98" */
    pos2D[ 47] = sym(pos2D[119-1]); /* "C48" <- "C119" */
    pos2D[ 48] = sym(pos2D[113-1]); /* "C49" <- "C113" */
    pos2D[ 49] = sym(pos2D[101-1]); /* "C50" <- "C101" */
    pos2D[ 50] = sym(pos2D[ 97-1]); /* "C51" <- "C97" */
    pos2D[ 51] = sym(pos2D[ 92-1]); /* "C52" <- "C92" */
    pos2D[ 52] = sym(pos2D[ 86-1]); /* "C53" <- "C86" */
    pos2D[ 53] = sym(pos2D[ 79-1]); /* "C54" <- "C79" */
    pos2D[ 55] = sym(pos2D[107-1]); /* "C56" <- "C107" */
    pos2D[ 56] = sym(pos2D[100-1]); /* "C57" <- "C100" */
    pos2D[ 57] = sym(pos2D[ 96-1]); /* "C58" <- "C96" */
    pos2D[ 58] = sym(pos2D[ 91-1]); /* "C59" <- "C91" */
    pos2D[ 59] = sym(pos2D[ 85-1]); /* "C60" <- "C85" */
    pos2D[ 60] = sym(pos2D[ 78-1]); /* "C61" <- "C78" */
    pos2D[ 62] = sym(pos2D[ 99-1]); /* "C63" <- "C99" */
    pos2D[ 63] = sym(pos2D[ 95-1]); /* "C64" <- "C95" */
    pos2D[ 64] = sym(pos2D[ 90-1]); /* "C65" <- "C90" */
    pos2D[ 65] = sym(pos2D[ 84-1]); /* "C66" <- "C84" */
    pos2D[ 66] = sym(pos2D[ 77-1]); /* "C67" <- "C77" */
    pos2D[ 67] = sym(pos2D[ 94-1]); /* "C68" <- "C94" */
    pos2D[ 68] = sym(pos2D[ 89-1]); /* "C69" <- "C89" */
    pos2D[ 69] = sym(pos2D[ 83-1]); /* "C70" <- "C83" */
    pos2D[ 70] = sym(pos2D[ 76-1]); /* "C71" <- "C76" */
    pos2D[ 72] = sym(pos2D[ 88-1]); /* "C73" <- "C88" */
    pos2D[ 73] = sym(pos2D[ 82-1]); /* "C74" <- "C82" */
    pos2D[126] = sym(pos2D[126-1]); /* "C127" <- "C126" */
    pos2D[127] = sym(pos2D[125-1]); /* "C128" <- "C125" */

    if (ne == 129)
      { /* Add the voltage reference "CZ": */
        pos2D[128] = (r2_t){{ 320, 320 }}; /* "CZ" */
      }
    
    /* Rescale to unit disk: */
    for (int32_t ie = 0; ie < ne; ie++) 
      { for (int32_t j = 0; j < 2; j++) 
          { double z_raw = pos2D[ie].c[j];
            double z_rel = (z_raw - ctr_ref.c[j])/rad_ref.c[j];
            if (j == 1) { z_rel = -z_rel; }
            pos2D[ie].c[j] = z_rel;
          }
      }
  }
 
r3_t neuromat_eeg_geom_3D_from_2D(r2_t *p, r2_t *rad2, r3_t *rad3)
  {
    r2_t q2 = (*p);
    if (rad2 != NULL) { r2_unweigh(&q2, rad2, &q2); }
    double R2 = r2_norm_sqr(&q2);
    double Den = 1 + R2;
    r3_t q3 = (r3_t){{ 2*q2.c[0]/Den, 2*q2.c[1]/Den, (R2 - 1)/Den }};
    /* Paranoia: */
    (void)r3_dir(&q3, &q3);
    if (rad3 != NULL) { r3_weigh(&q3, rad3, &q3); }
    return q3;
  }

r2_t neuromat_eeg_geom_2D_from_3D(r3_t *p, r3_t *rad3, r2_t *rad2)
  {
    r3_t q3 = (*p);
    if (rad3 != NULL) { r3_unweigh(&q3, rad3, &q3); }
    double Den = 1 + q3.c[2];
    r2_t q2 = (r2_t){{ q3.c[0]/Den, q3.c[1]/Den }};
    if (rad2 != NULL) { r2_weigh(&q2, rad2, &q2); }
    return q2;
  }

r2_t neuromat_eeg_geom_disk_from_ellipse(r2_t *p, r2_t *ctr, r2_t *rad)
  {
    return (r2_t) 
      {{  (p->c[0] - ctr->c[0])/rad->c[0],
          (p->c[1] - ctr->c[1])/rad->c[1]
      }};
  }

r2_t neuromat_eeg_geom_ellipse_from_disk(r2_t *p, r2_t *ctr, r2_t *rad)
  {
    return (r2_t) 
      {{  ctr->c[0] + p->c[0]*rad->c[0],
          ctr->c[1] + p->c[1]*rad->c[1]
      }};
  }

void neuromat_eeg_geom_map_many_disk_to_ellipse(int32_t np, r2_t p[], r2_t *ctr, r2_t *rad, r2_t q[])
  {
    for (int32_t ip = 0; ip < np; ip++) 
      { q[ip] = neuromat_eeg_geom_ellipse_from_disk(&(p[ip]), ctr, rad); }
  }

void neuromat_eeg_geom_map_many_ellipse_to_disk(int32_t np, r2_t p[], r2_t *ctr, r2_t *rad, r2_t q[])
  {
    for (int32_t ip = 0; ip < np; ip++) 
      { q[ip] = neuromat_eeg_geom_disk_from_ellipse(&(p[ip]), ctr, rad); }
  }
