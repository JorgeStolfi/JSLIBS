/* See {neuromat_{eeg_|}image.h}. */
/* Last edited on 2013-11-29 01:15:38 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <r2.h>
#include <r3.h>
#include <affirm.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_geom.h>

r2_t *neuromat_eeg_geom_get_schematic_2D_points(int ne)
  {
    r2_t *pos = NULL;
    if (ne == 20)
      { pos = neuromat_eeg_geom_get_20_schematic_2D_points(); }
    else if (ne == 128)
      { pos = neuromat_eeg_geom_get_128_schematic_2D_points(); }
    else
      { demand(FALSE, "invalid electrode count"); }
    return pos;
  }

r2_t *neuromat_eeg_geom_get_20_schematic_2D_points(void)
  {
    int ne = 20;
    r2_t *pos = notnull(malloc(ne*sizeof(r2_t)), "no mem");
    
    /* Positions in arbitrary coordinates:  */
    r2_t ctr_ref = (r2_t){{ 334, 304 }};
    r2_t rad_ref = (r2_t){{ 250, 250 }};
    pos[ 0] = (r2_t){{ 156, 193 }}; /* F7 */ 
    pos[ 1] = (r2_t){{ 121, 304 }}; /* T3 */ 
    pos[ 2] = (r2_t){{ 162, 429 }}; /* T5 */ 
    pos[ 3] = (r2_t){{ 265, 102 }}; /* Fp1 */
    pos[ 4] = (r2_t){{ 244, 204 }}; /* F3 */ 
    pos[ 5] = (r2_t){{ 224, 304 }}; /* C3 */ 
    pos[ 6] = (r2_t){{ 248, 415 }}; /* P3 */ 
    pos[ 7] = (r2_t){{ 263, 502 }}; /* O1 */ 
    pos[ 8] = (r2_t){{ 512, 193 }}; /* F8 */ 
    pos[ 9] = (r2_t){{ 537, 304 }}; /* T4 */ 
    pos[10] = (r2_t){{ 506, 429 }}; /* T6 */ 
    pos[11] = (r2_t){{ 403, 102 }}; /* Fp2 */
    pos[12] = (r2_t){{ 424, 204 }}; /* F4 */ 
    pos[13] = (r2_t){{ 444, 304 }}; /* C4 */ 
    pos[14] = (r2_t){{ 418, 415 }}; /* P4 */ 
    pos[15] = (r2_t){{ 405, 502 }}; /* O2 */ 
    pos[16] = (r2_t){{ 334, 207 }}; /* Fz */ 
    pos[17] = (r2_t){{ 334, 304 }}; /* Cz */ 
    pos[18] = (r2_t){{ 334, 412 }}; /* Pz */ 
    pos[19] = (r2_t){{ 334, 509 }}; /* Oz */ 
    assert(ne == 20);
    
    /* Rescale to unit disk: */
    int ie,j;
    for (ie = 0; ie < ne; ie++) 
      { for (j = 0; j < 2; j++) 
          { double z_raw = pos[ie].c[j];
            double z_rel = (z_raw - ctr_ref.c[j])/rad_ref.c[j];
            if (j == 1) { z_rel = -z_rel; }
            pos[ie].c[j] = z_rel;
          }
      }

    return pos;
  }

r2_t *neuromat_eeg_geom_get_128_schematic_2D_points(void)
  {
    int ne = 128;
    r2_t *pos = notnull(malloc(ne*sizeof(r2_t)), "no mem");
    
    /* Positions in arbitrary coordinates:  */
    r2_t ctr_ref = (r2_t){{ 320, 320 }};
    r2_t rad_ref = (r2_t){{ 320, 320 }};
    
    auto r2_t sym(r2_t p);
      /* Mirrors {p} across the midline. */
    
    r2_t sym(r2_t p) 
      { return (r2_t){{ 640 - p.c[0], p.c[1] }}; }
    
    /* Midline: */
    pos[ 16] = (r2_t){{ 320,  62 }}; /* "C17" medial */
    pos[ 14] = (r2_t){{ 320, 121 }}; /* "C15" medial */
    pos[ 15] = (r2_t){{ 320, 160 }}; /* "C16" medial */
    pos[ 10] = (r2_t){{ 320, 189 }}; /* "C11" medial */
    pos[  5] = (r2_t){{ 320, 256 }}; /* "C6" medial */
    pos[ 54] = (r2_t){{ 320, 352 }}; /* "C55" medial */
    pos[ 61] = (r2_t){{ 320, 420 }}; /* "C62" medial */
    pos[ 71] = (r2_t){{ 320, 450 }}; /* "C72" medial */
    pos[ 74] = (r2_t){{ 320, 513 }}; /* "C75" medial */
    pos[ 80] = (r2_t){{ 320, 577 }}; /* "C81" medial */
    /* Right half: */
    pos[  0] = (r2_t){{ 477, 162 }}; /* "C1" -> "C32" */
    pos[  1] = (r2_t){{ 432, 165 }}; /* "C2" -> "C26" */
    pos[  2] = (r2_t){{ 396, 176 }}; /* "C3" -> "C23" */
    pos[  3] = (r2_t){{ 364, 198 }}; /* "C4" -> "C19" */
    pos[  4] = (r2_t){{ 340, 226 }}; /* "C5" -> "C12" */
    pos[  7] = (r2_t){{ 430, 127 }}; /* "C8" -> "C25" */
    pos[  8] = (r2_t){{ 393, 141 }}; /* "C9" -> "C22" */
    pos[  9] = (r2_t){{ 355, 164 }}; /* "C10" -> "C18" */
    pos[ 13] = (r2_t){{ 369, 100 }}; /* "C14" -> "C21" */
    pos[ 75] = (r2_t){{ 348, 479 }}; /* "C76" -> "C71" */
    pos[ 76] = (r2_t){{ 358, 442 }}; /* "C77" -> "C67" */
    pos[ 77] = (r2_t){{ 360, 409 }}; /* "C78" -> "C61" */
    pos[ 78] = (r2_t){{ 359, 371 }}; /* "C79" -> "C54" */
    pos[ 79] = (r2_t){{ 348, 330 }}; /* "C80" -> "C31" */
    pos[ 81] = (r2_t){{ 358, 541 }}; /* "C82" -> "C74" */
    pos[ 82] = (r2_t){{ 378, 503 }}; /* "C83" -> "C70" */
    pos[ 83] = (r2_t){{ 389, 464 }}; /* "C84" -> "C66" */
    pos[ 84] = (r2_t){{ 392, 424 }}; /* "C85" -> "C60" */
    pos[ 85] = (r2_t){{ 388, 384 }}; /* "C86" -> "C53" */
    pos[ 86] = (r2_t){{ 378, 340 }}; /* "C87" -> "C37" */
    pos[ 87] = (r2_t){{ 397, 565 }}; /* "C88" -> "C73" */
    pos[ 88] = (r2_t){{ 421, 519 }}; /* "C89" -> "C69" */
    pos[ 89] = (r2_t){{ 433, 474 }}; /* "C90" -> "C65" */
    pos[ 90] = (r2_t){{ 432, 431 }}; /* "C91" -> "C59" */
    pos[ 91] = (r2_t){{ 425, 395 }}; /* "C92" -> "C52" */
    pos[ 92] = (r2_t){{ 409, 350 }}; /* "C93" -> "C42" */
    pos[ 93] = (r2_t){{ 468, 527 }}; /* "C94" -> "C68" */
    pos[ 94] = (r2_t){{ 474, 476 }}; /* "C95" -> "C64" */
    pos[ 95] = (r2_t){{ 475, 425 }}; /* "C96" -> "C58" */
    pos[ 96] = (r2_t){{ 460, 391 }}; /* "C97" -> "C51" */
    pos[ 97] = (r2_t){{ 436, 361 }}; /* "C98" -> "C47" */
    pos[ 98] = (r2_t){{ 523, 469 }}; /* "C99" -> "C63" */
    pos[ 99] = (r2_t){{ 518, 415 }}; /* "C100" -> "C57" */
    pos[100] = (r2_t){{ 498, 375 }}; /* "C101" -> "C50" */
    pos[101] = (r2_t){{ 474, 346 }}; /* "C102" -> "C46" */
    pos[102] = (r2_t){{ 442, 324 }}; /* "C103" -> "C41" */
    pos[103] = (r2_t){{ 412, 304 }}; /* "C104" -> "C36" */
    pos[104] = (r2_t){{ 378, 297 }}; /* "C105" -> "C30" */
    pos[105] = (r2_t){{ 339, 293 }}; /* "C106" -> "C7" */
    pos[106] = (r2_t){{ 558, 393 }}; /* "C107" -> "C56" */
    pos[107] = (r2_t){{ 506, 326 }}; /* "C108" -> "C45" */
    pos[108] = (r2_t){{ 471, 297 }}; /* "C109" -> "C40" ("T7"?) */
    pos[109] = (r2_t){{ 437, 279 }}; /* "C110" -> "C35" */
    pos[110] = (r2_t){{ 398, 269 }}; /* "C111" -> "C29" */
    pos[111] = (r2_t){{ 352, 267 }}; /* "C112" -> "C13" */
    pos[112] = (r2_t){{ 567, 328 }}; /* "C113" -> "C49" */
    pos[113] = (r2_t){{ 534, 287 }}; /* "C114" -> "C44" */
    pos[114] = (r2_t){{ 496, 261 }}; /* "C115" -> "C39" */
    pos[115] = (r2_t){{ 459, 250 }}; /* "C116" -> "C34" */
    pos[116] = (r2_t){{ 416, 241 }}; /* "C117" -> "C28" */
    pos[117] = (r2_t){{ 374, 242 }}; /* "C118" -> "C20" */
    pos[118] = (r2_t){{ 592, 234 }}; /* "C119" -> "C48" */
    pos[119] = (r2_t){{ 554, 237 }}; /* "C120" -> "C43" */
    pos[120] = (r2_t){{ 512, 215 }}; /* "C121" -> "C38" */
    pos[121] = (r2_t){{ 471, 207 }}; /* "C122" -> "C33" */
    pos[122] = (r2_t){{ 430, 208 }}; /* "C123" -> "C27" */
    pos[123] = (r2_t){{ 394, 218 }}; /* "C124" -> "C24" */
    pos[124] = (r2_t){{ 511, 155 }}; /* "C125" -> "C128" */
    pos[125] = (r2_t){{ 357,  40 }}; /* "C126" -> "C127" */
    /* Left half: */
    pos[  6] = sym(pos[106-1]); /* "C7" <- "C106" */
    pos[ 11] = sym(pos[  5-1]); /* "C12" <- "C5" */
    pos[ 12] = sym(pos[112-1]); /* "C13" <- "C112" */
    pos[ 17] = sym(pos[ 10-1]); /* "C18" <- "C10" */
    pos[ 18] = sym(pos[  4-1]); /* "C19" <- "C4" */
    pos[ 19] = sym(pos[118-1]); /* "C20" <- "C118" */
    pos[ 20] = sym(pos[ 14-1]); /* "C21" <- "C14" */
    pos[ 21] = sym(pos[  9-1]); /* "C22" <- "C9" */
    pos[ 22] = sym(pos[  3-1]); /* "C23" <- "C3" */
    pos[ 23] = sym(pos[124-1]); /* "C24" <- "C124" */
    pos[ 24] = sym(pos[  8-1]); /* "C25" <- "C8" */
    pos[ 25] = sym(pos[  2-1]); /* "C26" <- "C2" */
    pos[ 26] = sym(pos[123-1]); /* "C27" <- "C123" */
    pos[ 27] = sym(pos[117-1]); /* "C28" <- "C117" */
    pos[ 28] = sym(pos[111-1]); /* "C29" <- "C111" */
    pos[ 29] = sym(pos[105-1]); /* "C30" <- "C105" */
    pos[ 30] = sym(pos[ 80-1]); /* "C31" <- "C80" */
    pos[ 31] = sym(pos[  1-1]); /* "C32" <- "C1" */
    pos[ 32] = sym(pos[122-1]); /* "C33" <- "C122" */
    pos[ 33] = sym(pos[116-1]); /* "C34" <- "C116" */
    pos[ 34] = sym(pos[110-1]); /* "C35" <- "C110" */
    pos[ 35] = sym(pos[104-1]); /* "C36" <- "C104" */
    pos[ 36] = sym(pos[ 87-1]); /* "C37" <- "C87" */
    pos[ 37] = sym(pos[121-1]); /* "C38" <- "C121" */
    pos[ 38] = sym(pos[115-1]); /* "C39" <- "C115" */
    pos[ 39] = sym(pos[109-1]); /* "C40" <- "C109" ("T7"?) */
    pos[ 40] = sym(pos[103-1]); /* "C41" <- "C103" */
    pos[ 41] = sym(pos[ 93-1]); /* "C42" <- "C93" */
    pos[ 42] = sym(pos[120-1]); /* "C43" <- "C120" */
    pos[ 43] = sym(pos[114-1]); /* "C44" <- "C114" */
    pos[ 44] = sym(pos[108-1]); /* "C45" <- "C108" */
    pos[ 45] = sym(pos[102-1]); /* "C46" <- "C102" */
    pos[ 46] = sym(pos[ 98-1]); /* "C47" <- "C98" */
    pos[ 47] = sym(pos[119-1]); /* "C48" <- "C119" */
    pos[ 48] = sym(pos[113-1]); /* "C49" <- "C113" */
    pos[ 49] = sym(pos[101-1]); /* "C50" <- "C101" */
    pos[ 50] = sym(pos[ 97-1]); /* "C51" <- "C97" */
    pos[ 51] = sym(pos[ 92-1]); /* "C52" <- "C92" */
    pos[ 52] = sym(pos[ 86-1]); /* "C53" <- "C86" */
    pos[ 53] = sym(pos[ 79-1]); /* "C54" <- "C79" */
    pos[ 55] = sym(pos[107-1]); /* "C56" <- "C107" */
    pos[ 56] = sym(pos[100-1]); /* "C57" <- "C100" */
    pos[ 57] = sym(pos[ 96-1]); /* "C58" <- "C96" */
    pos[ 58] = sym(pos[ 91-1]); /* "C59" <- "C91" */
    pos[ 59] = sym(pos[ 85-1]); /* "C60" <- "C85" */
    pos[ 60] = sym(pos[ 78-1]); /* "C61" <- "C78" */
    pos[ 62] = sym(pos[ 99-1]); /* "C63" <- "C99" */
    pos[ 63] = sym(pos[ 95-1]); /* "C64" <- "C95" */
    pos[ 64] = sym(pos[ 90-1]); /* "C65" <- "C90" */
    pos[ 65] = sym(pos[ 84-1]); /* "C66" <- "C84" */
    pos[ 66] = sym(pos[ 77-1]); /* "C67" <- "C77" */
    pos[ 67] = sym(pos[ 94-1]); /* "C68" <- "C94" */
    pos[ 68] = sym(pos[ 89-1]); /* "C69" <- "C89" */
    pos[ 69] = sym(pos[ 83-1]); /* "C70" <- "C83" */
    pos[ 70] = sym(pos[ 76-1]); /* "C71" <- "C76" */
    pos[ 72] = sym(pos[ 88-1]); /* "C73" <- "C88" */
    pos[ 73] = sym(pos[ 82-1]); /* "C74" <- "C82" */
    pos[126] = sym(pos[126-1]); /* "C127" <- "C126" */
    pos[127] = sym(pos[125-1]); /* "C128" <- "C125" */
    assert(ne == 128);
    
    /* Rescale to unit disk: */
    int ie,j;
    for (ie = 0; ie < ne; ie++) 
      { for (j = 0; j < 2; j++) 
          { double z_raw = pos[ie].c[j];
            double z_rel = (z_raw - ctr_ref.c[j])/rad_ref.c[j];
            if (j == 1) { z_rel = -z_rel; }
            pos[ie].c[j] = z_rel;
          }
      }

    return pos;
  }
 
r3_t neuromat_eeg_geom_3D_from_2D(r2_t *p)
  {
    double R2 = r2_norm_sqr(p);
    double Den = 1 + R2;
    r3_t q = (r3_t){{ 2*p->c[0]/Den, 2*p->c[1]/Den, (R2 - 1)/Den }};
    /* Paranoia: */
    (void)r3_dir(&q, &q);
    return q;
  }

r2_t neuromat_eeg_geom_2D_from_3D(r3_t *p)
  {
    double Den = 1 + p->c[2];
    return (r2_t){{ p->c[0]/Den, p->c[1]/Den }};
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

void neuromat_eeg_geom_map_many_disk_to_ellipse(int np, r2_t p[], r2_t *ctr, r2_t *rad, r2_t q[])
  {
    int ip;
    for (ip = 0; ip < np; ip++) 
      { q[ip] = neuromat_eeg_geom_ellipse_from_disk(&(p[ip]), ctr, rad); }
  }

void neuromat_eeg_geom_map_many_ellipse_to_disk(int np, r2_t p[], r2_t *ctr, r2_t *rad, r2_t q[])
  {
    int ip;
    for (ip = 0; ip < np; ip++) 
      { q[ip] = neuromat_eeg_geom_disk_from_ellipse(&(p[ip]), ctr, rad); }
  }
