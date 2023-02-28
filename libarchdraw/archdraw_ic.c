/* See {archdraw_ic.h} */
/* Last edited on 2023-02-20 18:20:30 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>

#include <affirm.h>
#include <frgb.h>
/* #include <epswr.h> */

/* #include <archdraw.h> */
#include <archdraw_ic.h>

/* INTERNAL PROTOTYPES */

void adrw_ic_plot_histogram_bar_pair
  ( char *dir,
    char *prefix,
    int32_t type,
    char *type_tag,
    adrw_unit_style_t *style,
    double mods,
    double area
   );
   /* Writes two EPS files with histogram bars for the number of modules {mods}
     and the total area {area} for the space type with index {type} and tag {type_tag},
     using the color of the specified {style}. 
     
     The files will be called "{prefix}_{type_tag}_mods.eps" and  "{prefix}_{type_tag}_area.eps",
     respectively. */

frgb_t *adrw_ic_define_type_fill_colors(void)
  {
    int32_t ntypes = adrw_ic_space_type_MAX + 1;
    frgb_t *color = (frgb_t *)notnull(malloc(ntypes*sizeof(frgb_t)), "no mem");
    
    color[adrw_ic_space_type_ADM] = (frgb_t){{ 0.000f, 0.600f, 0.550f }}; /* Y = 0.4150 */
    color[adrw_ic_space_type_AUD] = (frgb_t){{ 0.650f, 0.350f, 1.000f }}; /* Y = 0.5050 */
    color[adrw_ic_space_type_DOC] = (frgb_t){{ 0.000f, 0.520f, 1.000f }}; /* Y = 0.4120 */
    color[adrw_ic_space_type_INF] = (frgb_t){{ 0.150f, 0.680f, 0.000f }}; /* Y = 0.4530 */
    
    color[adrw_ic_space_type_BAN] = (frgb_t){{ 0.790f, 0.910f, 0.890f }}; /* Y = 0.8720 */
    color[adrw_ic_space_type_HAL] = (frgb_t){{ 0.900f, 0.850f, 0.750f }}; /* Y = 0.8550 */
    color[adrw_ic_space_type_REU] = (frgb_t){{ 0.820f, 0.820f, 0.820f }}; /* Y = 0.8200 */
    color[adrw_ic_space_type_DEP] = (frgb_t){{ 0.930f, 0.780f, 0.530f }}; /* Y = 0.8000 */
    color[adrw_ic_space_type_POS] = (frgb_t){{ 0.850f, 0.850f, 0.650f }}; /* Y = 0.8300 */
    color[adrw_ic_space_type_SRV] = (frgb_t){{ 1.000f, 0.700f, 0.850f }}; /* Y = 0.8050 */

    color[adrw_ic_space_type_ETC] = (frgb_t){{ 0.900f, 0.900f, 0.900f }}; /* Y = 0.9000 */
    color[adrw_ic_space_type_ARR] = (frgb_t){{ 1.000f, 0.300f, 0.000f }}; /* Y = 0.4800 */

    color[adrw_ic_space_type_NEX] = (frgb_t){{ 0.200f, 0.200f, 0.200f }}; /* Y = 0.2000 */
    color[adrw_ic_space_type_PIL] = (frgb_t){{ 0.000f, 0.300f, 0.100f }}; /* Y = 0.1900 */
    color[adrw_ic_space_type_SIT] = (frgb_t){{ 0.500f, 0.100f, 0.700f }}; /* Y = 0.2800 */

    return color;
  }

char **adrw_ic_define_type_tags(void)
  { int32_t ntypes = adrw_ic_space_type_MAX + 1;
    char **tag = (char **)notnull(malloc(ntypes*sizeof(char *)), "no mem");
    tag[adrw_ic_space_type_ADM] = "ADM";
    tag[adrw_ic_space_type_AUD] = "AUD";
    tag[adrw_ic_space_type_BAN] = "BAN";
    tag[adrw_ic_space_type_DEP] = "DEP";
    tag[adrw_ic_space_type_DOC] = "DOC";
    tag[adrw_ic_space_type_ETC] = "ETC";
    tag[adrw_ic_space_type_HAL] = "HAL";
    tag[adrw_ic_space_type_INF] = "INF";
    tag[adrw_ic_space_type_POS] = "POS";
    tag[adrw_ic_space_type_NEX] = "NEX";
    tag[adrw_ic_space_type_PIL] = "PIL";
    tag[adrw_ic_space_type_REU] = "REU";
    tag[adrw_ic_space_type_SIT] = "SIT";
    tag[adrw_ic_space_type_SRV] = "SRV";
    tag[adrw_ic_space_type_ARR] = "ARR";
    return tag;
  } 

bool_t *adrw_ic_define_printable_types(void)
  { int32_t ntypes = adrw_ic_space_type_MAX + 1;
    bool_t *prt = (bool_t *)notnull(malloc(ntypes*sizeof(bool_t)), "no mem");
    int32_t type;
    for (type = 0; type < ntypes; type++) { prt[type] = TRUE; }
    prt[adrw_ic_space_type_ETC] = FALSE;
    prt[adrw_ic_space_type_HAL] = FALSE;
    prt[adrw_ic_space_type_NEX] = FALSE;
    prt[adrw_ic_space_type_PIL] = FALSE;
    prt[adrw_ic_space_type_SIT] = FALSE;
    prt[adrw_ic_space_type_ARR] = FALSE;
    return prt;
  } 

bool_t *adrw_ic_define_movable_types(void)
  { int32_t ntypes = adrw_ic_space_type_MAX + 1;
    bool_t *mov = (bool_t *)notnull(malloc(ntypes*sizeof(bool_t)), "no mem");
    int32_t type;
    for (type = 0; type < ntypes; type++) { mov[type] = FALSE; }
    mov[adrw_ic_space_type_ADM] = TRUE;
    mov[adrw_ic_space_type_AUD] = TRUE;
    mov[adrw_ic_space_type_DOC] = TRUE;
    mov[adrw_ic_space_type_INF] = TRUE;
    return mov;
  } 

int32_t *adrw_ic_define_color_key_types(void)
  { int32_t nkeys = 10;
    int32_t *tp = (int32_t *)notnull(malloc((nkeys+1)*sizeof(int32_t)), "no mem");
    tp[ 0] = adrw_ic_space_type_DOC;
    tp[ 1] = adrw_ic_space_type_ADM;
    tp[ 2] = adrw_ic_space_type_INF;
    tp[ 3] = adrw_ic_space_type_AUD;
    tp[ 4] = adrw_ic_space_type_REU;
    tp[ 5] = adrw_ic_space_type_SRV;
    tp[ 6] = adrw_ic_space_type_DEP;
    tp[ 7] = adrw_ic_space_type_POS;
    tp[ 8] = adrw_ic_space_type_BAN;
    tp[ 9] = adrw_ic_space_type_HAL;
    tp[10] = -1;
    /* tp[  ] = adrw_ic_space_type_ETC; */
    /* tp[  ] = adrw_ic_space_type_NEX; */
    /* tp[  ] = adrw_ic_space_type_PIL; */
    /* tp[  ] = adrw_ic_space_type_SIT; */
    /* tp[  ] = adrw_ic_space_type_ARR; */
    return tp;
  } 

adrw_unit_style_t **adrw_ic_define_type_styles(void)
  {
    int32_t ntypes = adrw_ic_space_type_MAX + 1;
    adrw_unit_style_t **style = (adrw_unit_style_t **)notnull(malloc(ntypes*sizeof(adrw_unit_style_t *)), "no mem");
    frgb_t *color = adrw_ic_define_type_fill_colors();
    frgb_t black = (frgb_t){{0,0,0}};
    int32_t i;
    for (i = 0; i <= adrw_ic_space_type_MAX; i++)
      { style[i] = adrw_make_unit_style(&(color[i]),&black,0.10,NULL,0.0); }
    return style;
  }

void adrw_ic_plot_histogram_bars
  ( char *dir,
    char *prefix,
    adrw_building_t *B,
    bool_t select[], 
    char *type_tag[], 
    adrw_unit_style_t *style[]
  )
  { fprintf(stderr, "--- plotting histogram bars ---\n");
    int32_t ntypes = adrw_ic_space_type_MAX+1;
    double modules[ntypes];
    double area[ntypes];
    adrw_compute_building_stats(B, ntypes, NULL, modules, area, NULL);
    double tot_modules = 0, tot_area = 0;
    int32_t type;
    for (type = 0; type < ntypes; type++)
      { if ((select != NULL) && (! select[type])) { continue; }
        fprintf(stderr, "type %2d modules = %6.1f area = %6.1f\n", type, modules[type], area[type]);
        
        adrw_ic_plot_histogram_bar_pair(dir, prefix, type, type_tag[type], style[type], modules[type], area[type]);

        tot_modules += modules[type];
        tot_area += area[type];
      }
      
    type = adrw_ic_space_type_SIT;
    adrw_ic_plot_histogram_bar_pair(dir, prefix, type, "TUT", style[type], tot_modules, tot_area);

    type = adrw_ic_space_type_ETC;
    adrw_ic_plot_histogram_bar_pair(dir, prefix, type, "TCO", style[type], modules[type], area[type]);
  }

void adrw_ic_plot_histogram_bar_pair
  ( char *dir,
    char *prefix,
    int32_t type,
    char *type_tag,
    adrw_unit_style_t *style,
    double mods,
    double area
   )
   { adrw_plot_histogram_bar(dir, prefix, type_tag, "mods", type, style, mods, 100.0);
     adrw_plot_histogram_bar(dir, prefix, type_tag, "area", type, style, area/10000, 2400.0);
   }
