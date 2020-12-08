#ifndef archdraw_ic_H
#define archdraw_ic_H

/* Common defs for IC architectural drawings */
/* Last edited on 2009-09-18 01:45:25 by stolfi */

#define _GNU_SOURCE
#include <frgb.h>
#include <r3.h>
#include <vec.h>
#include <epswr.h>

#include <archdraw.h>

typedef enum 
  { 
    adrw_ic_space_type_ADM, /* Administracao/Of.manut. */
    adrw_ic_space_type_AUD, /* Auditório. */
    adrw_ic_space_type_BAN, /* Banheiros. */
    adrw_ic_space_type_DEP, /* Depósito/Limpeza. */
    adrw_ic_space_type_DOC, /* Docentes/Posdocs/Visitantes. */
    adrw_ic_space_type_HAL, /* Corredores/Saguões/Escadas/Elevadores. */
    adrw_ic_space_type_INF, /* Analistas/Of.eletronica. */
    adrw_ic_space_type_POS, /* Laboratórios/Escr.Doutorandos. */
    adrw_ic_space_type_NEX, /* Inexistente. */
    adrw_ic_space_type_PIL, /* Pilares. */
    adrw_ic_space_type_REU, /* Reunião/Estar/Lazer. */
    adrw_ic_space_type_SIT, /* Assentos. */
    adrw_ic_space_type_SRV, /* Servidoras/Impressoras/Xerox/Rede. */
    adrw_ic_space_type_ARR, /* Setas. */
    adrw_ic_space_type_ETC  /* Etc. (must be last). */
  } adrw_ic_space_type_t;
  /* Code for office usage type. */
  
#define adrw_ic_space_type_MAX (adrw_ic_space_type_ETC)
  /* The largest {space_type_t} value. */

char **adrw_ic_define_type_tags(void);
  /* Returns an array, indexed by {adrw_ic_space_type_t}, with a three-letter 
    abbrev of each usage type. */

frgb_t *adrw_ic_define_type_fill_colors(void);
  /* Returns an array, indexed by {adrw_ic_space_type_t}, with the fill color to use
    for each space type. */

adrw_unit_style_t **adrw_ic_define_type_styles(void);
  /* Returns an array, indexed by {adrw_ic_space_type_t}, with the style to use
    when plotting each space type. */

bool_t *adrw_ic_define_printable_types(void);
  /* Returns an array, indexed by {adrw_ic_space_type_t}, with TRUE iff
    the usage type should appear in the building units printout. */

bool_t *adrw_ic_define_movable_types(void);
  /* Returns an array, indexed by {adrw_ic_space_type_t}, with TRUE iff
    the usage type should be included in the total area to be moved to IC-4a. */

int *adrw_ic_define_color_key_types(void);
  /* Returns an array with the types which should be displayed in 
    the color key, in the proper order.  The end of the array 
    is inicated by a {-1} value. */

void adrw_ic_plot_histogram_bars
  ( char *prefix, 
    adrw_building_t *B, 
    bool_t select[], 
    char *type_tag[], 
    adrw_unit_style_t *style[]
  );
  /* Plots the histogram bars for modules and areas in building {B}.
    Each bar will a separate EPS file, called "{prefix}_{ttag}_{vtag}"
    where {ttag} is the tag of the space type, and {vtag} is
    "mods" for number of modules, or "area" for area. */


#endif
