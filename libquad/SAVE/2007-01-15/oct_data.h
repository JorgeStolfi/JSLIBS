/* Last edited on 2007-01-15 16:19:22 by stolfi */

/* Aborted plan to have separate data fields for each arc. */
/* Better to have a single ID number for each edge record + tumble bits. */

/* CLIENT DATA FIELDS

  For each edge of the map, the structure includes four 64-bit data
  fields that can be set independently by the client. These fields are
  assumed to be associated with the extremities of that edge and of
  its dual edge. Reading and writing are always done relative to an
  arc {e}, by the procedures {oct_{OP}_{X}data(e)} where {OP}
  specifies the operation ('{get}' or '{set}') and {X} specifies the
  data field ('{o}', '{d}', '{r}' or '{l}') --- relative to {e}. */

typedef uint64_t oct_data_t; 
  /* A data value associated with an arc. */

oct_data_t oct_get_odata(oct_arc_t e);
oct_data_t oct_set_odata(oct_arc_t e);
  /* Read/write the data associated to the origin extremity of {e}. */

oct_data_t oct_get_ddata(oct_arc_t e);
oct_data_t oct_set_ddata(oct_arc_t e);
  /* Read/write the data associated to the destination extremity of {e}. */

oct_data_t oct_get_ldata(oct_arc_t e);
oct_data_t oct_set_ldata(oct_arc_t e);
  /* Read/write the data associated to the left bank of {e}. */

oct_data_t oct_get_rdata(oct_arc_t e);
oct_data_t oct_set_rdata(oct_arc_t e);
  /* Read/write the data associated to the right bank of {e}. */

/* Bit tweaking of data fields

  Note that the same data field can be read out in eight different
  ways: {get_odata(e),get_odata(flip(e)),get_ldata(rot(e))}, etc.
  
  If the same arc {e} is used for writing and reading a field, the
  the data is preserved in full: after {set_rdata(e,X)}, the call
  {get_rdata(e)} will return {X}. However, if a data field is written
  relative to an arc {e} and read out relative to another arc {f}, the
  {get} procedure may complement either of the last two bits of the
  data value, depending on the relation between {e} and {f}.
  
  Specifically, bit 0 (the units bit) of the result will be
  complemented if {e} and {f} have opposite circular orientations; and
  bit 1 will be complemented if {e} and {f} differ in primal/dual
  character.  Thus, for example, after {set_odata(e,X)}, we will
  have
  
    {get_odata(e) == get_ddata(sym(e)) == X} 
    {get_odata(flip(e)) == get_ddata(vlip(e)) == X^1}
    {get_ldata(rot(e)) == get_rdata(tor(e)) == X^2}
    {get_ldata(flip(tor(e)) == get_rdata(flip(rot(e))) == X^3}
    
  Similar relations hold for the other operators; for example, after
  {set_ldata(e,X)}, one will get
  
    {get_ldata(e) == get_rdata(sym(e)) == X},
    {get_ldata(vlip(e)) == get_rdata(flip(e)) == X^1},
    
  etc. */

