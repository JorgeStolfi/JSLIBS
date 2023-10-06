The data structures for mesh topology: half-edgem quad-edge, oct-edge.
Uses low-order bits of address to indicate sym, flip and rot.

  hedge.h, hedge.c
  
    The half-edge data structure for orientable meshes. 

  quad.h, quad.c
  
    The quad-edge data structure for orientable maps.
    
  oct.h, oct.c 
  
    The oct-edge (full quad-edge) data structure for non-orientable maps.

  haf_draw.h, haf_draw.c 
  
    Procedures to draw the mesh and pointers from an half-edge data structure.