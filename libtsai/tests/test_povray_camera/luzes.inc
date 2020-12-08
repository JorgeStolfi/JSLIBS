// Luzes do LIV-V
// Last edited on 2007-07-22 17:46:02 by stolfi

#macro luzes(lux)
  #if (lux > 0)
    
    light_source {
      < 1200, 1200, 3300 >, 
      color rgb 0.20 * (lux) * <1.0,1.0,1.0>
    }
    light_source {
      < 1200, 3000, 3300 >, 
      color rgb 0.20 * (lux) * <1.0,1.0,1.0>
    }
    light_source {
      < 1200, 4800, 3300 >, 
      color rgb 0.20 * (lux) * <1.0,1.0,1.0>
    }

    light_source {
      < 3000, 1200, 3300 >, 
      color rgb 0.22 * (lux) * <1.0,1.0,1.0>
    }
    light_source {
      < 3000, 3000, 3300 >, 
      color rgb 0.22 * (lux) * <1.0,1.0,1.0>
    }
    light_source {
      < 3000, 4800, 3300 >, 
      color rgb 0.22 * (lux) * <1.0,1.0,1.0>
    }
    
    light_source {
      < 4800, 1200, 3300 >, 
      color rgb 0.22 * (lux) * <1.0,1.0,1.0>
    }
    light_source {
      < 4800, 3000, 3300 >, 
      color rgb 0.22 * (lux) * <1.0,1.0,1.0>
    }
    light_source {
      < 4800, 4800, 3300 >, 
      color rgb 0.22 * (lux) * <1.0,1.0,1.0>
    }
  #end
#end
