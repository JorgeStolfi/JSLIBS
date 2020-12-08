/* See {fn1_functions.h}. */
/* Last edited on 2007-01-04 02:59:35 by stolfi */

#include <fn1_functions.h>

#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <affirm.h>
#include <bool.h>

#include <fn1_f1.h>
#include <fn1_f2.h>
#include <fn1_f3.h>
#include <fn1_f4.h>
#include <fn1_f5.h>

#include <fn1_g1.h>
#include <fn1_g2.h>
#include <fn1_g3.h>
#include <fn1_g4.h>
#include <fn1_g5.h>
#include <fn1_g6.h>
#include <fn1_g7.h>
#include <fn1_g8.h>
#include <fn1_g10.h>
#include <fn1_g11.h>
#include <fn1_g12.h>
#include <fn1_g14.h>
#include <fn1_gasqrt.h>
#include <fn1_gbadia.h>
#include <fn1_gbadmul.h>
#include <fn1_gdiv.h>
#include <fn1_gexp.h>
#include <fn1_giaboom.h>
#include <fn1_ginv.h>
#include <fn1_glog4.h>
#include <fn1_gsin3.h>
#include <fn1_gsqr.h>
#include <fn1_gsqrt.h>

#include <string.h>

fn1_data_t fn1_from_tag(char *tag)
  {
    if (strcmp(tag, "f1"                ) == 0) { return fn1_f1_get_data();       }
    if (strcmp(tag, "f2"                ) == 0) { return fn1_f2_get_data();       }
    if (strcmp(tag, "f3"                ) == 0) { return fn1_f3_get_data();       }
    if (strcmp(tag, "f4"                ) == 0) { return fn1_f4_get_data();       }
    if (strcmp(tag, "f5"                ) == 0) { return fn1_f5_get_data();       }     

    if (strcmp(tag, "g1"                ) == 0) { return fn1_g1_get_data();       }
    if (strcmp(tag, "g2"                ) == 0) { return fn1_g2_get_data();       }
    if (strcmp(tag, "g3"                ) == 0) { return fn1_g3_get_data();       }
    if (strcmp(tag, "g4"                ) == 0) { return fn1_g4_get_data();       }
    if (strcmp(tag, "g5"                ) == 0) { return fn1_g5_get_data();       }
    if (strcmp(tag, "g6"                ) == 0) { return fn1_g6_get_data();       }
    if (strcmp(tag, "g7"                ) == 0) { return fn1_g7_get_data();       }
    if (strcmp(tag, "g8"                ) == 0) { return fn1_g8_get_data();       }
    if (strcmp(tag, "g10"               ) == 0) { return fn1_g10_get_data();      }
    if (strcmp(tag, "g11"               ) == 0) { return fn1_g11_get_data();      }
    if (strcmp(tag, "g12"               ) == 0) { return fn1_g12_get_data();      }
    if (strcmp(tag, "g14"               ) == 0) { return fn1_g14_get_data();      }
    if (strcmp(tag, "gasqrt"            ) == 0) { return fn1_gasqrt_get_data();   }
    if (strcmp(tag, "gbadia"            ) == 0) { return fn1_gbadia_get_data();   }
    if (strcmp(tag, "gbadmul"           ) == 0) { return fn1_gbadmul_get_data();  }
    if (strcmp(tag, "gdiv"              ) == 0) { return fn1_gdiv_get_data();     }
    if (strcmp(tag, "gexp"              ) == 0) { return fn1_gexp_get_data();     }
    if (strcmp(tag, "giaboom"           ) == 0) { return fn1_giaboom_get_data();  }
    if (strcmp(tag, "ginv"              ) == 0) { return fn1_ginv_get_data();     }
    if (strcmp(tag, "glog4"             ) == 0) { return fn1_glog4_get_data();    }
    if (strcmp(tag, "gsin3"             ) == 0) { return fn1_gsin3_get_data();    }
    if (strcmp(tag, "gsqr"              ) == 0) { return fn1_gsqr_get_data();     }
    if (strcmp(tag, "gsqrt"             ) == 0) { return fn1_gsqrt_get_data();    }
    demand(FALSE, "bad function tag");
    return fn1_g1_get_data(); /* To pacify the compiler. */
  }
