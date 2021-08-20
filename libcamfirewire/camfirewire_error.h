#ifndef camfirewire_error_H
#define camfirewire_error_H

/* Last edited on 2021-07-24 16:30:17 by jstolfi */

#include <dc1394/dc1394.h>

dc1394error_t camfirewire_error_test(dc1394error_t err, char* error_text);
  /*A wrapper for this defined return error function of dc1394*/

dc1394error_t camfirewire_error_test_camera_critical
  ( dc1394error_t err, 
    dc1394camera_t* camera_handler,
    char* error_text
  ) ;
  /*A wrapper for erro function with support to a shutdown callback*/

#endif
