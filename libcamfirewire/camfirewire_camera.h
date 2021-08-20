/* Last edited on 2021-07-24 16:29:25 by jstolfi */

#ifndef camfirewire_camera_H
#define camfirewire_camera_H

#define _GNU_SOURCE_
#include <stdio.h>

#include <camfirewire_session.h> 
#include <dc1394/dc1394.h>

struct camfirewire_camera_t{
  dc1394camera_t* camera_handler;

  dc1394video_modes_t camera_modes;
  dc1394framerates_t camera_framerates;
};

camfirewire_camera_t* camfirewire_camera_init(camfirewire_session_t* cs,int n);
  /*Get the {n} camera from the session. If there is no camera available
    with id {n}, returns null, also, alocates the buffer space of
    {buffer_size} for such camera */

void camfirewire_list_camera_features(camfirewire_camera_t* cc, FILE* arq);
  /* Prints in ${arq} the camera capabilities */

void camfirewire_camera_reset(camfirewire_camera_t* cc, bool_t use_default_settings);
  /*Reset the camera handlers. If {use_default_settings} is TRUE resets 
    its configurations to default */

void camfirewire_camera_release(camfirewire_camera_t* cc);
  /*Close camera handlers and release used frame buffers and structures*/

#endif
