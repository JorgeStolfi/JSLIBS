#ifndef __LIBCAMFIREWIRE_SESSION_H__
#define __LIBCAMFIREWIRE_SESSION_H__

#include <dc1394/dc1394.h>
#include <stdio.h>
#include <bool.h>

struct camfirewire_session_t{
    dc1394_t* session_handler;
    dc1394camera_list_t* session_devices;
    
};


typedef struct camfirewire_session_t camfirewire_session_t;
typedef struct camfirewire_camera_t camfirewire_camera_t ;

camfirewire_session_t* camfirewire_session_init(void);
/*Initializes a dc1394 session and should be always the first step when handling such cameras,
it also lists all the devices available*/

void camfirewire_session_resease(camfirewire_session_t* cs);
/*Releases a dc1394 session. should be the last step of a software wich handles the camera*/

int camfirewire_session_get_num_devices(camfirewire_session_t* cs);
/*returns the number of available devices*/

#endif

