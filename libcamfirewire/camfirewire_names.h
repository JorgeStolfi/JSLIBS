#ifndef __LIBCAMFIREWIRE_NAMES_H__
#define __LIBCAMFIREWIRE_NAMES_H__

#include <dc1394/dc1394.h>

char* camfirewire_video_mode_name(dc1394video_mode_t video_mode);
/*Given a video mode it will give its name*/
char* camfirewire_framerate_name(dc1394framerate_t framerate);
/*Given a framerate it will give its name*/
char* camfirewire_color_coding_name(dc1394color_coding_t color_coding);
/*Given a color coding it will give its name*/
char* camfirewire_iso_speed_name(dc1394speed_t iso_speed);
/*Given a iso speed it will give its name*/
char* camfirewire_operation_mode_name(dc1394operation_mode_t mode);
/*Given a operation mode it will give its name*/
#endif