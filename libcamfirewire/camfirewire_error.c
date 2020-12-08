#include <camfirewire_error.h>
#include <assert.h>
#include <bool.h>

void camfirewire_error_camera_shutdown(dc1394camera_t *camera);
/*Shuts down the camera handler and exits the program with assertion error*/

dc1394error_t camfirewire_error_test(dc1394error_t err, char* error_text){
  DC1394_ERR_RTN(err,error_text);
  return err;
}

dc1394error_t camfirewire_error_test_camera_critical(dc1394error_t err, dc1394camera_t* camera_handler,char* error_text) {
  DC1394_ERR_CLN_RTN(err,camfirewire_error_camera_shutdown(camera_handler),error_text);
  return err;
}

void camfirewire_error_camera_shutdown(dc1394camera_t *camera)
{
    dc1394_video_set_transmission(camera, DC1394_OFF);
    dc1394_capture_stop(camera);
    dc1394_camera_free(camera);
    fprintf(stderr,"Camera Panic! : Shutting down program");
    assert(FALSE);
}
