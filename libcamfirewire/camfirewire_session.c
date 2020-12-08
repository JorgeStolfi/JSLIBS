/* See {camfirewire_session.h}. */
/* Last edited on 2020-12-07 23:55:09 by jstolfi */

#include <camfirewire_session.h>
#include <camfirewire_names.h>
#include <camfirewire_error.h>
#include <assert.h>
#include <affirm.h>
#include <stdlib.h>

camfirewire_session_t* camfirewire_session_init(void){
  
  camfirewire_session_t* cs = (camfirewire_session_t*)malloc(sizeof(camfirewire_session_t));
  dc1394error_t err;
  cs->session_handler = dc1394_new ();
  if( cs->session_handler  != NULL) {
    fprintf(stderr,"camfirewire_init_session : Could not init dc1394 session !");
    
  }
  assert( cs->session_handler  != NULL);
  err = dc1394_camera_enumerate (cs->session_handler, &(cs->session_devices));
  camfirewire_error_test(err,"camfirewire_init_session : Failed to enumerate cameras");
  /*Prints out the Device ID */
  int  num_devices = cs->session_devices->num;
  int i;
  fprintf(stderr,"-------------------------\n");
  if(num_devices == 0){
    fprintf(stderr,"No devices available. \n");
  }
  for(i = 0; i < num_devices; i++){
    fprintf(stderr,"Device [%d]\n GUID %lx UNIT %d\n",i,cs->session_devices->ids[i].guid,cs->session_devices->ids[i].unit);
  }
  fprintf(stderr,"-------------------------\n");
  return cs;
}

void camfirewire_session_release(camfirewire_session_t* cs){
  
  dc1394_camera_free_list(cs->session_devices);
  dc1394_free(cs->session_handler);
  
  free(cs);
}

int camfirewire_session_get_num_devices(camfirewire_session_t* cs){
  return cs->session_devices->num;
}


