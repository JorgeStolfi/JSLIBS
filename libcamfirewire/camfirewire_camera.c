/* See {camfirewire_camera.h}. */
/* Last edited on 2020-12-08 13:57:44 by jstolfi */

#define _GNU_SOURCE_
#include <stdio.h>
#include <stdlib.h>

#include <camfirewire_names.h>
#include <camfirewire_error.h>

#include <camfirewire_camera.h>

camfirewire_camera_t* camfirewire_camera_init(camfirewire_session_t* cs,int n){
  
  int num_devices = cs->session_devices->num;
  if((n < 0) || (n >= num_devices) ){
    fprintf(stderr,"camfirewire_camera_init : Invalid camera number - %d\n",n);
    return NULL;
  }
  
  camfirewire_camera_t* cc = (camfirewire_camera_t*) malloc(sizeof(camfirewire_camera_t));
  cc->camera_handler =  dc1394_camera_new (cs->session_handler, cs->session_devices->ids[n].guid);
  if( cc->camera_handler == NULL) {
    free(cc);
    fprintf(stderr,"camfirewire_camera_init : Failure accesing camera number - %d\n",n);
    dc1394_log_error("Failed to initialize camera with guid %llx", cs->session_devices->ids[n].guid);
    return NULL;
  }
  
  dc1394error_t err = dc1394_video_get_supported_modes(cc->camera_handler, &(cc->camera_modes));
  camfirewire_error_test(err ,"camfirewire_camera_init : Could not get list of modes");

  uint32_t selected_mode = cc->camera_modes.modes[0];
  
  /*Now we will set up the parameters to full compatibility mode*/
  
  err = dc1394_video_set_iso_speed(cc->camera_handler, DC1394_ISO_SPEED_400);
  camfirewire_error_test_camera_critical(err,cc->camera_handler,"camfirewire_camera_init : Could not set iso speed");
  
  err = dc1394_video_set_mode(cc->camera_handler, selected_mode);
  camfirewire_error_test_camera_critical(err,cc->camera_handler,"camfirewire_camera_init : Could not set video mode\n");

  err = dc1394_video_set_framerate(cc->camera_handler, DC1394_FRAMERATE_7_5);
  camfirewire_error_test_camera_critical(err,cc->camera_handler,"camfirewire_camera_init : Could not set framerate\n");

  err = dc1394_capture_setup(cc->camera_handler,4, DC1394_CAPTURE_FLAGS_DEFAULT);
  camfirewire_error_test_camera_critical(err,cc->camera_handler,"camfirewire_camera_init : Could not setup camera-\n");
  
  return cc;
}



void camfirewire_camera_reset(camfirewire_camera_t *cc,bool_t use_default_settings)
{
    dc1394error_t err;
    
    dc1394_video_set_transmission(cc->camera_handler, DC1394_OFF);
    dc1394_capture_stop(cc->camera_handler);
    
    if(use_default_settings){
	uint32_t selected_mode = cc->camera_modes.modes[0];
	 
	err = dc1394_video_set_iso_speed(cc->camera_handler, DC1394_ISO_SPEED_400);
	camfirewire_error_test_camera_critical(err,cc->camera_handler,"camfirewire_camera_reset : Could not set iso speed");
    
	err = dc1394_video_set_mode(cc->camera_handler, selected_mode);
	camfirewire_error_test_camera_critical(err,cc->camera_handler,"camfirewire_camera_reset: Could not set video mode\n");
    
	err = dc1394_video_set_framerate(cc->camera_handler, DC1394_FRAMERATE_7_5);
	camfirewire_error_test_camera_critical(err,cc->camera_handler,"camfirewire_camera_reset: Could not set framerate\n");

	err = dc1394_capture_setup(cc->camera_handler,4, DC1394_CAPTURE_FLAGS_DEFAULT);
	camfirewire_error_test_camera_critical(err,cc->camera_handler,"camfirewire_camera_reset: Could not setup camera-\n");
    }
    
    free(cc);
    
}

void camfirewire_camera_release(camfirewire_camera_t *cc)
{

    dc1394_video_set_transmission(cc->camera_handler, DC1394_OFF);
    dc1394_capture_stop(cc->camera_handler);
    dc1394_camera_free(cc->camera_handler);
    
    free(cc);
    
}



void camfirewire_list_camera_features(camfirewire_camera_t* cc, FILE* arq){
     fprintf(arq,"------CAM-INFORMATIONS------\n");
     dc1394_camera_print_info(cc->camera_handler,arq);
//   fprintf(arq,"GUID: %llx\n",cc->camera_handler->guid);
//   fprintf(arq,"UNIT: %d\n",cc->camera_handler->unit);
//   fprintf(arq,"VENDOR: %s\n",cc->camera_handler->vendor);
//   fprintf(arq,"MODEL: %s\n",cc->camera_handler->model);
    fprintf(arq,"------CURRENT-SETTINGS-------\n");
    dc1394error_t err;
    
    dc1394framerate_t framerate;
    float fps;
    err = dc1394_video_get_framerate(cc->camera_handler, &framerate);
    if (err != 0) { fprintf(stderr, "{dc1394_video_get_framerate} err = %d\n", err); }
    err = dc1394_framerate_as_float(framerate,&fps);
    if (err != 0) { fprintf(stderr, "{dc1394_framerate_as_float} err = %d\n", err); }
    
    dc1394video_mode_t video_mode;
    uint32_t h,w;
    err = dc1394_video_get_mode(cc->camera_handler,&video_mode);
    if (err != 0) { fprintf(stderr, "{dc1394_video_get_mode} err = %d\n", err); }
    dc1394_get_image_size_from_video_mode(cc->camera_handler, video_mode, &w, &h);
    
    dc1394operation_mode_t mode;
    err = dc1394_video_get_operation_mode(cc->camera_handler, &mode);
    if (err != 0) { fprintf(stderr, "{dc1394_video_get_operation_mode} err = %d\n", err); }
    
    dc1394speed_t speed;
    err = dc1394_video_get_iso_speed(cc->camera_handler,&speed);
    if (err != 0) { fprintf(stderr, "{dc1394_video_get_iso_speed} err = %d\n", err); }
    
    fprintf(arq,"ISO SPEED: %s \n",camfirewire_iso_speed_name(speed));
    fprintf(arq,"VIDEO MODE: %s (%dx%d)\n",camfirewire_video_mode_name(video_mode),w,h);
    fprintf(arq,"FRAMERATE: %3.2f \n",fps);
    fprintf(arq,"OPERATION MODE: %s \n",camfirewire_operation_mode_name(mode));
    fprintf(arq,"------VIDEO-MODES------------\n");
    fprintf(arq,"Supported Modes: \n\n");
    int i;
    for(i = 0; i < cc->camera_modes.num; i++){
      
      video_mode = cc->camera_modes.modes[i];
      dc1394_get_image_size_from_video_mode(cc->camera_handler, video_mode, &w, &h);
      fprintf(arq,"%s - %dx%d ",camfirewire_video_mode_name(video_mode),w,h);
      
      if(!dc1394_is_video_mode_scalable(video_mode)){
	fprintf(arq,"fps (");
	dc1394framerates_t framerate_list;
	err = dc1394_video_get_supported_framerates(cc->camera_handler, video_mode,&framerate_list);
        if (err != 0) { fprintf(stderr, "{dc1394_video_get_supported_framerates} err = %d\n", err); }
	for(int j= 0; j < framerate_list.num; j++){
	  framerate = framerate_list.framerates[j];
	  float fps;
	  err = dc1394_framerate_as_float(framerate,&fps);
          if (err != 0) { fprintf(stderr, "{dc1394_framerate_as_float} err = %d\n", err); }
	  fprintf(arq," %3.2f ",fps);
	}
	fprintf(arq,")");
      }else{
	fprintf(arq," Scalable");
      }
      fprintf(arq,"\n");
      
    }

  fprintf(arq,"------FEATURES---------------\n");
  
  dc1394featureset_t featureset;
  dc1394_feature_get_all(cc->camera_handler,&featureset);
  dc1394_feature_print_all(&featureset, arq);

  fprintf(arq,"-----------------------------\n");
}


