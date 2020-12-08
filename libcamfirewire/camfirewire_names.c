#include <camfirewire_names.h>

char* video_mode_names[] = { 
    "DC1394_VIDEO_MODE_160x120_YUV444",
    "DC1394_VIDEO_MODE_320x240_YUV422",
    "DC1394_VIDEO_MODE_640x480_YUV411",
    "DC1394_VIDEO_MODE_640x480_YUV422",
    "DC1394_VIDEO_MODE_640x480_RGB8",
    "DC1394_VIDEO_MODE_640x480_MONO8",
    "DC1394_VIDEO_MODE_640x480_MONO16",
    "DC1394_VIDEO_MODE_800x600_YUV422",
    "DC1394_VIDEO_MODE_800x600_RGB8",
    "DC1394_VIDEO_MODE_800x600_MONO8",
    "DC1394_VIDEO_MODE_1024x768_YUV422",
    "DC1394_VIDEO_MODE_1024x768_RGB8",
    "DC1394_VIDEO_MODE_1024x768_MONO8",
    "DC1394_VIDEO_MODE_800x600_MONO16",
    "DC1394_VIDEO_MODE_1024x768_MONO16",
    "DC1394_VIDEO_MODE_1280x960_YUV422",
    "DC1394_VIDEO_MODE_1280x960_RGB8",
    "DC1394_VIDEO_MODE_1280x960_MONO8",
    "DC1394_VIDEO_MODE_1600x1200_YUV422",
    "DC1394_VIDEO_MODE_1600x1200_RGB8",
    "DC1394_VIDEO_MODE_1600x1200_MONO8",
    "DC1394_VIDEO_MODE_1280x960_MONO16",
    "DC1394_VIDEO_MODE_1600x1200_MONO16",
    "DC1394_VIDEO_MODE_EXIF",
    "DC1394_VIDEO_MODE_FORMAT7_0",
    "DC1394_VIDEO_MODE_FORMAT7_1",
    "DC1394_VIDEO_MODE_FORMAT7_2",
    "DC1394_VIDEO_MODE_FORMAT7_3",
    "DC1394_VIDEO_MODE_FORMAT7_4",
    "DC1394_VIDEO_MODE_FORMAT7_5",
    "DC1394_VIDEO_MODE_FORMAT7_6",
    "DC1394_VIDEO_MODE_FORMAT7_7",
    "UNKNOW_VIDEO_MODE"
};
    
char* framerate_names[] = {
    "DC1394_FRAMERATE_1_875",
    "DC1394_FRAMERATE_3_75",
    "DC1394_FRAMERATE_7_5",
    "DC1394_FRAMERATE_15",
    "DC1394_FRAMERATE_30",
    "DC1394_FRAMERATE_60",
    "DC1394_FRAMERATE_120",
    "DC1394_FRAMERATE_240",
    "UNKNOW_FRAMERATE"
};

char* color_coding_names[] = {
    "DC1394_COLOR_CODING_MONO8",
    "DC1394_COLOR_CODING_YUV411",
    "DC1394_COLOR_CODING_YUV422",
    "DC1394_COLOR_CODING_YUV444",
    "DC1394_COLOR_CODING_RGB8",
    "DC1394_COLOR_CODING_MONO16",
    "DC1394_COLOR_CODING_RGB16",
    "DC1394_COLOR_CODING_MONO16S",
    "DC1394_COLOR_CODING_RGB16S",
    "DC1394_COLOR_CODING_RAW8",
    "DC1394_COLOR_CODING_RAW16",
    "UNKNOW_COLOR_CODING"
};

char* iso_speed_names[] = {
    "DC1394_ISO_SPEED_100",
    "DC1394_ISO_SPEED_200",
    "DC1394_ISO_SPEED_400",
    "DC1394_ISO_SPEED_800",
    "DC1394_ISO_SPEED_1600",
    "DC1394_ISO_SPEED_3200",
    "UNKNOW_ISO_SPEED"
};

char* operation_mode_names[] = {
  "DC1394_OPERATION_MODE_LEGACY",
  "DC1394_OPERATION_MODE_1394B",
  "UNKNOW_OPERATION_MODE"
};

char* camfirewire_property_name(uint32_t value, uint32_t min,uint32_t num,char* array_name[]);
/*Generic function to treat this name handling properly*/

char* camfirewire_property_name(uint32_t value, uint32_t min,uint32_t num,char* array_name[]){
  int ind = value - min;
  if( (ind < 0) || (ind >= num) ){
    ind = num;
  }
  return array_name[ind];
}

char* camfirewire_video_mode_name(dc1394video_mode_t video_mode){
    return camfirewire_property_name(video_mode,DC1394_VIDEO_MODE_MIN,DC1394_VIDEO_MODE_NUM,video_mode_names);
}

char* camfirewire_framerate_name(dc1394framerate_t framerate){
  return camfirewire_property_name(framerate,DC1394_FRAMERATE_MIN,DC1394_FRAMERATE_NUM,framerate_names);
}

char* camfirewire_color_coding_name(dc1394color_coding_t color_coding){
  return camfirewire_property_name(color_coding,DC1394_COLOR_CODING_MIN,DC1394_COLOR_CODING_NUM,color_coding_names);
}

char* camfirewire_iso_speed_name(dc1394speed_t iso_speed){
  return camfirewire_property_name(iso_speed,DC1394_ISO_SPEED_MIN,DC1394_ISO_SPEED_NUM,iso_speed_names);
}

char* camfirewire_operation_mode_name(dc1394operation_mode_t mode){
  return camfirewire_property_name(mode,DC1394_OPERATION_MODE_MIN,DC1394_OPERATION_MODE_NUM,operation_mode_names);
}

