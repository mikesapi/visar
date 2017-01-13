#ifndef VISAR_PARAMS
#define VISAR_PARAMS

#include "botCalib.h"
#include "camCalib.h"

#include <opencv/cv.h>

class Params
{
public:
  enum CaptureType
  {
    CT_WEBCAM,
    CT_VIDEO,
    CT_IMAGESEQ,
  };

public:
    //capture
    CaptureType captureType;
    std::string captureSource;

    bool debug;
    int debug_delay_ms;

    int log_post_thres_zero_position;

    //processing resolution
    double proc_H;
    double proc_W;

    bool refresh; //refresh the window names and positions

    std::string robot_name;
    const static std::string valid_robot_names[3];
    int robot_id;

    int webcamid;

    CamCalib camera;
    BotCalib bot;

public:
    Params();
    void init_specific_robot();
};

#endif
