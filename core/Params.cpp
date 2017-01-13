#include "RobotController.h"
#include "Params.h"

#include <opencv/highgui.h>

//Params constructor initialises parameters to default values
const std::string Params::valid_robot_names[3] = {"dummy","rover","visar"};

Params::Params()
: captureType(CT_VIDEO),
  captureSource("data/eng_stat_obst.avi"),
  debug(false),
  debug_delay_ms(10),
  log_post_thres_zero_position(300),
  proc_H(120),
  proc_W(160),
  refresh(true), // refresh window names and positions.
  robot_name("dummy"),
  robot_id(0),
  webcamid(0)
{}

//initialise the parameters of a specific robot
void Params::init_specific_robot()
{
  int listsize = sizeof valid_robot_names / sizeof valid_robot_names[0];
  for(int i=0; i<listsize; ++i)
  {
    if(valid_robot_names[i].compare(robot_name)==0)
    {
      robot_id=i;
      break;
    }
  }

  if((robot_id < 0) || (robot_id > 2)) throw std::runtime_error("Invalid robot id");
  else
  {
    camera.cx = proc_W/2;
    camera.cy = proc_H/2;
    camera.height_from_ground = 400.0; //mm
    camera.fx = 3.6; //mm
    camera.fy = 3.6; //mm
    camera.pix_size_x = 0.0014*(2592/proc_W); //mm
    camera.pix_size_y = 0.0014*(1944/proc_H); //mm
    camera.tilt_angle = 0.0;
    camera.pan_angle = 0.0;

    bot.min_dist_to_obst = 4000; //mm
    bot.min_ang_between_obst = 15;
    bot.display_scale_factor = 0.08;
    bot.max_depth_threshold = 1e6;
  }

  if(robot_id == RobotController::ID::ROVER)
  {
    // Params common to all bots.
    captureSource = "tcpclientsrc host=192.168.2.127 port=5000  ! gdpdepay !  rtph264depay ! ffdec_h264 ! ffmpegcolorspace ! appsink";

    camera.height_from_ground = 240.0; //mm

    bot.min_dist_to_obst = 2000; //mm
  }
  else if(robot_id == RobotController::ID::VISAR)
  {
    throw std::runtime_error("WIP..");
  }

}
