#ifndef VANILLA_BOTCALIB
#define VANILLA_BITCALIB

struct botCalib {
  
    //structure collecting the parameters related to the robot configuration
    double min_dist_to_obst;
    double min_ang_between_obst;
    double max_depth_threshold;
    double display_scale_factor;
};
typedef struct botCalib BotCalib;

#endif
