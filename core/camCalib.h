#ifndef VANILLA_CAMCALIB
#define VANILLA_CAMCALIB

struct camCalib
{
    //structure collecting the parameters related to the robot camera configuration
    double height_from_ground;
    double fy, fx;
    double tilt_angle, pan_angle;
    double cx, cy;
    double pix_size_x, pix_size_y;
};
typedef struct camCalib CamCalib;

#endif

