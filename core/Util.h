#ifndef VANILLA_UTIL
#define VANILLA_UTIL

#include <opencv/cv.h>
#include <opencv/highgui.h>

struct Util
{

static double rad2deg(double rad);
static double randdouble();
static double randdouble(double min, double max);
static void save_image(IplImage *im, const char *dir, const char *name, int frame_no);

};

#endif
