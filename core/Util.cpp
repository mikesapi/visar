#include "Util.h"

#include <cmath>
#include <cstdlib>

double Util::rad2deg(double rad)
{
  return rad * (180.0 / M_PI);
}

double Util::randdouble()
{
  return rand() / (double(RAND_MAX) + 1);
}

double Util::randdouble(double min, double max)
{
  if (min>max)
  {
    return randdouble()*(min-max)+max;
  }
  else
  {
    return randdouble()*(max-min)+min;
  }
}

void Util::save_image(IplImage *im, const char *dir, const char *name, int frame_no)
{
  char path[250];
  sprintf(path,"%s/%s_%d.png", dir, name, frame_no);
  cv::Mat image(im);
  std::cout << "\nWriting image: " << path << '\n' << std::flush;
  cv::imwrite(path, image, std::vector<int>(0));
}

