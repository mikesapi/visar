#ifndef VISAR_PIPELINE
#define VISAR_PIPELINE

#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv/cxcore.h>

#include "structures.h"

class SegmentationPipeline
{
public:
  IplImage *source_img;
  IplImage *gray_img;
  IplImage *gbs_img;
  IplImage *contour;

  CvRect SafeRegion;

  Features features;
  Statistics stats;
  Model statsmodel;

  bool debug;
  bool dynamic;
  int frameNo;
  bool refresh;
  std::string saveOutputDir;

  int imageWidth;
  int imageHeight;
  int windowWidth;
  int windowHeight;

public:
  SegmentationPipeline(int width, int height, bool dynamicFlag_, bool debugFlag_, std::string saveOutputDir = "");
  ~SegmentationPipeline();

public:
  void display(IplImage *im, const char *name, int x, int y);
  void process_frame(IplImage *input, IplImage *output);
  void reinit_stats();

private:
  CvRect get_safe_region(int imageWidth, int imageHeight);
};

#endif
