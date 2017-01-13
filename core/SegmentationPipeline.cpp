#include "SegmentationPipeline.h"

#include "init_structures.h"
#include "img_proc_fcns.h"
#include "Util.h"

SegmentationPipeline::SegmentationPipeline(int width, int height, bool dynamic_, bool debug_, std::string saveOutputDir_)
: debug(debug_),
  dynamic(dynamic_),
  frameNo(-1), //Katramados videos start from -1
  refresh(true),
  saveOutputDir(saveOutputDir_),
  imageWidth(width),
  imageHeight(height)
{
  windowWidth = 1.02 * width;
  windowHeight = 1.3 * height + 20;

  CvSize FrameSize = cvSize(width, height);

  source_img = cvCreateImage(FrameSize, 8, 3);
  gray_img   = cvCreateImage(FrameSize, 8, 1);
  gbs_img    = cvCreateImage(FrameSize, 8, 1);
  contour = cvCreateImage(FrameSize, 8, 3);

  // Initialise the safe region binary mask image
  SafeRegion = get_safe_region(imageWidth, imageHeight);


  init_model(FrameSize, SafeRegion, &statsmodel);

  if(debug) display(statsmodel.safeRegionMask, "safe region 2", 0 , 1);

  init_stats(FrameSize, &stats, 1);
  init_features(FrameSize, &features);
}

SegmentationPipeline::~SegmentationPipeline()
{
  cvReleaseImage(&source_img);
  cvReleaseImage(&gray_img);
  cvReleaseImage(&gbs_img);
  cvReleaseImage(&contour);

  release_features(&features);
  release_stats(&stats);
  release_model(&statsmodel);
}

void SegmentationPipeline::reinit_stats()
{
  const bool allocate(false);
  init_stats(cvSize(imageWidth, imageHeight), &stats, allocate);
}

void SegmentationPipeline::process_frame(IplImage *input, IplImage *output)
{
  cvResize(input, source_img, CV_INTER_AREA);
  if(debug) display(source_img, "source", 0 , 0);

  //Convert Color to gray RGB[A]->Gray: Y<-0.299*R + 0.587*G + 0.114*B
  cvCvtColor(source_img, gray_img, CV_BGR2GRAY);
  if(debug) display(gray_img, "gray", 1, 0);

  cvCvtColor(source_img, features.hsv, 	CV_BGR2HSV);
  cvCvtColor(source_img, features.YCrCb,CV_BGR2YCrCb);
  cvCvtColor(source_img, features.lab, 	CV_BGR2Lab);
  //Divides a multi-channel array into separate single-channel arrays
  cvCvtPixToPlane(features.hsv, features.hue, features.sat, 0, 0);
  cvCvtPixToPlane(features.YCrCb, 0, features.Cr, features.Cb, 0);
  cvCvtPixToPlane(features.lab, 0, features.a, 0, 0);

  //D = (A + B + 2*C)/4
  //illumination invariant color channel combination
  combine_channels(features.Cr, features.Cb, features.a, features.iic);
  cvNormalize(features.iic, features.iic, 0, 255, CV_MINMAX);
  convertGray2LBP(gray_img, features.lbp); //Convert to LBP
  getEdgeMagOri(features.sat, features.mag, features.ang32); //Extract edge magnitudes and orientation

  if(debug)
  {
    display(features.hue, "Hue", 2, 0);
    display(features.sat, "Sat", 3, 0);
    display(features.iic, "IIC", 4, 0);
    display(features.mag, "mag", 5, 0);
    display(features.ang32, "ang32", 6, 0);
    display(features.lbp, "lbp", 7, 0);
  }

  GraphBasedSegmentation(gbs_img, gray_img);
  SuperPixelStats(gbs_img, gray_img, &stats);  

  if(debug) display(stats.prior_img, "prior", 0, 2);

  if(dynamic)
  {
    UpdatePrior(gbs_img, &stats, &features);
  }

  GetModel(gray_img, &features, &statsmodel, dynamic);

  FeatureAnalysis(&features, &statsmodel, &stats, gbs_img, dynamic);
  ProbAnalysis2(&features, &stats, gbs_img);

  if(debug) display(features.post1, "posterior", 0, 4);

  UpdateParams(features.bin_class_result, &stats, &features, dynamic);

  // Not sure what this is for.
  //cvMerge(features.bin_class_result, features.bin_class_result, features.bin_class_result, NULL, result_img);
  
  // Post process the binary segmentation result.
  cvDilate(features.bin_class_result, features.bin_class_result, NULL, 1);
  cvErode(features.bin_class_result, features.bin_class_result, NULL, 3);
  cvSmooth(features.bin_class_result, features.bin_class_result, CV_MEDIAN, 3, 3);
  FindObstacleBoundary(features.bin_class_result);

  if(debug)
  {
    display(features.post1, "posterior", 0, 4);
    display(features.bin_class_result, "binary segmentation result", 1, 4);
  }

  Contours(features.bin_class_result);
  DrawContours(contour, CV_RGB(255,0,0), SafeRegion);
  cvAddWeighted(contour, 0.5, source_img, 1, 0, source_img);
  display(source_img, "result display", 0, 0);
  ReleaseContours();

  cvCopy(features.bin_class_result, output);
  refresh = false;
  std::cout << "Frame: " << ++frameNo << '\n';

  if(!saveOutputDir.empty() && (frameNo % 25 == 0) && (frameNo > 0))
  {
    Util::save_image(source_img, saveOutputDir.c_str(), "source", frameNo);
    Util::save_image(features.bin_class_result, saveOutputDir.c_str(), "binaryImage", frameNo);
  }
}

void SegmentationPipeline::display(IplImage *im, const char *name, int x, int y)
{
  if(refresh)
  {
    cvNamedWindow(name);
    cvResizeWindow(name, imageWidth, imageHeight);
    cvMoveWindow(name, int(round(x*windowWidth)),
                       int(round(y*windowHeight)));
  }

  cvShowImage(name, im);
}

CvRect SegmentationPipeline::get_safe_region(int imageWidth, int imageHeight)
{
  return cvRect(cvRound(imageWidth/3.0f),
                imageHeight - cvRound(imageHeight/7.0f),
                cvRound(imageWidth/3.0f),
                cvRound(imageHeight/8.0f)
                );
}
