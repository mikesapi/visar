#include <cstdio> 

#include "capture.h"
#include "init_structures.h"

CvCapture *capture = 0;

CvSize initVideoCapture(const char *source)
{
  printf("Capturing from %s\n", source);
  capture = cvCreateFileCapture(source);

  int frames = (int) cvGetCaptureProperty(
          capture,
          CV_CAP_PROP_FRAME_COUNT
          );

  printf("no of frames in video = %d\n\n", frames);

  if(!frames)
  {
    closeCapture();
    throw std::runtime_error("Failed to initialise camera capture");
  }

  return cvSize((int)cvGetCaptureProperty(capture, CV_CAP_PROP_FRAME_WIDTH),
                (int)cvGetCaptureProperty(capture, CV_CAP_PROP_FRAME_HEIGHT));
}

CvSize initWebcamCapture(int webcamid, const CvSize& frameSize)
{
  capture = cvCaptureFromCAM(webcamid);
  cvSetCaptureProperty(capture, CV_CAP_PROP_FRAME_WIDTH, frameSize.width);
  cvSetCaptureProperty(capture, CV_CAP_PROP_FRAME_HEIGHT, frameSize.height);

  if(!capture)
  {
    closeCapture();
    std::runtime_error("Failed to initialise camera capture");
  }

  return cvSize((int)cvGetCaptureProperty(capture, CV_CAP_PROP_FRAME_WIDTH),
                (int)cvGetCaptureProperty(capture, CV_CAP_PROP_FRAME_HEIGHT));
}

void closeCapture()
{
  cvReleaseCapture(&capture);
}

bool NextFrame(IplImage **frame)
{
  *frame = cvQueryFrame( capture );

  if(!*frame)
  {
    fprintf(stderr, "failed to get a video frame\n");
    return false;
  }

  return true;
}
