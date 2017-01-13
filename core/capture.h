#ifndef VISAR_CAPTURE
#define VISAR_CAPTURE

#include <opencv/cv.h>
#include <opencv/highgui.h>

CvSize initVideoCapture(const char *source);
CvSize initWebcamCapture(int webcamid, const CvSize& frameSize);
void closeCapture();
bool NextFrame(IplImage** frame);

#endif
