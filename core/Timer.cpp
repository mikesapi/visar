#include "Timer.h"

#include <opencv/cv.h>

Timer::Timer()
{
  startTime = (double)cvGetTickCount();
}

void Timer::end(const std::string& tag)
{
  endTime = (double)cvGetTickCount();

  double diffTime = endTime - startTime;
  
  double diffTimeMs = diffTime/((double)cvGetTickFrequency()*1000.0);
  std::cout << "Time for " << tag << ": " << diffTimeMs << "ms" << std::flush << std::endl;
}
