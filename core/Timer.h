#ifndef VISAR_TIMER
#define VISAR_TIMER

#include <string>

class Timer
{
private:
  double startTime;
  double endTime;

public:
  Timer();
  void end(const std::string& tag);
};

#endif
