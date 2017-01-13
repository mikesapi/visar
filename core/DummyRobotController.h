#ifndef VISAR_DUMMYROBOTCONTROLLER
#define VISAR_DUMMYROBOTCONTROLLER

#include "RobotController.h"

class DummyRobotController : public RobotController
{
public:
  DummyRobotController(const OverwriteQueue_Ptr& owqueue_, double minSafeDistanceToObstacle, double minSafeAngleBetweenObstacles);
  virtual ~DummyRobotController();

public:
  /* Override */
  virtual void control_robot(double angleToBestHeadingDirection, double distanceToObstacleInBestHeadingDirection, double angleBetweenBestObstacleFreeSegment) const;

  /* Override */
  virtual std::string construct_command(double turnAngle, double speed) const;
  virtual double translate_turn_angle(double turnAngle) const;
  virtual double translate_speed(double speed) const;
};

#endif
