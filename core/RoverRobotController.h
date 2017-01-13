#ifndef VISAR_ROVERROBOTCONTROLLER
#define VISAR_ROVERROBOTCONTROLLER

#include "RobotController.h"

class RoverRobotController : public RobotController
{
public:
  RoverRobotController(const OverwriteQueue_Ptr& owqueue_, double minSafeDistanceToObstacle_, double minSafeAngleBetweenObstacles_);
  virtual ~RoverRobotController();

public:
  /* Override */
  virtual void control_robot(double angleToBestHeadingDirection, double distanceToObstacleInBestHeadingDirection, double angleBetweenBestObstacleFreeSegment) const;

  /* Override */
  virtual std::string construct_command(double turnAngle, double speed) const;
  virtual double translate_turn_angle(double turnAngle) const;
  virtual double translate_speed(double speed) const;
};

#endif
