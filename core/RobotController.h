#ifndef VISAR_ROBOTCONTROLLER
#define VISAR_ROBOTCONTROLLER

#include "OverwriteQueue.h"

class RobotController
{
public:
  enum ID
  {
    DUMMY,
    ROVER,
    VISAR
  };

public:
  mutable bool move;
  OverwriteQueue_Ptr owqueue;
  double minSafeDistanceToObstacle;
  double minSafeAngleBetweenObstacles;

public:
  RobotController(const OverwriteQueue_Ptr& owqueue_, double minSafeDistanceToObstacle_, double minSafeAngleBetweenObstacles_);
  virtual ~RobotController() = 0;

public:
  virtual void control_robot(double angleToBestHeadingDirection, double distanceToObstacleInBestHeadingDirection, double angleBetweenBestObstacleFreeSegment) const = 0;

  virtual std::string construct_command(double turnAngle, double speed) const = 0;

  void write_command(const std::string& command) const;

  double calculate_speed(double distanceToObstacleInBestHeadingDirection, double angleBetweenBestObstacleFreeSegment) const;

  virtual double translate_turn_angle(double turnAngle) const = 0;
  virtual double translate_speed(double speed) const = 0;
};

#endif
