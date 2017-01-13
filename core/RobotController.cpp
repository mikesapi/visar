#include "RobotController.h"

#include <iostream>

RobotController::RobotController(const OverwriteQueue_Ptr& owqueue_, double minSafeDistanceToObstacle_, double minSafeAngleBetweenObstacles_)
: move(false),
  owqueue(owqueue_),
  minSafeDistanceToObstacle(minSafeDistanceToObstacle_),
  minSafeAngleBetweenObstacles(minSafeAngleBetweenObstacles_)
{
  std::cout << "HELLO FROM ";
}

RobotController::~RobotController()
{
  std::cout << "SAYS BYE\n";
}

void RobotController::write_command(const std::string& command) const
{
  if(owqueue) owqueue->Reset(command);
  else std::cout << "Command: " << command << '\n';
}

double RobotController::calculate_speed(double distanceToObstacleInBestHeadingDirection, double angleBetweenBestObstacleFreeSegment) const
{
  double a = distanceToObstacleInBestHeadingDirection;
  double b = minSafeDistanceToObstacle;

  double speed = 0.0;
  if(angleBetweenBestObstacleFreeSegment > minSafeAngleBetweenObstacles)
  {
    if(a >= 3*b)
    {
      speed = 0.15;
    }
    else if((a < 3*b) && (a >= 1.5*b))
    {
      speed = 0.06;
    }
    else if((a < 1.5*b) && (a >= 0.8*b))
    {
      speed = 0.02;
    }
  }

  return speed;
}
