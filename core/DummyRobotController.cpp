#include "DummyRobotController.h"

#include <iostream>

DummyRobotController::DummyRobotController(const OverwriteQueue_Ptr& owqueue_, double minSafeDistanceToObstacle_, double minSafeAngleBetweenObstacles_)
: RobotController(owqueue_, minSafeDistanceToObstacle_, minSafeAngleBetweenObstacles_)
{
  std::cout << "DUMMY\n";
}

DummyRobotController::~DummyRobotController()
{
  std::cout << "DUMMY ";
}

void DummyRobotController::control_robot(double angleToBestHeadingDirection, double distanceToObstacleInBestHeadingDirection, double angleBetweenBestObstacleFreeSegment) const
{
  double speed = calculate_speed(distanceToObstacleInBestHeadingDirection, angleBetweenBestObstacleFreeSegment);
  write_command(construct_command(angleToBestHeadingDirection, speed));
}

std::string DummyRobotController::construct_command(double turnAngle, double speed) const
{
  char command[200];
  sprintf(command,"Angle=%0.2f, Speed=%0.2f", translate_turn_angle(turnAngle), translate_speed(speed));
  return std::string(command);
}

double DummyRobotController::translate_turn_angle(double turnAngle) const
{
  return turnAngle;
}

double DummyRobotController::translate_speed(double speed) const
{
  return speed;
}
