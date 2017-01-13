#include "RoverRobotController.h"

#include <iostream>
#include <cmath>

RoverRobotController::RoverRobotController(const OverwriteQueue_Ptr& owqueue_, double minSafeDistanceToObstacle_, double minSafeAngleBetweenObstacles_)
: RobotController(owqueue_, minSafeDistanceToObstacle_, minSafeAngleBetweenObstacles_)
{
  std::cout << "ROVER\n";
}

RoverRobotController::~RoverRobotController()
{
  std::cout << "ROVER ";
}

void RoverRobotController::control_robot(double angleToBestHeadingDirection, double distanceToObstacleInBestHeadingDirection, double angleBetweenBestObstacleFreeSegment) const
{
  static double turn(0.0);
  static double speed(0.0);

  // Temporal smoothing
  turn = 0.3*turn + 0.7*angleToBestHeadingDirection;

  // Speed range -1/+1
  // Turn angle range -90/+90 degrees

  if(move)
  {
    speed = 0.3*speed + 0.7*calculate_speed(distanceToObstacleInBestHeadingDirection, angleBetweenBestObstacleFreeSegment);

    write_command(construct_command(turn, speed));

    if(turn < 0.01 && speed < 0.01) move = false;
  }
  else
  {
    if(angleBetweenBestObstacleFreeSegment < minSafeAngleBetweenObstacles)
    {
      std::cout << "Rotation routine not yet implemented\n";
      //rotation_routine();
    }
    else if(angleBetweenBestObstacleFreeSegment >= minSafeAngleBetweenObstacles)
    {
      move = true;
    }
  }
}

std::string RoverRobotController::construct_command(double turnAngle, double speed) const
{
  char command[200];
  sprintf(command,"%ds%dt", int(translate_turn_angle(turnAngle)), int(translate_speed(speed)));
  return std::string(command);
}

double RoverRobotController::translate_turn_angle(double turnAngle) const
{
  return round(50 - turnAngle);
}

double RoverRobotController::translate_speed(double speed) const
{
  //scaled = ((mtx - MIN)./RANGE).*(MAXVAL-MINVAL) + MINVAL;
  int offset = 0;
  int scale = 5;
  return round(((( speed * scale + offset)-(-1))/2)*(75-25)+25);
}
