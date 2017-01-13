/*
 *  Copyright (C) <2014>  <Michael Sapienza>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "core/Util.h"
#include "core/DummyRobotController.h"
#include "core/RoverRobotController.h"
#include "core/Params.h"
#include "core/ProgramInputParser.h"
#include "core/SegmentationPipeline.h"
#include "core/init_structures.h"
#include "core/img_proc_fcns.h"
#include "core/capture.h"
#include "core/OverwriteQueue.h"
#include "core/CommunicationLoop.h"
#include "core/Timer.h"

int loop_through(int j, int keypress_id, int size_of_collection = 99)
{
  //control test image sequence
  if(keypress_id=='w' || keypress_id == 1048695)
  {
    j--;
    j = j < 0 ? size_of_collection : j;
  }
  else if(keypress_id == 'e' || keypress_id == 1048677)
  {
    j++;
    j = j > size_of_collection ? 0 : j;
  }

  return j;
}

void help()
{
  std::vector<std::string> optionIds = {"-h","-d","-r","-i","-c","-imageDataset"};
  std::vector<std::string> optionDescriptions =
  {"help",
   "set debug flag",
   "select specific robot [name: dummy, rover]",
   "input video path [path: foo/bar.avi]",
   "capture from camera [optional id: {0, 1}]",
   "use the static traversability image dataset [optional path: foo/bar]"};

  std::cout << "Command line options:\n";
  for(size_t i = 0; i < optionIds.size(); ++i)
  {
    std::cout << optionIds[i] << " -> " << optionDescriptions[i] << '\n';
  }
}

std::string get_option_helper(const ProgramInputParser& pip, const std::string& optionId)
{
  std::string option;
  if(pip.option_id_exists(optionId))
  {
    option = pip.get_option(optionId);
    if(option.empty())
    {
      help();
      throw std::runtime_error("option not found for option id: " + optionId);
    }
    else
    {
      std::cout << "Found option " << option << " for id: " << optionId<< std::endl;
    }
  }

  return option;
}

int set_params_from_input_args(int argc, char **argv, Params *p)
{
  ProgramInputParser pip(argc, argv);

  if(pip.option_id_exists("-h"))
  {
    help();
    return 0;
  }

  p->debug = pip.option_id_exists("-d");

  std::string robotOptId("-r");
  if(pip.option_id_exists(robotOptId))
  {
    p->robot_name = get_option_helper(pip, robotOptId);
  }
  p->init_specific_robot();

  std::string inputVideoOptId("-i");
  if(pip.option_id_exists(inputVideoOptId))
  {
    p->captureSource = get_option_helper(pip, inputVideoOptId);
  }

  std::string cameraOptId("-c");
  if(pip.option_id_exists(cameraOptId))
  {
    p->captureType = Params::CT_WEBCAM; p->webcamid = 0;
    std::string option = pip.get_option(cameraOptId);
    if(!option.empty())
    {
      p->webcamid = std::stoi(option);
    }
  }

  std::string imageDataOptId("-imageDataset");
  if(pip.option_id_exists(imageDataOptId))
  {
    p->captureType = Params::CT_IMAGESEQ;
    p->captureSource = "static_traversability_dataset";

    std::string option = pip.get_option(imageDataOptId);
    if(!option.empty())
    {
      p->captureSource = option;
    }
  }

  return 1;
}

//#define USECOMMS
//#define CRANFIELDDATASETMODE

typedef std::shared_ptr<RobotController> RobotController_Ptr;

Params p;
Params *pPtr = &p;

int main(int argc, char** argv)
{
  if(argc > 1)
  {
    if(!set_params_from_input_args(argc, argv, pPtr)) return 0;
  }

  bool dynamicFlag(true);
  if(p.captureType == Params::CT_IMAGESEQ) dynamicFlag = false;
 
  int key=0; //holds user-input keystroke

  OverwriteQueue_Ptr overwriteQueue(nullptr);
#ifdef USECOMMS
  net::Address sendSocketAddress(127,0,0,1,8888);

  if(p.robot_id == RobotController::ID::ROVER)
  {
    sendSocketAddress = net::Address(192,168,2,6,8888);
  }

  // Create socket if communicating with robot via UDP
  net::Socket socket;
  net::Socket *socketPtr = &socket;
  net::SocketUtil::InitSocketUdp(socketPtr, sendSocketAddress.GetPort());

  // Create a buffer to be able to store the control messages.
  // The buffer is there so that the communication loop can send messages to a robot at a constant rate,
  // whilst receiving new updates from the vision at a variable rate.
  overwriteQueue.reset(new OverwriteQueue<std::string>(5));

  // The communication loop takes care of sending the control messages in the buffer.
  CommunicationLoop commloop(sendSocketAddress);
  std::thread communicationThread(&CommunicationLoop::Start, &commloop, overwriteQueue, socketPtr);
#endif

  // The robot takes care of translating the output of the computer vision to each respective robot.
  RobotController_Ptr robotController;
  if(p.robot_id == RobotController::ID::DUMMY)
  {
    robotController.reset(new DummyRobotController(overwriteQueue, p.bot.min_dist_to_obst, p.bot.min_ang_between_obst));
  }
  else if(p.robot_id == RobotController::ID::ROVER)
  {
    robotController.reset(new RoverRobotController(overwriteQueue, p.bot.min_dist_to_obst, p.bot.min_ang_between_obst));
  }
  else
  {
    throw std::runtime_error("No robot selected\n");
  }

  Boundary Bound;
  Boundary *BoundPtr = &Bound;

   //LOOPS for generating results on datasets
#ifdef CRANFIELDDATASETMODE
  char cranfieldDatasetVideoNames[8][30] =
  {{"CloudyDry"},{"CloudyWet"},{"CloudyMuddy"},{"SunnyWet"},{"ComplexScene"},{"Fence"},{"Shadows"},{"Snow"}};

  for(int id=0; id<8; ++id)
  {

  std::string datadir = "/media/mikesapi/DATASETDISK/ms-workspace/traversabilitydetection";
  //std::string datadir = "/home/mikesapi/media/data/traversabilitydetection";
  std::string inputdir =  datadir + '/' + "input";
  p.captureSource = inputdir + '/' + cranfieldDatasetVideoNames[id] + ".avi";

  const std::string resultsdir(datadir + '/' + "results");
  std::string saveOutputDir = resultsdir + '/' + cranfieldDatasetVideoNames[id];

  SegmentationPipeline pipeline(160, 120, dynamicFlag, p.debug, saveOutputDir);
#else

  SegmentationPipeline pipeline(160, 120, dynamicFlag, p.debug);
#endif

  // T represents the log posterior threshold on a opencv trackbar, 0 => 300
  int T=0;
  p.log_post_thres_zero_position = T*10 + 300;

  CvSize frameSize = cvSize(cvRound(p.proc_W), cvRound(p.proc_H));

  CvSize captureSize;
  if(p.captureType == Params::CT_WEBCAM)
  {
    captureSize = initWebcamCapture(p.webcamid, frameSize);
  }
  else if(p.captureType == Params::CT_VIDEO)
  {
    // capture video from file or network stream
    captureSize = initVideoCapture(p.captureSource.c_str());
  }

  IplImage *input_image = cvCreateImageHeader(captureSize, 8, 3);

  init_images_img_proc(frameSize);
  init_boundary(frameSize, BoundPtr);

#ifndef CRANFIELDDATASETMODE
  printf( "\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n"
  "This software runs a vision-based autonomous guidance algorithm\n"
  "for a mobile robot equipped with a low-quality monocular camera.\n"
  "The vision system allows a mobile robot to autonomously guide itself\n"
  "past static or dynamic obstacles in both indoor or outdoor natural environments\n"
  "in a real-time, reactive manner.\n\n"
  "To begin, make sure all parameters are correct and press 'Enter'\n\n"
  "To exit, click inside one of the display windows,\n"
  "then press the 'q' key\n\n\n"
  "***WARNING***\n"
  "--Prolonged use of this software may cause injuries--"
  "\n*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n");

  cvWaitKey(0);
#endif

  //MAIN PROCESSING LOOP
  while(1)
  {
      //--START-->CAPTURE FRAME//
      if(p.captureType == Params::CT_IMAGESEQ)
      {
        static int j = 0;
        static int change = 0;

        if(change != j)
        {
          pipeline.reinit_stats();
          change=j;
        }

        j = loop_through(j, key); //control test image sequence

        char buffer [250];
        sprintf(buffer,"%s/source_%d.bmp", p.captureSource.c_str(), j);

        input_image = cvLoadImage(buffer);
        if(!input_image) throw std::runtime_error("Could not load image: " + std::string(buffer));
      }
      else
      {
        if(!NextFrame(&input_image)) break; //captures a video frame from webcam or video file
      }

    Timer processingLoop;

      // Traversability segmentation pipeline.
    Timer timeSegmentation;
      pipeline.process_frame(input_image, Bound.Bimg);
    timeSegmentation.end("timeSegmentation");

#ifndef CRANFIELDDATASETMODE
    Timer timeBoundaryInterpretation;
      // Extracting depth and steering angle from the segmentation.
      ExtractBoundary(frameSize, BoundPtr);
      CalculateDistances(frameSize, BoundPtr, p.camera, p.bot);
      InterpretDepthArray(frameSize, BoundPtr, p.bot);
    timeBoundaryInterpretation.end("BoundaryInterpretation");

    Timer timeRobotController;
      // Interpreting the information extracted from image sensor for a specific robot.
      robotController->control_robot(Util::rad2deg(Bound.angleToBestHeadingDirection), Bound.distanceToObstacleInBestHeadingDirection, Util::rad2deg(Bound.angleBetweenBestObstacleFreeSegment));
    timeRobotController.end("RobotController");
#endif

      // Termination of processing loop.
    processingLoop.end("Entire Processing Loop");
    std::cout << '\n';

      if(!p.debug) 	key = cvWaitKey(10);
      else 		key = cvWaitKey(p.debug_delay_ms);

      if(key == 'q' || key == 'x' || key == 1048689 || key == 1048603)
      break;
      p.refresh=0;
  }

#ifdef USECOMMS
  communicationThread.join();
#endif

#ifdef CRANFIELDDATASETMODE
  } //end for id
#endif

  // Release resources.
  closeCapture();
  release_images_img_proc();
  release_boundary(BoundPtr);
  cvDestroyAllWindows();

#ifdef USECOMMS
  if(socketPtr->IsOpen()) net::ShutdownSockets();
#endif

  return 0;
}//end main


