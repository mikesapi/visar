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

#ifndef VISAR_COMMUNICATIONLOOP
#define VISAR_COMMUNICATIONLOOP

#include "OverwriteQueue.h"

#include <iostream>

#include <thread>

#include "../comms/Net.h"

class CommunicationLoop
{
private:
  int comm_time_;
  bool stop_;
  int watch_dog_;
  int wait_limit_;

  net::Address sendSocketAddress_;
  
public:
  CommunicationLoop(const net::Address& sendSocketAddress)
  : comm_time_(50),
    stop_(false),
    watch_dog_(5),
    wait_limit_(5),
    sendSocketAddress_(sendSocketAddress)
    {}

  void Send(const OverwriteQueue_Ptr& overwriteQueue, net::Socket *socPtr)
  {
    static int initial_sleep_until_setup = 1;
    if(initial_sleep_until_setup)
    {
      std::this_thread::sleep_for(std::chrono::milliseconds(2000)); 
      initial_sleep_until_setup = 0;
    }
    bool empty = overwriteQueue->IsEmpty();
    
    if(empty) --watch_dog_;
    else
    {

      std::string command = overwriteQueue->Pop();

      std::cout << "Robot command: \"" << command << "\"";
      std::cout << ", sending to: " << sendSocketAddress_ << '\n';
      socPtr->Send(sendSocketAddress_,
                   command.c_str(), sizeof(command.c_str()) );

      watch_dog_ = watch_dog_ > wait_limit_ ? wait_limit_ : watch_dog_ + 1;
    }

    if(watch_dog_>0)
    {
      std::this_thread::sleep_for(std::chrono::milliseconds(comm_time_));
    }
    else stop_=true;
  }

  void Start(const OverwriteQueue_Ptr& overwriteQueue, net::Socket *socPtr)
  {
    while(!stop_)
    {
      Send(overwriteQueue,socPtr);
    }
    if(stop_) std::cout << "warning communication has stopped!!!!\n" << std::endl;
  }
};


#endif
