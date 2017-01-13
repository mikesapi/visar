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

//ms_overwrite_safe_buffer.h
#ifndef __MS_OVERWRITE_SAFE_BUFFER__
#define __MS_OVERWRITE_SAFE_BUFFER__

#include <memory>
#include <mutex>
#include <queue>

template <class Dtype>
class OverwriteQueue
{
private:
  std::mutex mtx;
  std::queue<Dtype> qbuf_;
  int size_;
  
public:
  OverwriteQueue(int size);

public:
  Dtype Pop();
  void Reset(const Dtype& contents);
  bool IsEmpty();

private:
  void Fill(const Dtype& contents);
  void Empty();
};

typedef std::shared_ptr<OverwriteQueue<std::string> > OverwriteQueue_Ptr;
typedef std::shared_ptr<const OverwriteQueue<std::string> > OverwriteQueue_CPtr;

#endif
