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

#include <iostream>
#include <string>
#include <cassert>

#include "OverwriteQueue.h"

template <class Dtype>
OverwriteQueue<Dtype>::OverwriteQueue(int size):size_(size){}
template class OverwriteQueue<int>;
template class OverwriteQueue<std::string>;


template< class Dtype >
bool OverwriteQueue<Dtype>::IsEmpty()
{
  std::lock_guard<std::mutex> lck(mtx);
  return qbuf_.empty();
}
template bool OverwriteQueue<int>::IsEmpty();
template bool OverwriteQueue<std::string>::IsEmpty();


template <class Dtype>
Dtype OverwriteQueue<Dtype>::Pop()
{
  std::lock_guard<std::mutex> lck(mtx);
  assert(!qbuf_.empty());

  Dtype contents = qbuf_.front();
  qbuf_.pop();
  return contents;
}
template int OverwriteQueue<int>::Pop();
template std::string OverwriteQueue<std::string>::Pop();


template< class Dtype >
void OverwriteQueue<Dtype>::Reset(const Dtype& contents)
{
  std::lock_guard<std::mutex> lck (mtx);
  Empty();
  Fill(contents);
}
template void OverwriteQueue<int>::Reset(const int& contents);
template void OverwriteQueue<std::string>::Reset(const std::string& contents);


template< class Dtype >
void OverwriteQueue<Dtype>::Fill(const Dtype& contents)
{
  assert(qbuf_.empty()==true);
  for(int i=0;i<size_;++i)
  {
    qbuf_.push(contents);
  }
}
template void OverwriteQueue<int>::Fill(const int &contents);
template void OverwriteQueue<std::string>::Fill(const std::string& contents);

template< class Dtype >
void OverwriteQueue<Dtype>::Empty()
{
  while(!qbuf_.empty())
  {
    qbuf_.pop();
  }
}
template void OverwriteQueue<int>::Empty();
template void OverwriteQueue<std::string>::Empty();
