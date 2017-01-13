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

#ifndef INIT_STRUCTURES
#define INIT_STRUCTURES

#include <iostream>

#include "structures.h"

void init_stats(CvSize Img_Size, Statistics *S, bool allocate);
void release_stats(Statistics *S);

void init_boundary(CvSize S, Boundary *B);
void release_boundary(Boundary *B);

void init_model(const CvSize& S, const CvRect& SafeRegion, Model *M);
void release_model(Model *M);

void init_features(CvSize S, Features *F);
void release_features(Features *F);

#endif
