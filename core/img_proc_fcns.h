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

#ifndef IMG_PROC_FCNS
#define IMG_PROC_FCNS

#include "init_structures.h"

void init_images_img_proc(const CvSize& S);
void release_images_img_proc(void);

void init_windows_img_proc(CvSize S);

void GetModel(IplImage *gray_img, Features *F, Model *StatModelPtr, bool dynamic);

void GraphBasedSegmentation(IplImage *seg, IplImage *gray_source);
void SuperPixelStats(IplImage *gbs, IplImage *gray, Statistics *sts);

//COLOUR
void ColorHistAnalysis(IplImage *src, IplImage *hsv, IplImage *col_cue, IplImage *col_prob, CvRect Safe);
//LBP
void convertGray2LBP(IplImage *gray, IplImage *LBP);
void LBPAnalysis(IplImage* LBP, IplImage* LBP_CUE, IplImage* LBP_PROB, CvRect Safe);

void Contours(const IplImage *binary);
void DrawContours(IplImage* contour, CvScalar color, CvRect S);
void ReleaseContours(void);

void getEdgeMagOri(IplImage* gray, IplImage* mag32, IplImage* ang32);
void ProbAnalysis2(Features *F, Statistics* sts, IplImage* gbs);
void UpdatePrior(IplImage *gbs, Statistics *sts, Features *F);
//void UpdateModel(Model* Current, Model* Prev);
//void UpdatePrevModel(Model* Current, Model* Prev);
void FeatureAnalysis(Features *F, Model* M, Statistics *sts, IplImage *gbs, bool dynamic);

bool CheckConvergence(Statistics* S, int em);

void UpdateParams(IplImage* bin, Statistics *S, Features *F, bool dynamic);

CvRect BoundingRect(const IplImage *im);

void combine_channels(IplImage* Cr, IplImage* Cb, IplImage* a, IplImage* iic);

void FindObstacleBoundary(IplImage* Out);
void ExtractBoundary(CvSize S, Boundary *B);
void CalculateDistances(CvSize S, Boundary *B, CamCalib camera, BotCalib bot);
void InterpretDepthArray(CvSize S, Boundary *B, BotCalib bot);



#endif
