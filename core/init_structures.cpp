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

#include "init_structures.h"

#define FEATURECOUNT 6

void init_stats(CvSize Img_Size, Statistics * S, bool allocate)
{
  static int range = 256;
  S->no_features = 6;
  S->nos = 0; //number of segments
  S->img_w = Img_Size.width;
  S->img_h = Img_Size.height;

  if(allocate)
  {
    S->id = new int[range]();
    S->size = new int[range]();
    S->gray_id = new int[range]();
    S->mean = new CvScalar[range]();
    S->stdDev = new CvScalar[range]();
    S->box = new CvRect[range]();

    S->priorTrue = new double[range]();
    S->priorFalse = new double[range]();
    S->P_GtgF = new double[range]();
    S->P_GfgF = new double[range]();

    S->P_FgGt = new double[range*S->no_features];
    S->P_FgGf = new double[range*S->no_features];

    S->prior_img = cvCreateImage(Img_Size,32,1);
  }
  else
  {
    for(int i=0; i<range; ++i)
    {
      S->id[i] = 0;
      S->size[i] = 0;
      S->gray_id[i] = 0;
      S->mean[i] = cvScalar(0,0,0,0);
      S->stdDev[i] = cvScalar(0,0,0,0);
      S->box[i] = cvRect(0,0,0,0);

      S->priorTrue[i] = 0.;
      S->priorFalse[i] = 0.;
      S->P_GtgF[i] = 0.;
      S->P_GfgF[i] = 0.;

      for (int j=0; j<S->no_features; j++)
      {
        S->P_FgGt[S->no_features*i + j] = 0.;
        S->P_FgGf[S->no_features*i + j] = 0.;
      }
    }

    cvZero(S->prior_img);
  }


  unsigned int iseed = (unsigned int)time(NULL);
  srand (iseed);

  double min = 0.01;
  //double max = 40.01;

  if(allocate)
  {
    S->L1 = new double[S->no_features];
    S->L0 = new double[S->no_features];
    S->gmax = new double[S->no_features];
    S->Z1 = new double[S->no_features];
    S->Z0 = new double[S->no_features];
    S->G_score = new double[S->no_features];
  }

  for (int i=0; i<S->no_features; i++)
  {
    S->L1[i] = min;
    S->L0[i] = min;
    S->gmax[i] = 2.;
    S->Z1[i] = (1/S->L1[i])*(1-exp(-S->L1[i]*S->gmax[i]));
    S->Z0[i] = (1/S->L0[i])*(exp(S->L0[i]*S->gmax[i])-1);
    S->G_score[i] = 0.;

    //S->L1[i] = randdouble(min, max);
    //S->gmax[i] = randdouble(1, 10);
    //printf("L1[%d] = %0.5f\n", i, S->L1[i]);
  }


  //Histogram Initialization
  static int dim_40 	= 40;
  static int range_40 	= 40;
  float range_40_arr[] = {float(0.),float(range_40-1)};
  float* range_40_ptr = range_40_arr;

  if(allocate)
  {
    for(int i=0;i<S->no_features;i++)
    {
      S->H_G1[i] = cvCreateHist(1, &dim_40, CV_HIST_ARRAY, &range_40_ptr, 1);
      S->H_G0[i] = cvCreateHist(1, &dim_40, CV_HIST_ARRAY, &range_40_ptr, 1);
      S->H_G1_DISP[i] = cvCreateHist(1, &dim_40, CV_HIST_ARRAY, &range_40_ptr, 1);
      S->H_G0_DISP[i] = cvCreateHist(1, &dim_40, CV_HIST_ARRAY, &range_40_ptr, 1);
    }
  }
  else
  {
    for(int i=0;i< S->no_features; i++)
    {
      cvClearHist(S->H_G1[i]);
      cvClearHist(S->H_G0[i]);
      cvClearHist(S->H_G1_DISP[i]);
      cvClearHist(S->H_G0_DISP[i]);
    }
  }

  //Histogram Initialization
  static int dim_9	= 9;
  static int dim_32 	= 32;

  static int range_256 	= 256;
  static int range_181 	= 181;

  float range_256_arr[] = {float(0.),float(range_256-1)};
  float range_181_arr[] = {float(0.),float(range_181-1)};
  float range_2pi_arr[] = {-CV_PI,CV_PI};

  float* range_256_ptr = range_256_arr;
  float* range_181_ptr = range_181_arr;
  float* range_2pi_ptr = range_2pi_arr;

  if(allocate){
    S->H_SF[0] = cvCreateHist( 1, &dim_32, CV_HIST_ARRAY, &range_256_ptr, 1 );
    S->H_SF[1] = cvCreateHist( 1, &dim_9, CV_HIST_ARRAY, &range_2pi_ptr, 1 );
    S->H_SF[2] = cvCreateHist( 1, &dim_32, CV_HIST_ARRAY, &range_181_ptr, 1 );
    S->H_SF[3] = cvCreateHist( 1, &dim_32, CV_HIST_ARRAY, &range_256_ptr, 1 );
    S->H_SF[4] = cvCreateHist(1, &dim_32, CV_HIST_ARRAY, &range_256_ptr, 1);
    S->H_SF[5] = cvCreateHist(1, &dim_32, CV_HIST_ARRAY, &range_256_ptr, 1);
  }
  else
  {
    for(int i=0;i<S->no_features;i++)
    {
      cvClearHist(S->H_SF[i]);
    }
  }
}

void release_stats(Statistics *S)
{
  delete [] S->id;
  delete [] S->size;
  delete [] S->gray_id;
  delete [] S->mean;
  delete [] S->stdDev;
  delete [] S->box;
  delete [] S->priorTrue;
  delete [] S->priorFalse;
  delete [] S->P_GtgF;
  delete [] S->P_GfgF;
  delete [] S->P_FgGt;
  delete [] S->P_FgGf;
  cvReleaseImage(&S->prior_img);

  delete [] S->L1;
  delete [] S->L0;
  delete [] S->gmax;
  delete [] S->Z1;
  delete [] S->Z0;
  delete [] S->G_score;

  for(int i = 0; i < S->no_features; ++i)
  {
    cvReleaseHist(&S->H_G1[i]);
    cvReleaseHist(&S->H_G0[i]);
    cvReleaseHist(&S->H_G1_DISP[i]);
    cvReleaseHist(&S->H_G0_DISP[i]);
    cvReleaseHist(&S->H_SF[i]);
  }
}


void init_boundary(CvSize S, Boundary * B)
{
    B->Bimg = cvCreateImage(S,8,1);

    B->boundaryPixels = new CvPoint[S.width];
    B->obstaclePolarCoord = new CvPoint2D32f[S.width];
    B->obstacleCartesianCoord = new CvPoint2D32f[S.width];
    B->obstacleDisplayCoord = new CvPoint[S.width];
    B->frontiers = new CvPoint3D32f[S.width];

    const int maxSegments = (S.width+1)/2;
    B->Median_Angle = new CvPoint3D32f[maxSegments];

    B->RobotDisplayCoord = cvPoint(200,390);
}

void release_boundary(Boundary *B)
{
  cvReleaseImage(&B->Bimg);
  delete [] B->boundaryPixels;
  delete [] B->obstaclePolarCoord;
  delete [] B->obstacleCartesianCoord;
  delete [] B->obstacleDisplayCoord;
  delete [] B->frontiers;
  delete [] B->Median_Angle;
}

void init_model(const CvSize& S, const CvRect& SafeRegion, Model *M)
{
    M->safeRegionMask = cvCreateImage(S,8,1);
    cvZero(M->safeRegionMask);

    cvSetImageROI( M->safeRegionMask, SafeRegion );
    cvSet(M->safeRegionMask,cvScalarAll(255),0);
    cvResetImageROI( M->safeRegionMask );

    const int n = 1; //number of model regions
    M->mean = new CvScalar[n]();
    M->stdDev = new CvScalar[n]();

    //Histogram Initialization
    static int dim_9	= 9;
    static int dim_32 	= 32;

    static int range_256 	= 256;
    static int range_181 	= 181;

    float range_256_arr[] = {float(0),float(range_256-1)};
    float range_181_arr[] = {float(0),float(range_181-1)};
    float range_2pi_arr[] = {-CV_PI,CV_PI};

    float* range_256_ptr = range_256_arr;
    float* range_181_ptr = range_181_arr;
    float* range_2pi_ptr = range_2pi_arr;

    M->dim = new int[FEATURECOUNT];
    for(int i=0; i<FEATURECOUNT; i++)
    {
      if(i==1) M->dim[i] = 9;
      else     M->dim[i] = 32;
    }

    M->histograms[0] = cvCreateHist(1, &dim_32, CV_HIST_ARRAY, &range_256_ptr);
    M->histograms[1] = cvCreateHist(1, &dim_9 , CV_HIST_ARRAY, &range_2pi_ptr);
    M->histograms[2] = cvCreateHist(1, &dim_32, CV_HIST_ARRAY, &range_181_ptr);
    M->histograms[3] = cvCreateHist(1, &dim_32, CV_HIST_ARRAY, &range_256_ptr);
    M->histograms[4] = cvCreateHist(1, &dim_32, CV_HIST_ARRAY, &range_256_ptr);
    M->histograms[5] = cvCreateHist(1, &dim_32, CV_HIST_ARRAY, &range_256_ptr);

    M->histograms_DISP[0] = cvCreateHist(1, &dim_32, CV_HIST_ARRAY, &range_256_ptr);
    M->histograms_DISP[1] = cvCreateHist(1, &dim_9 , CV_HIST_ARRAY, &range_2pi_ptr);
    M->histograms_DISP[2] = cvCreateHist(1, &dim_32, CV_HIST_ARRAY, &range_181_ptr);
    M->histograms_DISP[3] = cvCreateHist(1, &dim_32, CV_HIST_ARRAY, &range_256_ptr);
    M->histograms_DISP[4] = cvCreateHist(1, &dim_32, CV_HIST_ARRAY, &range_256_ptr);
    M->histograms_DISP[5] = cvCreateHist(1, &dim_32, CV_HIST_ARRAY, &range_256_ptr);
}

void release_model(Model *M)
{
  cvReleaseImage(&M->safeRegionMask);
  delete [] M->mean;
  delete [] M->stdDev;
  delete [] M->dim;

  for(int i=0; i<FEATURECOUNT; i++)
  {
    cvReleaseHist(&M->histograms[i]);
    cvReleaseHist(&M->histograms_DISP[i]);
  }
}

void init_features(CvSize S, Features * F)
{

    F->mag = cvCreateImage(S,8,1);
    F->ang32 = cvCreateImage(S,32,1);
    F->P_ang = cvCreateImage(S,32,1);

    F->hsv   = cvCreateImage( S, 8, 3);
    F->lab   = cvCreateImage( S, 8, 3);
    F->YCrCb   = cvCreateImage( S, 8, 3);
    F->hue   = cvCreateImage( S, 8, 1);
    F->sat   = cvCreateImage( S, 8, 1);
    F->val   = cvCreateImage( S, 8, 1);
    F->Cr   = cvCreateImage( S, 8, 1);
    F->a   = cvCreateImage( S, 8, 1);
    F->Cb   = cvCreateImage( S, 8, 1);
    F->iic   = cvCreateImage( S, 8, 1);
    F->P_hue   = cvCreateImage( S, 32, 1);
    F->P_sat   = cvCreateImage( S, 32, 1);
    F->P_val   = cvCreateImage( S, 32, 1);

    F->lbp  = cvCreateImage( S, 8, 1);
    F->P_lbp = cvCreateImage( S, 32, 1);

    F->post0 = cvCreateImage( S, 32, 1);
    F->post1 = cvCreateImage( S, 32, 1);
    F->post_ratio = cvCreateImage( S, 32, 1);

    F->bin_class_result = cvCreateImage( S, 8, 1);

    for(int i=0;i<FEATURECOUNT;i++){
        F->P_X1[i] = cvCreateImage( S, 32, 1);
        F->P_X0[i] = cvCreateImage( S, 32, 1);
    }

}

void release_features(Features *F)
{
  cvReleaseImage(&F->mag );
  cvReleaseImage(&F->ang32 );
  cvReleaseImage(&F->P_ang );
  cvReleaseImage(&F->hsv );
  cvReleaseImage(&F->lab );
  cvReleaseImage(&F->YCrCb );
  cvReleaseImage(&F->hue );
  cvReleaseImage(&F->sat );
  cvReleaseImage(&F->val );
  cvReleaseImage(&F->Cr );
  cvReleaseImage(&F->a );
  cvReleaseImage(&F->Cb );
  cvReleaseImage(&F->iic );
  cvReleaseImage(&F->P_hue );
  cvReleaseImage(&F->P_sat );
  cvReleaseImage(&F->P_val );
  cvReleaseImage(&F->lbp );
  cvReleaseImage(&F->P_lbp );
  cvReleaseImage(&F->post0 );
  cvReleaseImage(&F->post1 );
  cvReleaseImage(&F->post_ratio );
  cvReleaseImage(&F->bin_class_result );

  for(int i=0;i<FEATURECOUNT;i++){
    cvReleaseImage(&F->P_X1[i]);
    cvReleaseImage(&F->P_X0[i]);
  }
}

