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

#include "img_proc_fcns.h"

#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv/cxcore.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <stdio.h>

#include "../graphbasedsegmentation/image.h"
#include "../graphbasedsegmentation/misc.h"
#include "../graphbasedsegmentation/segment-image.h"


//Binary image zero & one

#define FEATURECOUNT 6

#define Window_W 1.02*p.proc_W //appriximate wht window width and hight as a function of the frame size
#define Window_H 1.3*(p.proc_H)+20

#define DISPLAY_IMAGE_NAME(img,NAME)		cvNamedWindow(NAME); cvShowImage(NAME, img);

#define DISPLAY_IMAGE_XY(R,img,X,Y)		if(R){cvNamedWindow(#img); cvMoveWindow(#img, int(round(X*Window_W)), int(round(Y*Window_H))) ;} cvShowImage(#img, img);

#define DISPLAY_IMAGE_XY_NAME(R,img,X,Y,NAME)		if(R){cvNamedWindow(NAME); cvMoveWindow(NAME, int(round(X*Window_W)), int(round(Y*Window_H))) ;} cvShowImage(NAME, img);

extern Params p;

const char * LBP_HIST 		= "LBP Image Histogram"; // window name
const char * iic_HIST 		= "IIC Image Histogram"; // window name
const char * HIST_HUE = "Hue Hist Model";
const char * HIST_SAT = "Sat Hist Model";
const char * HIST_GRAD_MAG = "Grag Mag Hist Model";
const char * HIST_ANGLE = "Edge Grad Angle";

const char * P1name[FEATURECOUNT] = {"MAG1 G","ANG1 G","HUE1 G","SAT1 G","LBP1 G","IIC1 G"};
const char * P0name[FEATURECOUNT] = {"MAG0 G","ANG0 G","HUE0 G","SAT0 G","LBP0 G","IIC0 G"};


// min/max values for Variance/Saturation mask
static int vmin = 10, vmax = 256, smin = 5;

CvHistogram *hist = NULL;

//Images for use as masks eg:superpixels
IplImage *mask[FEATURECOUNT];

//images to use as masks with colour analysis
IplImage *bin[FEATURECOUNT];

//Images to display histograms
IplImage *HistImgH, *HistImgS, *HistImgMag;

IplImage *sobel[FEATURECOUNT], *AngleHist_img; //Images to be used with edge analysis

IplImage *inv_prob[FEATURECOUNT], *temp[FEATURECOUNT], *PG[FEATURECOUNT], *PG_prev[FEATURECOUNT], *PG1_DISP[FEATURECOUNT], *PG0_DISP[FEATURECOUNT]; //Images for use when calculating probabilities

IplImage* DEPTH_MAP;//Display depth array

IplImage* LBPhist_img, *iichist_img;

IplImage* GhistImg;
CvHistogram *Ghist = NULL;

IplImage* GhistImg2;
CvHistogram *Ghist2 = NULL;
CvHistogram *Ghist2DISP = NULL;

//Display superpixel statistics
IplImage *stats_disp, *hist_temp;


inline void POLAR2CART(const CvPoint2D32f *P, CvPoint2D32f *C)
{
  C->x = (P->x)*(sin(P->y));
  C->y = (P->x)*(cos(P->y));
}

inline void CART2DISPLAY(CvPoint2D32f *C, CvPoint *D, CvPoint *Offset, double scale)
{
  D->x = -(C->x)*(scale) + Offset->x;
  D->y = -(C->y)*(scale) + Offset->y;
}

inline CvScalar hue2rgb(float hue)
{
    // hue to RGB conversion : coverts a given hue value to a RGB triplet for
    // display
    // parameters:
    // hue - hue value in range 0 to 180 (OpenCV implementation of HSV)
    // return value - CvScalar as RGB triple
    // taken from OpenCV 1.0 camshiftdemo.c example
    int rgb[3], p, sector;
    static const int sector_data[][3]= {{0,2,1}, {1,2,0}, {1,0,2}, {2,0,1}, {2,1,0}, {0,1,2}};
    hue *= 0.033333333333333333333333333333333f;
    sector = cvFloor(hue);
    p = cvRound(255*(hue - sector));
    p ^= sector & 1 ? 255 : 0;

    rgb[sector_data[sector][0]] = 255;
    rgb[sector_data[sector][1]] = 0;
    rgb[sector_data[sector][2]] = p;

    return cvScalar(rgb[2], rgb[1], rgb[0],0);
}

inline void PrintGHistogram(int hist_size, CvHistogram *Hist, IplImage* Hist_img, const char * Window, bool flag, int X, int Y)
{
  int bin_w = cvRound((double)Hist_img->width/hist_size);
  cvNormalizeHist(Hist, (10*Hist_img->height));

  if(flag==0)
  {
    cvSet(Hist_img, cvScalarAll(255));
    for(int i = 0; i < hist_size; i++ )
    {
      int val = cvRound( cvGetReal1D(Hist->bins,i));
      cvRectangle( Hist_img, cvPoint(i*bin_w,Hist_img->height),
                      cvPoint((i+1)*bin_w,Hist_img->height - val),
                      cvScalarAll(0), -1, 8, 0 );
    }
  }

  DISPLAY_IMAGE_XY_NAME(p.refresh,Hist_img,X,Y,Window)
}

inline void PrintHistogram(int hist_size, CvHistogram *Hist, IplImage* Hist_img, const char * Window, bool flag, int X, int Y)
{
    int bin_w = cvRound((double)Hist_img->width/hist_size);
    cvNormalizeHist(Hist, (Hist_img->height));

    //grab the min and max values and their indeces
    //cvGetMinMaxHistValue( Hist, NULL, &max_value, NULL, NULL);
    //scale the bin values so that they will fit in the image representation
    //cvScale( Hist->bins, Hist->bins, ((double)Hist_img->height - 10.)/(max_value), 0 );

    if(flag==0){
        cvSet( Hist_img, cvScalarAll(255), 0 );
        for(int i = 0; i < hist_size; i++ ){
            int val = cvRound( cvGetReal1D(Hist->bins,i));
            cvRectangle( Hist_img, cvPoint(i*bin_w,Hist_img->height),
                            cvPoint((i+1)*bin_w,Hist_img->height - val),
                            cvScalarAll(0), -1, 8, 0 );
        }

    }

    else if(flag==1){
        cvZero( Hist_img );
        for(int i = 0; i < hist_size; i++ ){
            int val = cvRound( cvGetReal1D(Hist->bins,i));
            CvScalar color = hue2rgb(i*180.f/hist_size);
            cvRectangle( Hist_img, cvPoint(i*bin_w,Hist_img->height),
                            cvPoint((i+1)*bin_w,Hist_img->height - val),
                            color, -1, 8, 0 );
        }

    }

    if (p.debug){
        cvShowImage(Window, Hist_img);
        DISPLAY_IMAGE_XY_NAME(p.refresh,Hist_img,X,Y,Window)
    }
}

void init_images_img_proc(const CvSize& S)
{
  //View_Histogram Graph Based Segmentation
  int GBShist_size = 256;			// size of histogram (number of bins)
  float range_0[] = {0,256};
  float* ranges[] = {range_0};
  hist = cvCreateHist(1, &GBShist_size, CV_HIST_ARRAY, ranges, 1);

  //Histogram Analysis
  CvSize HistSize = cvSize(128,100);
  HistImgH = cvCreateImage( HistSize, 8, 3 );
  HistImgS = cvCreateImage( HistSize, 8, 3 );
  HistImgMag = cvCreateImage( HistSize, 8, 3 );
  AngleHist_img = cvCreateImage(HistSize, 8, 1);
  LBPhist_img = cvCreateImage(HistSize, 8, 1);
  iichist_img = cvCreateImage(HistSize, 8, 1);

  GhistImg = cvCreateImage(cvSize(512,100), 8, 1);
  GhistImg2 = cvCreateImage(cvSize(120,100), 8, 1);

  int i;
  for(i=0;i<FEATURECOUNT;++i)
  {
    bin[i] = cvCreateImage( S, 8, 1);
    mask[i] = cvCreateImage( S, 8, 1);
    sobel[i] = cvCreateImage( S, IPL_DEPTH_8U, 1);
    inv_prob[i] = cvCreateImage( S, 32, 1);
    temp[i] = cvCreateImage( S, 32, 1);
    PG[i] = cvCreateImage( S, 32, 1);
    PG_prev[i] = cvCreateImage( S, 32, 1);
    cvZero(PG_prev[i]);
    PG1_DISP[i] = cvCreateImage( cvSize(120,100), 8, 3);
    PG0_DISP[i] = cvCreateImage( cvSize(120,100), 8, 3);
  }

  static int vrange = 256;

  static int dim_40 	= 40;
  static int range_40 	= 40;
  float range_40_arr[] = {float(0),float(range_40-1)};
  float* range_40_ptr = range_40_arr;

  float range_log_arr[] = {-30,30};
  float* range_log_ptr = range_log_arr;
  Ghist = cvCreateHist( 1, &vrange, CV_HIST_ARRAY, &range_log_ptr, 1 );
  Ghist2 = cvCreateHist( 1, &dim_40, CV_HIST_ARRAY, &range_40_ptr, 1 );
  Ghist2DISP = cvCreateHist( 1, &dim_40, CV_HIST_ARRAY, &range_40_ptr, 1 );

  DEPTH_MAP = cvCreateImage(cvSize(400,400), 8, 3);

  stats_disp = cvCreateImage( S, 8, 3);
  hist_temp = cvCreateImage( cvSize(120,100), 8, 3);
}

void release_images_img_proc( void )
{
    cvReleaseImage( &stats_disp );

    cvReleaseImage( &AngleHist_img );
    cvReleaseImage( &LBPhist_img );

    cvReleaseImage( &DEPTH_MAP );

    cvReleaseHist( &hist );

    //Histogram Analysis
    cvReleaseImage( &HistImgH );
    cvReleaseImage( &HistImgS );
    cvReleaseImage( &HistImgMag );

    for(int i=0;i<4;++i){
        cvReleaseImage( &bin[i]  );
        cvReleaseImage( &mask[i] );
        cvReleaseImage( &sobel[i]);
        cvReleaseImage( &inv_prob[i]);
        cvReleaseImage( &temp[i]);
    }

    cvDestroyAllWindows();
}

CvMemStorage* storage = cvCreateMemStorage();
CvSeq* first_contour = NULL;
CvSeq* poly_result = NULL;

void Contours(const IplImage *binarymask)
{
  IplImage *bcopy = cvCreateImage(cvSize(binarymask->width,binarymask->height), 8, 1);
  cvCopy(binarymask, bcopy);

  cvFindContours(
    bcopy,
    storage,
    &first_contour,
    sizeof(CvContour),
    CV_RETR_CCOMP,
    CV_CHAIN_APPROX_NONE,
    cvPoint(0,0)
  );

  cvReleaseImage(&bcopy);
}

void DrawContours(IplImage* contour, CvScalar color, CvRect S)
{
  cvZero(contour);
  for(; first_contour != 0; first_contour = first_contour->h_next)
  {
    cvDrawContours(
    contour,
    first_contour,
    color,
    color,
    2,
    CV_FILLED,
    8
    );
  }

  cvRectangle(contour, cvPoint(S.x,S.y),
  cvPoint(S.x+S.width,S.y+S.height), cvScalarAll(255), 1, 8, 0);
}

void ReleaseContours(void)
{
  if(first_contour != NULL)
  {
    cvClearSeq(first_contour);
  }

  cvClearMemStorage(storage);
}

void GraphBasedSegmentation(IplImage* seg, IplImage* gray_source)
{
  static int int_sigma = 5;
  static int int_k = 100;
  static int min_size = 100;

  static int num_ccs;
  static float sigma;
  static float k;

  sigma = (float)int_sigma*0.1;
  k = (float)int_k;
  SegmentImage(seg, gray_source, sigma, k, min_size, &num_ccs);

  if(p.debug)
  {
    const char *imageName = "seg";
    DISPLAY_IMAGE_XY_NAME(p.refresh, seg, 1,1, imageName);

    if(p.refresh)
    {
      cvCreateTrackbar( "sigma", imageName, &int_sigma , 100, NULL );
      cvCreateTrackbar( "k", imageName, &int_k, 1000, NULL );
      cvCreateTrackbar( "min_size", imageName, &min_size, 400, NULL );
    }
  }
}

void GradMagAng32(IplImage* X, IplImage* Y, IplImage* G, IplImage* A)
{
    static int step       = A->widthStep/sizeof(float);
    static int channels   = A->nChannels;

    float * XPixelData = (float *)(X->imageData);
    float * YPixelData = (float *)(Y->imageData);
    float * GPixelData = (float *)(G->imageData);
    float * APixelData = (float *)(A->imageData);

    for (int y = 0; y < X->height; y++) {
        for (int x = 0; x < X->width; x++) {

            //int X_index = y*step+x*channels+0;
            float X_pixel = XPixelData[y*step+x*channels+0];

            //int Y_index = y*step+x*channels+0;
            float Y_pixel = YPixelData[y*step+x*channels+0];

            //int MAG = cvRound( sqrt( X_pixel*X_pixel + Y_pixel*Y_pixel ) );
            float MAG = abs(X_pixel) + abs(Y_pixel);
            GPixelData[y*step+x*channels+0] = 255 - MAG;
            GPixelData[y*step+x*channels+0] = MAG;
            //calculate edge orientations

            float ANG = atan2( Y_pixel , X_pixel);

            //if(ANG < 0) printf("This is smaller than zero: %0.2f\n", ANG);
            APixelData[y*step+x*channels+0] = ANG;

        }
    }

}

CvRect BoundingRect(const IplImage* im)
{
  unsigned char * imPixelData = (unsigned char *)(im->imageData);
  int width = im->width;
  int height = im->height;

  int xmin=width+1, xmax=-1, ymin=height+1, ymax=-1;
  for(int y = 0; y < height; ++y)
  {
    for(int x = 0; x < width; ++x)
    {
      int index = (y*width+x)*im->nChannels;
      if(imPixelData[index])
      {
        xmin = xmin > x ? x : xmin;
        xmax = xmax < x ? x : xmax;
        ymin = ymin > y ? y : ymin;
        ymax = ymax < y ? y : ymax;
      }
    }
  }

  int rwidth = xmax-xmin;
  int rheight = ymax-ymin;
  return cvRect(xmin, ymin, rwidth, rheight);
}

void combine_channels(IplImage* Cr, IplImage* Cb, IplImage* a, IplImage* iic){
    unsigned char * CrPixelData = (unsigned char *)(Cr->imageData);
    unsigned char * CbPixelData = (unsigned char *)(Cb->imageData);
    unsigned char * aPixelData = (unsigned char *)(a->imageData);
    unsigned char * iicPixelData = (unsigned char *)(iic->imageData);

    cvZero(iic);

    for (int y = 0; y < iic->height; y++) {
        for (int x = 0; x < iic->width; x++) {

            //int Cr_index = (y*Cr.width+x);
            //int Cb_index = (y*Cb.width+x);
            //int a_index = (y*a.width+x);
            int index = (y*iic->width+x);
            iicPixelData[index] = cvRound((CrPixelData[index] + CbPixelData[index] + 2*(aPixelData[index]))/4);
        }
    }

}

void getEdgeMagOri(IplImage* gray, IplImage* mag, IplImage* ang32){

    cvSobel(gray, temp[0], 1, 0, CV_SCHARR); //Find edges in x direction
    cvConvertScaleAbs(temp[0], sobel[1], 0.2, 0); //convert 16bit image to 8-bit

    cvSobel(gray, temp[1], 0, 1, CV_SCHARR);
    cvConvertScaleAbs(temp[1], sobel[2], 0.2, 0); //convert 16bit image to 8-bit

    GradMagAng32( temp[0], temp[1], temp[2], ang32 );


    cvNormalize( temp[2],temp[2], 0, 255, CV_MINMAX);
    cvNormalize(ang32,temp[3], 0, 255, CV_MINMAX);
    cvAbsDiffS(temp[2], temp[2], cvScalar(255));
    cvConvertScaleAbs(temp[2], mag, 1, 0);
    cvConvertScaleAbs(temp[3], sobel[0], 1, 0);

    if(p.debug)
    {
      DISPLAY_IMAGE_XY(p.refresh, sobel[0], 4,0);
      DISPLAY_IMAGE_XY(p.refresh, mag, 3,0);
    }

}

double GetBit(CvSize S, unsigned char * PixelData, int p_index, int pc_index, int py, int px)
{
  if ( px<0 || py<0 || px>=S.width || py>=S.height) //check that image point indexes are within range
          return 0;
  else
  {
      int p_val = PixelData[p_index]; // pixel value
      int p_c = PixelData[pc_index];  // centre pixel value
      if (p_val>=p_c)
          return 1;
      else
          return 0;
  }
}

void convertGray2LBP(IplImage* gray, IplImage* LBP){

    CvSize S = cvSize(gray->width, gray->height);
    unsigned char * GrayPixelData = (unsigned char *)(gray->imageData);
    //const int * imgData = (const int *)(gray->imageData);
    unsigned char * LbpPixelData = (unsigned char *)(LBP->imageData);

    cvZero(LBP);

    for (int y = 0; y < gray->height; y++) {
        for (int x = 0; x < gray->width; x++) {

            int p0x,p1x,p2x,p3x,p4x,p5x,p6x,p7x;
            int p0y,p1y,p2y,p3y,p4y,p5y,p6y,p7y;

            p0x=x-1;
            p0y=y-1;
            p1x=x;
            p1y=y-1;
            p2x=x+1;
            p2y=y-1;
            p3x=x+1;
            p3y=y;
            p4x=x+1;
            p4y=y+1;
            p5x=x;
            p5y=y+1;
            p6x=x-1;
            p6y=y+1;
            p7x=x-1;
            p7y=y;

            int pc_index = (y*gray->width+x); //centre pixel index

            int p0_index = (p0y*gray->width+p0x); //p0 index..
            int p1_index = (p1y*gray->width+p1x);
            int p2_index = (p2y*gray->width+p2x);
            int p3_index = (p3y*gray->width+p3x);

            int p4_index = (p4y*gray->width+p4x);
            int p5_index = (p5y*gray->width+p5x);
            int p6_index = (p6y*gray->width+p6x);
            int p7_index = (p7y*gray->width+p7x); //..p7 index

            double b0 = 128*GetBit(S, GrayPixelData, pc_index, p0_index, p0y, p0x);
            double b1 = 64*GetBit(S, GrayPixelData, pc_index, p1_index, p1y, p1x);
            double b2 = 32*GetBit(S, GrayPixelData, pc_index, p2_index, p2y, p2x);
            double b3 = 16*GetBit(S, GrayPixelData, pc_index, p3_index, p3y, p3x);

            double b4 = 8*GetBit(S, GrayPixelData, pc_index, p4_index, p4y, p4x);
            double b5 = 4*GetBit(S, GrayPixelData, pc_index, p5_index, p5y, p5x);
            double b6 = 2*GetBit(S, GrayPixelData, pc_index, p6_index, p6y, p6x);
            double b7 = 1*GetBit(S, GrayPixelData, pc_index, p7_index, p7y, p7x);

            //int decimal = cvRound(b0+b1+b2+b3 + b4+b5+b6+b7);

            LbpPixelData[pc_index] = cvRound(b0+b1+b2+b3 + b4+b5+b6+b7); //pixel intensity is equal to the decimal number equivalent to the binary 8-bit pattern.

        }
    }

}


inline double GetPrior(int h, CvRect* R){
    //const static int hmax = h-1;
    const static double lambda = 3./(double)h;
    double height = (double)(h - (R->y + (R->height)/2)) ;

    return exp(-lambda*height);
}

inline double EXP_DIST_1(double z, double l, double g){
    return (1/z)*exp(-l*g);
}

inline double EXP_DIST_0(double z, double l, double g){
    return (1/z)*exp(l*g);
}

void draw_dist(int range, double gmax, double Z, double L, IplImage* Img, bool T){
    int bin_w = cvRound((double)Img->width/(range));

    double value;
    for(int g=0; g<cvRound(gmax); g++){


        if (T==1) {
        value = 50*(EXP_DIST_1(Z, L, g));
            }
        if (T==0){
        value = 50*(EXP_DIST_0(Z, L, g));
        }

        value = value > 0 ? value : 0.001;
        value = value > Img->height ? Img->height : value;
        //double value = 100*(value1+value0)/2;

                cvRectangle( Img, cvPoint(cvRound(g)*bin_w,Img->height),
                        cvPoint((cvRound(g)+1)*bin_w,cvRound(Img->height - value)),
                        CV_RGB(200,0,0), -1, 8, 0 );

    }

}

inline double Gstat(CvHistogram* H1, CvHistogram* H2, int size){

    double G = 0.;
    const static double SMALL = 0.000001;
    cvNormalizeHist(H1, 1);
    cvNormalizeHist(H2, 1);

    for(int i=0; i<size; i++){
        double M = cvGetReal1D(H1->bins,i); //sample
        double S = cvGetReal1D(H2->bins,i); //model

        if(M == 0) M = SMALL;
        if(S == 0) S = SMALL;
        //printf("M=%0.5f, S=0.5f\n", M, S);
        G += 2*(S*(log (S)) - S*(log (M)));
        //G += - S*(log (M));
        //printf("G=%0.5f\n",G);
        //cvWaitKey(0);
    }

    //printf("G=%0.4f\n", G);

    //return 1/exp(0.2*G);
    return G;
}

void SuperPixelStats(IplImage *gbs, IplImage *gray, Statistics *sts)
{
  int hist_size = 256;			// size of histogram (number of bins)
  int k = 0;

  cvCalcHist(&gbs, hist, 0, NULL);
  cvSet(stats_disp, cvScalarAll(255), 0);

  for(int i = 0; i < hist_size; i++)
  {
    //draw histogram, where bin size indicates superpixel size
    if(cvRound(cvGetReal1D(hist->bins,i)) > 0)
    {
      sts->id[k] = k;
      sts->gray_id[k] = i;
      cvCmpS(gbs, i, mask[0], CV_CMP_EQ); //mask out image segment
      sts->size[k] = cvCountNonZero(mask[0]); //count number of pixels in current segment
      cvAvgSdv(gray, &sts->mean[k], &sts->stdDev[k], mask[0]);

      sts->box[k] = BoundingRect(mask[0]);

      sts->priorTrue[k] = GetPrior(sts->img_h, &sts->box[k]);
      sts->priorFalse[k] =  1. - sts->priorTrue[k];

      if(p.debug)
      {
        //DISPLAY
        cvSet(stats_disp, cvScalarAll(sts->mean[k].val[0]), mask[0]);
        cvSet(sts->prior_img, cvScalar(sts->priorTrue[k]), mask[0]);

        cvRectangle(stats_disp,	cvPoint(sts->box[k].x,sts->box[k].y),
                    cvPoint(sts->box[k].x+sts->box[k].width,
                            sts->box[k].y+sts->box[k].height),
                    CV_RGB(255,0,0), 1, 8, 0);

        cvLine(stats_disp, cvPoint(sts->box[k].x + (sts->box[k].width)/2,
                                   sts->box[k].y +(sts->box[k].height)/2),
                           cvPoint(sts->img_w/2,sts->img_h),
                           CV_RGB(255,0,0), 1, 8, 0);
      }

      k++;
    }
  }
  sts->nos = k;

  if(p.debug)
  {
    DISPLAY_IMAGE_XY(p.refresh, stats_disp, 8, 4);
  }
}

void UpdatePrior(IplImage *gbs, Statistics *sts, Features* F)
{
  CvScalar mean1 = cvScalar(0,0,0,0);
  CvScalar mean0 = cvScalar(0,0,0,0);
  static double alpha = 0.6;
  static double beta = 1.- alpha;

  cvAvgSdv(F->post1, &mean1, 0, NULL);

  if(mean1.val[0] > 0)
  {
    for(int i = 0; i < sts->nos; i++)
    {
    cvCmpS(gbs, sts->gray_id[i], mask[0], CV_CMP_EQ); //mask out image segmen

    cvSmooth(F->post1, F->post1, CV_GAUSSIAN, 5, 5);
    cvSmooth(F->post0, F->post0, CV_GAUSSIAN, 5, 5);

    cvAvgSdv(F->post1, &mean1, 0, mask[0]);
    cvAvgSdv(F->post0, &mean0, 0, mask[0]);

    if(mean1.val[0] > 0)
    {
      double prior1 = alpha*sts->priorTrue[i] + beta*mean1.val[0];
      //sts->priorFalse[i] =  1. - sts->priorTrue[i];
      double prior0 = alpha*sts->priorFalse[i] + beta*mean0.val[0];

      sts->priorTrue[i] = prior1/(prior1+prior0);
      sts->priorFalse[i] = prior0/(prior1+prior0);
    }

    cvSet(sts->prior_img, cvScalar(sts->priorTrue[i]), mask[0]);
    }
  }
}

void GetModel(IplImage* gray, Features* F, Model* M, bool dynamic)
{
  cvInRangeS( F->hsv, cvScalar(0,   smin, min(vmin,vmax), 0),
                      cvScalar(180, 256,  max(vmin,vmax), 0), bin[0] );

  cvAnd(bin[0], M->safeRegionMask, bin[1]);

  int loop(-1);
  if(dynamic) loop = 200;

  static int j=0;

  int acc = 0;
  if(j < loop)
  {
    ++j;
    // This ensures that the histograms are cleared in the beginning when it is allocated.
    if(j > 1) acc=1;
  }
  else
  {
    acc=0;
    j=0;
  }


  cvCalcHist(&F->mag, 	M->histograms[0], acc, M->safeRegionMask );
  cvCalcHist(&F->ang32, M->histograms[1], acc, M->safeRegionMask );
  cvCalcHist(&F->hue, 	M->histograms[2], acc, bin[1] );
  cvCalcHist(&F->sat, 	M->histograms[3], acc, M->safeRegionMask );
  cvCalcHist(&F->lbp, 	M->histograms[4], acc, M->safeRegionMask );
  cvCalcHist(&F->iic, 	M->histograms[5], acc, M->safeRegionMask );

  for(int i = 0; i < FEATURECOUNT; ++i)
  {
    cvCopyHist(M->histograms[i], &M->histograms_DISP[i]);
  }

  if(p.debug)
  {
    int dim_9 = 9;
    int dim_32 = 32;
    int row = 1;
    PrintHistogram(dim_32, M->histograms_DISP[0], HistImgMag,    HIST_GRAD_MAG, 0, 5, row);
    PrintHistogram(dim_9,  M->histograms_DISP[1], AngleHist_img, HIST_ANGLE,    0, 6, row);
    PrintHistogram(dim_32, M->histograms_DISP[2], HistImgH,      HIST_HUE,      1, 2, row);
    PrintHistogram(dim_32, M->histograms_DISP[3], HistImgS,      HIST_SAT,      1, 3, row);
    PrintHistogram(dim_32, M->histograms_DISP[4], LBPhist_img,   LBP_HIST,      0, 7, row);
    PrintHistogram(dim_32, M->histograms_DISP[5], iichist_img,   iic_HIST,      0, 4, row);
  }
}

static inline void UpdateHistogram(CvHistogram* H1, CvHistogram*H2, int dim){

    cvNormalizeHist(H1, 1);
    cvNormalizeHist(H2, 1);

    for(int i = 0; i < dim; i++ )
    {

        double val1 =  cvGetReal1D(H1->bins,i);
        double val2 =  cvGetReal1D(H2->bins,i);

        //double max = val1 > val2 ? val1 : val2;
        //if (val2 > val1){

        double update = 0.2*val1 + 0.8*val2;
        //double update = val2;
        //update  = update > 0.5 ? 0.5 : update;
        cvSetReal1D(H1->bins, i, update);
        //}

    }
}

void FeatureAnalysis(Features *F, Model* M, Statistics *S, IplImage *gbs, bool dynamic)
{
    cvInRangeS( F->hsv, cvScalar(0,   smin, min(vmin,vmax), 0),
                cvScalar(180, 256,  max(vmin,vmax), 0), bin[0] );

    //const static int Ehist_size = 9;
    //const static int Chist_size = 32;
    //const static int Lhist_size = 32;

    for(int i=0;i<S->no_features;i++){
        S->L1[i] = S->L1[i] > 0 ? S->L1[i] : 0.01;
        S->L0[i] = S->L0[i] > 0 ? S->L0[i] : 0.01;
        S->gmax[i] = S->gmax[i] > 0 ? S->gmax[i] : 1.;
    }

    for(int i = 0; i < S->nos; i++ )
    {
        cvCmpS(gbs, S->gray_id[i], mask[0], CV_CMP_EQ); //mask out image segmen

        cvAnd(bin[0], mask[0], bin[1]);

        cvCalcHist(&F->mag,   S->H_SF[0], 0, mask[0]);
        cvCalcHist(&F->ang32, S->H_SF[1], 0, mask[0]); //NULL=mask
        cvCalcHist(&F->hue,   S->H_SF[2], 0, bin[1]); //NULL=mask
        cvCalcHist(&F->sat,   S->H_SF[3], 0, mask[0]);
        cvCalcHist(&F->lbp,   S->H_SF[4], 0, mask[0]); //NULL=mask
        cvCalcHist(&F->iic,   S->H_SF[5], 0, mask[0]);

        for(int j=0;j<S->no_features;j++)
        {
            S->G_score[j] = Gstat(M->histograms[j], S->H_SF[j], M->dim[j]);

            //0.9*exp(-lambda*G_scoreV);
            S->P_FgGt[S->no_features*i+j] = EXP_DIST_1(S->Z1[j], S->L1[j], S->G_score[j]);
            //1. - S->P_VgGt[i];
            S->P_FgGf[S->no_features*i+j] = EXP_DIST_0(S->Z0[j], S->L0[j], S->G_score[j]);

            cvSet(PG[j], cvScalar( S->G_score[j] ), mask[0]);

            cvSet(F->P_X1[j], cvScalar(S->P_FgGt[S->no_features*i+j]/(S->P_FgGt[S->no_features*i+j]+S->P_FgGf[S->no_features*i+j])), mask[0]);

            cvSet(F->P_X0[j], cvScalar(S->P_FgGf[S->no_features*i+j]/(S->P_FgGt[S->no_features*i+j]+S->P_FgGf[S->no_features*i+j])), mask[0]);
        }


    }


    if (p.debug){
        static int row = 2;
        DISPLAY_IMAGE_XY(p.refresh, F->P_X1[0], 5, row); //MAG
        DISPLAY_IMAGE_XY(p.refresh, F->P_X1[1], 6, row); //ANGLE
        DISPLAY_IMAGE_XY(p.refresh, F->P_X1[2], 2, row); //IMG_HUE
        DISPLAY_IMAGE_XY(p.refresh, F->P_X1[3], 3, row); //IMG_SAT
        DISPLAY_IMAGE_XY(p.refresh, F->P_X1[4], 7, row); //LBP_ANALYSIS
        DISPLAY_IMAGE_XY(p.refresh, F->P_X1[5], 4, row); //IIC_ANALYSIS
    }

    if(dynamic)
    {
      static bool flag(false);
      for(int i = 0; i < S->no_features; i++)
      {
        if(flag)
        {
          cvAddWeighted(PG[0], 0.5, PG_prev[0], 0.5, 0, PG[0]);
        }

        cvCopy(PG[0], PG_prev[0], NULL);
      }
      flag = true;
    }
}

void ProbAnalysis2(Features *F, Statistics* S, IplImage* gbs){

  const double smallprob(1e-4);

    double t = (double)p.log_post_thres_zero_position*0.1 - 30;

    for(int i = 0; i < S->nos; i++ )
    {
        cvCmpS(gbs, S->gray_id[i], mask[0], CV_CMP_EQ); //mask out image segment

        for(int j=0; j < S->no_features; j++)
        {
          S->P_FgGt[S->no_features*i+j] = S->P_FgGt[S->no_features*i+j] > smallprob ? S->P_FgGt[S->no_features*i+j] : smallprob;

        }

        double post1 = (S->priorTrue[i])*
                       (S->P_FgGt[S->no_features*i+0])*
                       (S->P_FgGt[S->no_features*i+1])*
                       (S->P_FgGt[S->no_features*i+2])*
                       (S->P_FgGt[S->no_features*i+3])*
                       (S->P_FgGt[S->no_features*i+4])*
                       (S->P_FgGt[S->no_features*i+5]);

        double post0 = (S->priorFalse[i])*
                       (S->P_FgGf[S->no_features*i+0])*
                       (S->P_FgGf[S->no_features*i+1])*
                       (S->P_FgGf[S->no_features*i+2])*
                       (S->P_FgGf[S->no_features*i+3])*
                       (S->P_FgGf[S->no_features*i+4])*
                       (S->P_FgGf[S->no_features*i+5]);


        S->P_GtgF[i] = post1 / (post1 + post0);
        S->P_GfgF[i] = post0 / (post1 + post0);

        cvSet(F->post1, cvScalar(S->P_GtgF[i]), mask[0]);
        cvSet(F->post0, cvScalar(S->P_GfgF[i]), mask[0]);


        cvSet(F->post_ratio, cvScalar(log(post1/post0)), mask[0]);
    }


    cvCalcHist( &F->post_ratio, Ghist, 0, 0);

    if(p.debug)
    {
      const char * POST_RATIO       = "Log Post Ratio";
      PrintHistogram(256, Ghist, GhistImg, POST_RATIO, 0, 0, 3);
      if(p.refresh)
      {
        cvCreateTrackbar("T", POST_RATIO, &p.log_post_thres_zero_position, 600, NULL );
      }
    }

    cvThreshold(F->post_ratio, F->bin_class_result, t, 255, CV_THRESH_BINARY);
    cvNormalize(F->post_ratio, inv_prob[0], 0, 1, CV_MINMAX);
}

void GetMeanHist(int dim, CvHistogram* H, double mean, double max){

    cvNormalizeHist(H,1);

    for(int i = 0; i < dim; i++ )
    {
        float* val = cvGetHistValue_1D(H,i);
        //increment the mean value
        mean += i*val[0];
        if(val[0]>0.){
            max = i;
        }
    }

}

void UpdateParams(IplImage* T, Statistics *S, Features *F, bool dynamic){

    //initializations
    static int GmaxInt;
    static int L1Int;
    static int L0Int;
    static int dim_40 = 40;
    static bool weighted_mean=0;
    static bool truncated = 1;

    static double gmax=5;
    static double mean1=0.;
    static double mean0=0.;
    static CvScalar mean;
    const static double min_mean1 = 0.025;
    //const static double min_mean1 = 0.05;
    const static double max_mean1 = 10;
    //const static double min_mean1 = 0.01;
    //const static double max_mean1 = 1000;
    //const static double min_mean0 = 0.01;
    //const static double max_mean1 = 100;
    const static double min_mean0 = 0.25;

    cvCopy(T, bin[4],0);
    cvNot(T,bin[3]);
    cvNormalize(F->post_ratio,inv_prob[0], 0, 1, CV_MINMAX);
    cvAbsDiffS(inv_prob[0], inv_prob[3], cvScalar(1));

    int loop(-1);
    if(dynamic) loop=200;

    static int j=0;

    int acc = 0;
    if( j < loop)
    {
      ++j;
      if(j > 1) acc=1;
    }
    else
    {
      acc=0;
      j=0;
    }

    for(int i=0;i<S->no_features;i++)
    {

        //cvCopyHist(PG[i], &PG_DISP[i]);
        //cvNormalizeHist
        //if(i==5){
        cvCopy(PG[i],temp[5],0);
        cvNormalize(temp[5],temp[5],0,39,CV_MINMAX);
        cvCalcHist( &temp[5], Ghist2, 1, NULL);
        cvCopyHist(Ghist2, &Ghist2DISP);
        if(p.debug)
        {
          PrintGHistogram(40, Ghist2DISP, GhistImg2, "Gstat", 0, 7, 3); 
        }
        //}
        //cvWaitKey();

        cvMinMaxLoc( PG[i], NULL, &gmax, NULL, NULL, NULL);
        //printf("gmax = %0.2f",gmax);
        //smooth the change of gmax and limit minimum gmax
        if(dynamic){
            S->gmax[i] = (0.5*gmax + 0.5*S->gmax[i]) > 0 ? (0.5*gmax + 0.5*S->gmax[i]) : 0.5;
        }else{
            S->gmax[i] = gmax > 0 ? gmax : 0.5;
        }

        if(weighted_mean){
            //cvMul(F->post1,PG[i],inv_prob[1],1);
            cvMul(inv_prob[0],PG[i],inv_prob[1],1);
            cvAvgSdv(inv_prob[1], &mean, NULL, NULL);
        }else {
            cvAvgSdv(PG[i], &mean, NULL, bin[4]);
        }
        mean.val[0] = mean.val[0] > min_mean1 ? mean.val[0] : min_mean1;
        mean.val[0] = mean.val[0] > max_mean1 ? max_mean1 : mean.val[0];
        if (mean1==0.) mean1 = mean.val[0];
        if(dynamic){
            mean1 = 0.5*mean1 + 0.5*mean.val[0];
        }else{
            mean1 = mean.val[0];
        }
        if(truncated){
            if(mean1<(gmax/2)) mean1 = mean1 - gmax/(exp(gmax/mean1)-1);
            if(mean1>=(gmax/2)) mean1 = (gmax/2);
        }

        S->L1[i] = (1./(mean1)) > 0 ? (1./(mean1)) : 1 ;

        //printf("%0.4f\n", mean.val[0]);
        cvSet( hist_temp, cvScalarAll(255), 0 );
        cvCalcHist( &PG[i], S->H_G1[i], acc, bin[4]);

        cvCopyHist(S->H_G1[i], &S->H_G1_DISP[i]);


        /*
        GetMeanHist(dim_40, S->H_G0_DISP[i], mean_trial, gmax);
        S->gmax[i] = cvRound(0.5*gmax + 0.5*S->gmax[i]) > 0 ? cvRound(0.5*gmax + 0.5*S->gmax[i]) : 0.5;
        mean_trial = mean_trial > min_mean1 ? mean_trial : min_mean1;
        mean_trial = mean_trial > max_mean1 ? max_mean1 : mean_trial;
        S->L1[i] = cvRound( ((1./(mean_trial))*1000.) ) > 0 ? cvRound( ((1./(mean_trial))*1000.) ) : 1;
        //draw_dist(40, S->gmax[i], S->Z1[i], S->Z0[i], I2D(S->L1[i]), I2D(S->L0[i]), PG1_DISP[i], 1);
        */
        draw_dist(40, S->gmax[i], S->Z1[i], S->L1[i], hist_temp, 1);
        cvAddWeighted( hist_temp, 0.5, PG1_DISP[i], 0.5, 0, PG1_DISP[i] );


        if(p.debug)
        {
          for(int i=0;i<S->no_features;i++)
          {
            GmaxInt = cvRound(S->gmax[i]);
            L1Int = cvRound(S->L1[i]*1000);
            L0Int = cvRound(S->L0[i]*1000);

            if(p.refresh)
            {
                cvNamedWindow(P1name[i]);
                cvCreateTrackbar( "L1", P1name[i], &L1Int, 40000, NULL );
                cvCreateTrackbar( "L0", P1name[i], &L0Int, 4000, NULL );
                cvCreateTrackbar( "gmax", P1name[i], &GmaxInt, 60, NULL );
            }
          }
        }


        cvSet( hist_temp, cvScalarAll(255), 0 );

        if(weighted_mean){
            cvMul(inv_prob[3],PG[i],inv_prob[2],1);
            cvAvgSdv(inv_prob[2], &mean, NULL, NULL);
        }else{
            cvAvgSdv(PG[i], &mean, NULL, bin[3]);
        }

        mean.val[0] = gmax-mean.val[0];
        mean.val[0] = mean.val[0] > min_mean0 ? mean.val[0] : min_mean0;

        if (mean0==0.) mean0 = mean.val[0];
        if(dynamic){
            mean0 = 0.5*mean0 + 0.5*mean.val[0];
        }else{
            mean0 = mean.val[0];
        }
        if(truncated){
            if(mean0<(gmax/2)) mean0 = mean0 - gmax/(exp(gmax/mean0)-1);
            if(mean0>=(gmax/2)) mean0 = (gmax/2);
        }
        //S->L0[i] = cvRound( ((1./(mean0))*1000.) ) > 0 ? cvRound( ((1./(mean0))*1000.) ) : 1;
        S->L0[i] = (1./(mean0)) > 0 ? (1./(mean0)) : 1 ;
        //printf("%d\n", S->L0[i]);
        cvCalcHist( &PG[i], S->H_G0[i], acc, bin[3]);

        cvCopyHist(S->H_G0[i], &S->H_G0_DISP[i]);

        /*
        GetMeanHist(dim_40, S->H_G0_DISP[i],mean_trial, gmax);
        S->gmax[i] = cvRound(0.5*gmax + 0.5*S->gmax[i]) > 0 ? cvRound(0.5*gmax + 0.5*S->gmax[i]) : 0.5;
        mean_trial = mean_trial > min_mean0 ? mean_trial : min_mean0;
        S->L0[i] = cvRound( ((1./(mean_trial))*1000.) ) > 0 ? cvRound( ((1./(mean_trial))*1000.) ) : 1;
        */
        //printf("mean=%0.4f\n", mean_trial);

        draw_dist(40, S->gmax[i], S->Z0[i], S->L0[i], hist_temp, 0);

        cvAddWeighted( hist_temp, 0.5, PG0_DISP[i], 0.5, 0, PG0_DISP[i] );

        if (p.debug){ /*cvShowImage( P0name[i], PG0_DISP[i]);*/ }

        S->Z1[i] = (1/S->L1[i])*(1-exp(-S->L1[i]*S->gmax[i]));
        S->Z0[i] = (1/S->L0[i])*(exp(S->L0[i]*S->gmax[i])-1);

    }

    if (p.debug){
      PrintHistogram(dim_40, S->H_G1_DISP[0], PG1_DISP[0], P1name[0], 0, 5, 3);
      PrintHistogram(dim_40, S->H_G1_DISP[1], PG1_DISP[1], P1name[1], 0, 6, 3);
      PrintHistogram(dim_40, S->H_G1_DISP[2], PG1_DISP[2], P1name[2], 0, 2, 3);
      PrintHistogram(dim_40, S->H_G1_DISP[3], PG1_DISP[3], P1name[3], 0, 3, 3);
      PrintHistogram(dim_40, S->H_G1_DISP[4], PG1_DISP[4], P1name[4], 0, 7, 3);
      PrintHistogram(dim_40, S->H_G1_DISP[5], PG1_DISP[5], P1name[5], 0, 4, 3);
    }
}

bool CheckConvergence(Statistics* S, int em){

    static double * prev1;
    static double * prev0;
    //static bool * converged;
    const static double tolerance = 0.005;
    int count = 0;

    if(em == 0){
        prev1 = new double[S->no_features];
        prev0 = new double[S->no_features];
        for (int i=0; i<S->no_features; i++){
        prev1[i] = 0.;
        prev0[i] = 0.;
        //converged[i] = 0;
        }

    }
    else{

        for (int i=0; i<S->no_features; i++){

            //if ( S->L1[i] == prev1[i] && S->L0[i] == prev0[i]) count ++;
            if( ((S->L1[i] <= prev1[i] + tolerance) && (S->L1[i] >= prev1[i] - tolerance)) &&
            ((S->L0[i] <= prev0[i] + tolerance) && (S->L0[i] >= prev0[i] - tolerance))
            ) count++;

            prev1[i] = S->L1[i];
            prev0[i] = S->L0[i];
        }

    }
    if(count == S->no_features) return 1;
        else return 0;
}

void FindObstacleBoundary(IplImage *Out)
{
  unsigned char * OutPixelData = (unsigned char *)(Out->imageData);

  static int y = (Out->height) - 1;

  bool flag0 = 0, flag1 = 0;

  for(int x = 0; x < cvRound((Out->width)/2); x++)
  {
    int x0 = cvRound((Out->width)/2) - (x+1);
    int x1 = cvRound((Out->width)/2) + (x+1);

    int Out_index_0 = (y*Out->width+x0)*Out->nChannels;
    int Out_index_1 = (y*Out->width+x1)*Out->nChannels;

    if(flag0 == 1 || OutPixelData[Out_index_0+0] != 255)
    {
      flag0 = 1;
      OutPixelData[Out_index_0+0] = 0;
    }
    if(flag1 == 1 || OutPixelData[Out_index_1+0] != 255)
    {
      flag1 = 1;
      OutPixelData[Out_index_1+0] = 0;
    }
  }

  for(int x = 0; x < Out->width; x++)
  {
    int flag = 0;
    for (int y = 0; y < Out->height; y++)
    {
      int y1 = (Out->height) - (y+1);

      int Out_index = (y1*Out->width+x)*Out->nChannels;

      if(flag == 1 || OutPixelData[Out_index+0] != 255)
      {
          flag = 1;
          OutPixelData[Out_index+0] = 0;
      }
    }
  }
}

void ExtractBoundary(CvSize S, Boundary *B)
{
  unsigned char * OutPixelData = (unsigned char *)(B->Bimg->imageData);

  for(int x = 0; x < S.width; ++x)
  {
    for(int y = 0; y < S.height; ++y)
    {
      int yreverse = S.height - y - 1;

      int pixelIndex = (yreverse * S.width + x);

      if(OutPixelData[pixelIndex] != 255)
      {
        B->boundaryPixels[x] = cvPoint(x, y);
        break;
      }
    }
  }
}

void CalculateDistances(CvSize S, Boundary *B, CamCalib camera, BotCalib bot)
{
  if(p.debug) cvZero(DEPTH_MAP);

  for(int i = 0; i < S.width; ++i)
  {
    double px = (double)(B->boundaryPixels[i].x);
    double py = (double)(B->boundaryPixels[i].y);

    double beta = camera.pan_angle + atan((camera.cx - px) * camera.pix_size_x / camera.fx);

    double alpha = atan((camera.cy - py) * camera.pix_size_y / camera.fy);
    alpha = alpha > 0 ? alpha : 0.0001;

    double depth = camera.height_from_ground / tan(camera.tilt_angle + alpha);
    depth = depth < 300000 ? depth : 300000;

    B->obstaclePolarCoord[i].x = (float)depth; // Distance.
    B->obstaclePolarCoord[i].y = beta; // Angle w.r.t heading direction of robot.

    if(p.debug)
    {
      POLAR2CART(&B->obstaclePolarCoord[i], &B->obstacleCartesianCoord[i]);
      //fit representation to 400x400 pixel image
      CART2DISPLAY(&B->obstacleCartesianCoord[i], &B->obstacleDisplayCoord[i], &B->RobotDisplayCoord, bot.display_scale_factor);
      cvCircle(DEPTH_MAP, B->obstacleDisplayCoord[i], 1, CV_RGB(255,0,0), 1, 8, 0);
    }
  }

  if(p.debug)
  {
    // draw robot
    cvCircle(DEPTH_MAP, B->RobotDisplayCoord, 2, CV_RGB(255,0,0), 2, 8, 0);

    // draw min distance to obstacles
    cvCircle(DEPTH_MAP, B->RobotDisplayCoord, bot.min_dist_to_obst * bot.display_scale_factor, CV_RGB(255,255,0), 2, 8, 0);

    // draw side lines
    cvLine(DEPTH_MAP, B->RobotDisplayCoord, B->obstacleDisplayCoord[0], CV_RGB(0,200,0),1,8,0);
    cvLine(DEPTH_MAP, B->RobotDisplayCoord, B->obstacleDisplayCoord[S.width-1], CV_RGB(0,200,0),1,8,0);

    // draw circles marking equidistant points.
    for(int i=0; i<10; i++)
    {
      cvCircle(DEPTH_MAP,  B->RobotDisplayCoord, i*30, CV_RGB(0,0,200), 1, 8, 0);
    }
  }
}

void InterpretDepthArray(CvSize S, Boundary *B, BotCalib bot)
{
  // A segment is a group of adjacent frontier points (without breaks) that are further that a specified distance threshold.
  int segmentCount=1, frontierCount=0;

  bool currentFlag = false, previousFlag = false;
  const int maxSegments = (S.width+1)/2;
  std::vector<int> segmentHistogram(maxSegments,0);

  // Loop over all columns in the image.
  for(int i=0; i < S.width; ++i)
  {
    // If the Euclidean distance to the boundary in that column is greater than the minimum distance threshold.
    if(abs(B->obstaclePolarCoord[i].x) > bot.min_dist_to_obst)
    {
      currentFlag = true; // Start a new segment of frontiers.

      B->frontiers[frontierCount].x = B->obstaclePolarCoord[i].x;
      B->frontiers[frontierCount].y = B->obstaclePolarCoord[i].y;
      B->frontiers[frontierCount].z = segmentCount;
      ++frontierCount;

      segmentHistogram[segmentCount] += 1;
    }
    else
    {
      currentFlag = false;
    }

    if(!currentFlag && previousFlag) ++segmentCount;

    previousFlag = currentFlag;
  }

  B->angleBetweenBestObstacleFreeSegment = 0.0;
  int best_segment = 0;
  int firstFrontierIndex = 0;
  for(int s = 1; s <= segmentCount; ++s)
  {
    if(segmentHistogram[s])
    {
      firstFrontierIndex += segmentHistogram[s - 1];
      int lastFrontierIndex = firstFrontierIndex + segmentHistogram[s] - 1;

      //For each segment, find the median angle in the frontiers.
      //This will be half way through the set of frontiers in each segment.
      int medianFrontierIndex = firstFrontierIndex + cvRound(segmentHistogram[s] / 2.0) - 1;

      B->Median_Angle[s].x = B->frontiers[medianFrontierIndex].y; // Angle from current heading horizontal direction
      B->Median_Angle[s].z = B->frontiers[medianFrontierIndex].x; // Euclidean Distance
      B->Median_Angle[s].y = B->frontiers[lastFrontierIndex].y - B->frontiers[firstFrontierIndex].y; // Angle between first and last frontiers in segment.

      if(B->angleBetweenBestObstacleFreeSegment < abs(B->Median_Angle[s].y))
      {
        B->angleBetweenBestObstacleFreeSegment = abs(B->Median_Angle[s].y);
        best_segment = s;
      }
    }
  }

  CvPoint2D32f distanceAndAngle;
  B->distanceToObstacleInBestHeadingDirection = distanceAndAngle.x = B->Median_Angle[best_segment].z; // Distance
  B->angleToBestHeadingDirection = distanceAndAngle.y = B->Median_Angle[best_segment].x; // Angle from vertical

  if(p.debug)
  {
    CvPoint2D32f bestMedianFrontierPoint;
    CvPoint bestMedianFrontierPointDisplay;

    POLAR2CART(&distanceAndAngle, &bestMedianFrontierPoint);
    CART2DISPLAY(&bestMedianFrontierPoint, &bestMedianFrontierPointDisplay, &B->RobotDisplayCoord, bot.display_scale_factor);

    cvLine(DEPTH_MAP, B->RobotDisplayCoord, bestMedianFrontierPointDisplay, CV_RGB(0,200,200), 1, 8, 0);

    DISPLAY_IMAGE_XY(p.refresh, DEPTH_MAP, 8, 0);
  }
}
