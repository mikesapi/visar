#ifndef VISAR_STRUCTURES
#define VISAR_STRUCTURES

#include <string>

#include "Params.h"

#include "../comms/Net.h" //header for UDP transmission
#include "OverwriteQueue.h"

/* boundary structure definition */
struct boundary
{
  IplImage *Bimg;
  CvPoint *boundaryPixels;
  CvPoint2D32f *obstaclePolarCoord;
  CvPoint2D32f *obstacleCartesianCoord;
  CvPoint *obstacleDisplayCoord;

  // The points (distance,angle,segmentId) which are further away than a specified threshold.
  CvPoint3D32f *frontiers;
  
  CvPoint RobotDisplayCoord;
  
  //steering direction
  CvPoint3D32f *Median_Angle; 
  
  double angleBetweenBestObstacleFreeSegment;
  double distanceToObstacleInBestHeadingDirection;
  double angleToBestHeadingDirection;
};
typedef struct boundary Boundary; 


#define FEATURECOUNT 6
struct statistics {
  
    int *id;
    int *size;
    int *gray_id;
    CvScalar* mean;
    CvScalar* stdDev;
    CvRect* box; 
    int no_features;


    double *priorTrue;
    double *priorFalse;
    IplImage* prior_img;

    double *P_FgGt, *P_FgGf;

    double *P_GtgF;//posterior
    double *P_GfgF;

    double *G_score;

    int nos; //nuber of segments
    int img_w;
    int img_h;

    double *L1, *L0;
    double *Z1, *Z0, *gmax;

    CvHistogram *H_SF[FEATURECOUNT];
    CvHistogram *H_G1[FEATURECOUNT], *H_G1_DISP[FEATURECOUNT];
    CvHistogram *H_G0[FEATURECOUNT], *H_G0_DISP[FEATURECOUNT];

};
typedef struct statistics Statistics;

struct model {

    //CvRect* box; //safe area
    IplImage* safeRegionMask; //safe area
    
    //Basic Stats
    CvScalar* mean;
    CvScalar* stdDev;

    CvHistogram *histograms[FEATURECOUNT], *histograms_DISP[FEATURECOUNT];
    int *dim;
};
typedef struct model Model;


struct features {

    //Edges
    IplImage *mag, *ang32, *P_ang;

    //Colour
    IplImage *hsv, *lab, *YCrCb;
    IplImage *hue, *sat, *val; 
    IplImage *Cr, *a, *Cb, *iic;
    IplImage *P_hue, *P_sat, *P_val;

    //LBP
    IplImage *lbp, *P_lbp;

    //Posterior
    IplImage *post0, *post1, *post_ratio, *P_X1[FEATURECOUNT], *P_X0[FEATURECOUNT];

    //Classification
    IplImage *bin_class_result;

};
typedef struct features Features;
#undef FEATURECOUNT

#endif
