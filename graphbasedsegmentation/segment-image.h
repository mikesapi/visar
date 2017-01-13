/*
Copyright (C) 2006 Pedro Felzenszwalb

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/

#ifndef SEGMENT_IMAGE
#define SEGMENT_IMAGE

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <cmath>
using namespace std;

#include "image.h"
#include "misc.h"
#include "filter.h"
#include "segment-graph.h"

// random color
rgb random_rgb(){ 
  rgb c;
  
  c.r = (uchar)random();
  c.g = (uchar)random();
  c.b = (uchar)random();

  return c;
}

// dissimilarity measure between pixels
static inline float diff(image<float> *r, image<float> *g, image<float> *b,
			 int x1, int y1, int x2, int y2) {
  return sqrt(square(imRef(r, x1, y1)-imRef(r, x2, y2)) +
	      square(imRef(g, x1, y1)-imRef(g, x2, y2)) +
	      square(imRef(b, x1, y1)-imRef(b, x2, y2)));
}

static inline float diff2(image<float> *r, 
			 int x1, int y1, int x2, int y2) {
  return sqrt(square(imRef(r, x1, y1)-imRef(r, x2, y2)) //+
	      //square(imRef(g, x1, y1)-imRef(g, x2, y2)) +
	      //square(imRef(b, x1, y1)-imRef(b, x2, y2))
	      );
}

/*
 * Segment an image
 *
 * Returns a color image representing the segmentation.
 *
 * im: image to segment.
 * sigma: to smooth the image.
 * c: constant for treshold function.
 * min_size: minimum component size (enforced by post-processing stage).
 * num_ccs: number of connected components in the segmentation.
 */

//image<float> *r;

void init_r(CvSize S){
 //r = new image<float>( S.width, S.height);
  //printf("Arrived here!\n");
}
//image<float> *smooth_r;

int* comp_arr;
int* segments;
uchar* graycolors;

void SegmentImage(IplImage* output, IplImage *im, float sigma, float c, int min_size,
			  int *num_ccs) {
  static int width = im->width;
  static int height = im->height;

  image<float> *r = new image<float>(width, height);
  //image<float> *g = new image<float>(width, height);
  //image<float> *b = new image<float>(width, height);

	char * pixelData = (char *)(im->imageData);
	//char * im_data = (char *)imPtr(input, 0, 0);


	  for (int y = 0; y < height; y++) {
    		for (int x = 0; x < width; x++) {
			int ind = (y*im->width+x)*im->nChannels;
      			imRef(r, x, y) = pixelData[ind+0];
    			}
  		}

  // smooth each color channel  
/*  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      imRef(r, x, y) = imRef(im, x, y).r;
      //imRef(g, x, y) = imRef(im, x, y).g;
      //imRef(b, x, y) = imRef(im, x, y).b;
    }
  }
*/
  image<float> *smooth_r = smooth(r, sigma);
  //image<float> *smooth_g = smooth(g, sigma);
  //image<float> *smooth_b = smooth(b, sigma);
  delete r;
  //delete g;
  //delete b;
 
  // build graph
  edge *edges = new edge[width*height*4];
  int num = 0;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      if (x < width-1) {
	edges[num].a = y * width + x;
	edges[num].b = y * width + (x+1);
	//edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y);
	edges[num].w = diff2(smooth_r, x, y, x+1, y);
	num++;
      }

      if (y < height-1) {
	edges[num].a = y * width + x;
	edges[num].b = (y+1) * width + x;
	//edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y);
	edges[num].w = diff2(smooth_r, x, y, x, y+1);
	num++;
      }

      if ((x < width-1) && (y < height-1)) {
	edges[num].a = y * width + x;
	edges[num].b = (y+1) * width + (x+1);
	//edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y);
	edges[num].w = diff2(smooth_r, x, y, x+1, y+1);
	num++;
      }

      if ((x < width-1) && (y > 0)) {
	edges[num].a = y * width + x;
	edges[num].b = (y-1) * width + (x+1);
	//edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y);
	edges[num].w = diff2(smooth_r, x, y, x+1, y-1);
	num++;
      }
    }
  }
  delete smooth_r;
  //delete smooth_g;
  //delete smooth_b;

  // segment
  universe *u = segment_graph(width*height, num, edges, c);
  
  // post process small components
  for (int i = 0; i < num; i++) {
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
      u->join(a, b);
  }
  delete [] edges;
  *num_ccs = u->num_sets();
	int length = u->num_sets();

  //image<rgb> *output = new image<rgb>(width, height);
//  IplImage* output = cvCreateImage( cvGetSize(im), IPL_DEPTH_8U, 1);


 // pick random colors for each component
//  rgb *colors = new rgb[width*height];
//  rgb *avg_colors = new rgb[width*height];
/*   for (int i = 0; i < width*height; i++){
    colors[i] = random_rgb();
    avg_colors[i] = random_rgb();
    }
*/
char * outputPixelData = (char *)(output->imageData); 
int comp1=0;
comp_arr = new int[width*height];
//int comp_arr[width*height];
int i = 0;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int comp = u->find(y * width + x);
//	avg_colors[1].r = comp;
//	avg_colors[1].g = comp;
//	avg_colors[1].b = comp;
	//uchar label = comp;

	if (comp != comp1){
	//printf("comp=%d, col=%d	", comp, avg_colors[1].r);
	comp_arr[i] =comp;
	i++;
	}
      	//imRef(output, x, y) = avg_colors[1];
	//imRef(output, x, y) = avg_colors[comp];

	//char * im_data = (char *)imPtr(input, 0, 0);


	 /*
			int index = (y*output->width+x)*output->nChannels;
      			//outputPixelData[index+0] = cvRound((comp/70000)*255);
			//outputPixelData[index+0] = avg_colors[1].r;
			outputPixelData[index+0] = comp;
*/
	
	comp1=comp;

    }
  }  


cvNormalize(output,output,0,255,CV_MINMAX,0);

//IplImage *output_img = cvCreateImage(cvGetSize(im), IPL_DEPTH_8U, 1);

//cvConvertScaleAbs(output, output_img, 255, 0);

//cvNamedWindow("Trial");
//cvShowImage("Trial", output_img);

vector<int> myvector (comp_arr, comp_arr+i);
vector<int>::iterator it;
sort (myvector.begin(), myvector.begin()+i);


 // print out content:
/*  cout << "myvector contains:";
  for (it=myvector.begin(); it!=myvector.end(); ++it)
    cout << " " << *it;

  cout << endl;
*/
  // using default comparison:
  it = unique (myvector.begin(), myvector.end()); 
                                                 

  myvector.resize( it - myvector.begin() );


segments = new int[length];
//int segments[length];

 // print out content:
  //cout << "myvector contains:";
i=0;
  for (it=myvector.begin(); it!=myvector.end(); ++it){
	segments[i] = *it;
    //cout << " " << segments[i];
	i++;
	}
  //cout << endl;
	//cout << endl;

	//cout << "vector size = " << myvector.size() << " num_ccs =" << *num_ccs << endl; 

//uchar graycolors[width*height];
graycolors = new uchar[width*height];
for(size_t i = 0; i < myvector.size(); i++){
graycolors[ segments[i] ] = i*(floor(255/myvector.size()));
//printf("col=%d, seg=%d\n", graycolors[segments[i]], segments[i]);

}

//delete[] segments;

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int comp = u->find(y * width + x);
 
			int index = (y*output->width+x)*output->nChannels;

			outputPixelData[index+0] = graycolors[comp];
			//printf("col= %d ", graycolors[comp]);


    }
  }  

//  delete [] colors;  
delete u;
	//cvReleaseImage(&output);
delete [] comp_arr;
comp_arr = NULL;

delete [] segments;
segments = NULL;

delete [] graycolors;
graycolors = NULL;

//~myvector();
}

#endif
