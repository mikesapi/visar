Implementation of the Ground Plane Segmentation algorithm as described in:

Michael Sapienza, Kenneth P. Camilleri, A generative traversability model for monocular robot self-guidance,
9th Int. Conf. on Informatics in Control, Automation and Robotics, 2012. (oral) [pdf - slides - video - BibTex]

This software was developed by Michael Sapienza within the department of 'Systems & Control Engineering (SCE)' at the University of Malta under the supervision of Prof. Kenneth Camilleri.

What to do to run this software (tested on Ubuntu Linux 14.04LTS):

* Open a terminal window

* Install the OPENCV Library (2.4.8)
> sudo apt-get install libopencv-dev

* Navigate to the directory in which you want to clone the repo.
> cd software/

* Clone the repo
> git clone git@github.com:mikesapi/visar.git

* Build the code with cmake
> cd visar
> mkdir build && cd build
> cmake ..
> make

* Run the executable on video input, and in debug mode
> ./bin/visar -i ../data/eng_stat_obst.avi -d


Contributers:
Michael Sapienza
Kenneth Camilleri
Kenneth Scerri
Cristian Roman
