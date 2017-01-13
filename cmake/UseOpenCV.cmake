find_package(OpenCV 2.4.8 REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
add_definitions(-DWITH_OPENCV)
