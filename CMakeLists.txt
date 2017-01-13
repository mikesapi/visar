cmake_minimum_required(VERSION 2.7)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -pthread")

project(visar)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/build/bin/")

set(targetname visar)

include(${PROJECT_SOURCE_DIR}/cmake/UseOpenCV.cmake)

set(main_source
main.cpp
)

set(core_sources
core/capture.cpp
core/OverwriteQueue.cpp
core/init_structures.cpp
core/img_proc_fcns.cpp
core/Util.cpp
core/Params.cpp
core/ProgramInputParser.cpp
core/SegmentationPipeline.cpp
core/RobotController.cpp
core/DummyRobotController.cpp
core/RoverRobotController.cpp
core/Timer.cpp
)

set(core_headers
core/capture.h
core/OverwriteQueue.h
core/CommunicationLoop.h
core/structures.h
core/init_structures.h
core/img_proc_fcns.h
core/Util.h
core/Params.h
core/ProgramInputParser.h
core/SegmentationPipeline.h
core/botCalib.h
core/camCalib.h
core/RobotController.h
core/DummyRobotController.h
core/RoverRobotController.h
core/Timer.h
)


set(comms_sources
comms/Net.cpp
)

set(comms_headers
comms/Net.h
)

set(graphbasedsegmentation_headers
graphbasedsegmentation/convolve.h
graphbasedsegmentation/disjoint-set.h
graphbasedsegmentation/filter.h
graphbasedsegmentation/image.h
graphbasedsegmentation/imconv.h
graphbasedsegmentation/imutil.h
graphbasedsegmentation/misc.h
graphbasedsegmentation/segment-graph.h
graphbasedsegmentation/segment-image.h
)

set(sources
${main_source}
${core_sources}
${comms_sources}
${graphbasedsegmentation_sources}
)

set(headers
${headers}
${core_headers}
${comms_headers}
${graphbasedsegmentation_headers}
)

add_executable(${targetname} ${sources} ${headers})

target_link_libraries(${targetname} pthread)
include(${PROJECT_SOURCE_DIR}/cmake/LinkOpenCV.cmake)
