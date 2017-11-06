cmake_minimum_required(VERSION 3.7)
project(pixy_roimux)

set(CMAKE_CXX_STANDARD 11)

list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT COMPONENTS Geom REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

set(Eigen3_DIR $ENV{Eigen3_CMAKEDIR})
find_package(Eigen3 REQUIRED CONFIG)
include_directories(${EIGEN3_INCLUDE_DIRS})

set(RapidJSON_DIR $ENV{RapidJSON_CMAKEDIR})
find_package(RapidJSON REQUIRED CONFIG)
include_directories(${RapidJSON_INCLUDE_DIRS})

find_library(GENFIT_LIBRARIES NAMES libgenfit2.so PATHS $ENV{GENFIT}/lib)
include_directories($ENV{GENFIT}/include)

set(SOURCE_FILES main.cpp ChargeData.cpp ChargeData.h RunParams.cpp RunParams.h Event.cpp Event.h ChargeHits.cpp ChargeHits.h NoiseFilter.cpp NoiseFilter.h PrincipalComponentsCluster.cpp PrincipalComponentsCluster.h KalmanFit.cpp KalmanFit.h)
add_executable(pixy_roimux ${SOURCE_FILES})

set(DEBUG_SOURCE_FILES debug.cpp ChargeData.cpp ChargeData.h RunParams.cpp RunParams.h Event.cpp Event.h ChargeHits.cpp ChargeHits.h NoiseFilter.cpp NoiseFilter.h PrincipalComponentsCluster.cpp PrincipalComponentsCluster.h KalmanFit.cpp KalmanFit.h)
add_executable(pixy_debug ${DEBUG_SOURCE_FILES})

set(DISPLAY_SOURCE_FILES genfitDisplay.cpp)
add_executable(pixy_display ${DISPLAY_SOURCE_FILES})

target_link_libraries(pixy_roimux ${ROOT_LIBRARIES} ${GENFIT_LIBRARIES})
target_link_libraries(pixy_debug ${ROOT_LIBRARIES} ${GENFIT_LIBRARIES})
target_link_libraries(pixy_display ${ROOT_LIBRARIES} ${GENFIT_LIBRARIES})

#install(TARGETS pixy_roimux DESTINATION bin)
