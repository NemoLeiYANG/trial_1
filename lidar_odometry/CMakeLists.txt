#v1************************************************************************************
# This is developed at Carnegie Mellon University in collaboration with Autel Robotics 
#                                                                                      
# PI:                                                                                  
# George Kantor                                                                        
#                                                                                      
# Authors:                                                                              
# Weizhao Shao                                                                         
# Cong Li                                                                              
# Srinivasan Vijayarangan                                                              
#                                                                                      
# Please refer to the contract document for details on license/copyright information.  
#**************************************************************************************
cmake_minimum_required(VERSION 2.8.3)
project(lidar_odometry)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x -g")

SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")

#GTSAM_USE_QUATERNIONS

find_package(catkin REQUIRED COMPONENTS
  geometry_msgs
  nav_msgs
  sensor_msgs
  roscpp
  rospy
  std_msgs
  tf)

find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)
find_package(OpenCV REQUIRED)
find_package(Ceres REQUIRED)
find_package( GTSAMCMakeTools )
find_package( GTSAM REQUIRED )
find_package(Boost REQUIRED COMPONENTS system)
find_package(yaml-cpp REQUIRED )

link_directories(/home/cmuautel/Documents/yaml-cpp/build/CMakeFiles/Export/lib)
include_directories(/home/cmuautel/Documents/yaml-cpp/include/yaml-cpp)
include_directories(
  include
	${catkin_INCLUDE_DIRS} 
	${EIGEN3_INCLUDE_DIR} 
	${PCL_INCLUDE_DIRS})

catkin_package(
  CATKIN_DEPENDS geometry_msgs nav_msgs roscpp rospy std_msgs
  DEPENDS EIGEN3 PCL OpenCV
  INCLUDE_DIRS include
)


add_executable(scan_to_matching src/scan_to_map_matching.cpp)
target_link_libraries(scan_to_matching ${catkin_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${CERES_LIBRARIES} gtsam yaml-cpp)

add_executable(lidarScanRegistration src/lidarScanRegistration.cpp)
target_link_libraries(lidarScanRegistration ${catkin_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${CERES_LIBRARIES} yaml-cpp)




