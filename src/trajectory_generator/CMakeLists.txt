cmake_minimum_required(VERSION 2.8.3)
project(trajectory_generator)

set(CMAKE_VERBOSE_MAKEFILE "false")
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

set(ADDITIONAL_CXX_FLAG "-Wall -O3 -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_CXX_FLAG}")

find_package(catkin REQUIRED COMPONENTS
  roscpp
  tf
  std_msgs
  nav_msgs
  visualization_msgs
  quadrotor_msgs
  plan_env
  bspline_opt
  ego_planner
)
catkin_package(
  INCLUDE_DIRS include
 LIBRARIES trajectory_generator
 CATKIN_DEPENDS plan_env bspline_opt ego_planner
#  DEPENDS system_lib
)

find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)
set(Eigen3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})

include_directories(
    include
    SYSTEM 
    third_party
  	${catkin_INCLUDE_DIRS}
  	${Eigen3_INCLUDE_DIRS}
    ${PCL_INCLUDE_DIRS} 
)

add_executable(trajectory_generator
	src/trajectory_generator_node.cpp
	src/trajectory_generator_waypoint.cpp
  src/Astar_searcher.cpp
  src/uniform_bspline.cpp
)

target_link_libraries(trajectory_generator
   ${catkin_LIBRARIES}
   ${PCL_LIBRARIES}
)

add_executable(trajectory_server 
  src/traj_server.cpp
  src/uniform_bspline.cpp
  src/Astar_searcher.cpp)
target_link_libraries(trajectory_server ${catkin_LIBRARIES})
#add_dependencies(traj_server ${${PROJECT_NAME}_EXPORTED_TARGETS})
