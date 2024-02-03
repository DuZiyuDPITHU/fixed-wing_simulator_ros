#include <algorithm>
#include <fstream>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseStamped.h>
#include <iostream>
#include <math.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>
#include <quadrotor_msgs/PolynomialTrajectory.h>
#include <quadrotor_msgs/Bspline.h>
#include <random>
#include <ros/console.h>
#include <ros/ros.h>
#include <sensor_msgs/PointCloud2.h>
#include <string>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>

// Useful customized headers
#include "Astar_searcher.h"
#include "backward.hpp"
#include "trajectory_generator_waypoint.h"
#include "uniform_bspline.h"
//#include "traj_server.h"

using namespace std;
using namespace Eigen;

TrajectoryGeneratorWaypoint *_trajGene = new TrajectoryGeneratorWaypoint();
AstarPathFinder *_astar_path_finder = new AstarPathFinder();

// Set the obstacle map
double _resolution, _inv_resolution, _path_resolution;
double _x_size, _y_size, _z_size;
Vector3d _map_lower, _map_upper;
int _max_x_id, _max_y_id, _max_z_id;
int vis_traj_id = 0;

// Param from launch file
double _vis_traj_width;
double _Vel, _Acc;
int _dev_order, _min_order;
double hovering_speed, hovering_radius;
int if_2d_search;
double h_reference;

int num_cp;
double dPI = 3.1415926535897;
ros::NodeHandle* nh_ptr;

// ros related
ros::Subscriber _map_sub, _pts_sub, _odom_sub;
ros::Publisher _traj_init_vis_pub, _traj_opt_vis_pub, _traj_vis_pub, _traj_pub, _path_vis_pub, _default_traj_pub;
ros::Publisher bspline_pub;

// for planning
Vector3d odom_pt, odom_vel, start_pt, target_pt, target_vel, target_acc, start_vel, start_acc(0, 0, 0);
int _poly_num1D;
MatrixXd _polyCoeff;
VectorXd _polyTime;
double time_duration;
double time_duration_default;
ros::Time time_traj_start;
bool has_odom = false;
bool has_target = false;
bool has_default = false;
UniformBspline cur_bspline_traj;
UniformBspline cur_bspline_vel;
UniformBspline cur_bspline_acc;
UniformBspline default_traj;
UniformBspline default_vel;
UniformBspline default_acc;
int hovering_pub_count = 0;
ros::Time time_hovering_start;

// for replanning
enum STATE {
  INIT,
  WAIT_TARGET,
  GEN_NEW_TRAJ,
  EXEC_TRAJ,
  REPLAN_TRAJ
} exec_state = STATE::INIT;
double no_replan_thresh, replan_thresh;
ros::Timer _exec_timer;
void execCallback(const ros::TimerEvent &e);

// declare
enum VIS_TYPE {
  VIS_INIT, 
  VIS_OPT, 
  VIS_TRAJ
} vis_type;

void changeState(STATE new_state, string pos_call);
void printState();
void visTrajectory(bool if_default, VIS_TYPE vis_type);
void visPath(MatrixXd nodes);
void trajOptimization(Eigen::MatrixXd path);
void rcvOdomCallback(const nav_msgs::Odometry::ConstPtr &odom);
void rcvWaypointsCallback(const nav_msgs::Path &wp);
void rcvPointCloudCallBack(const sensor_msgs::PointCloud2 &pointcloud_map);
void trajPublish(MatrixXd polyCoeff, VectorXd time);
void bspline_traj_pub(UniformBspline* uniform_bspline_ptr);
bool trajGeneration();
bool defaultTrajGeneration();

Vector3d zRotate(Vector3d origin_vec, double angle);
VectorXd timeAllocation(MatrixXd Path);
Vector3d getPos(double t_cur);
Vector3d getVel(double t_cur);

// change the state to the new state
void changeState(STATE new_state, string pos_call) {
  string state_str[5] = {"INIT", "WAIT_TARGET", "GEN_NEW_TRAJ", "EXEC_TRAJ",
                         "REPLAN_TRAJ"};
  int pre_s = int(exec_state);
  exec_state = new_state;
  cout << "[" + pos_call + "]: from " + state_str[pre_s] + " to " +
              state_str[int(new_state)]
       << endl;
}

void printState() {
  string state_str[5] = {"INIT", "WAIT_TARGET", "GEN_NEW_TRAJ", "EXEC_TRAJ",
                         "REPLAN_TRAJ"};
  cout << "[Clock]: state: " + state_str[int(exec_state)] << endl;
}

void rcvOdomCallback(const nav_msgs::Odometry::ConstPtr &odom) {
  odom_pt(0) = odom->pose.pose.position.x;
  odom_pt(1) = odom->pose.pose.position.y;
  odom_pt(2) = odom->pose.pose.position.z;

  odom_vel(0) = odom->twist.twist.linear.x;
  odom_vel(1) = odom->twist.twist.linear.y;
  odom_vel(2) = odom->twist.twist.linear.z;

  has_odom = true;
}

// Control the State changes
void execCallback(const ros::TimerEvent &e) {
  static int num = 0;
  num++;
  if (num == 100) {
    printState();
    if (!has_odom)
      cout << "no odom." << endl;
    if (!has_target)
      cout << "wait for goal." << endl;
    num = 0;
  }

  switch (exec_state) {
  case INIT: {
    if (!has_odom)
      return;
    if (!has_target)
      return;
    changeState(WAIT_TARGET, "STATE");
    break;
  }

  case WAIT_TARGET: {
    if (!has_target)
    {
      if (!default_traj.getFlag())
        return;
      if (hovering_pub_count==0)
      {
        time_hovering_start = ros::Time::now();
        bspline_traj_pub(&default_traj);
        visTrajectory(true, VIS_TRAJ);
        hovering_pub_count++;
        return;
      }
      else
      {
        ros::Time time_now = ros::Time::now();
        double t_cur = (time_now - time_hovering_start).toSec();
        if (t_cur > time_duration_default - 2e-3)
          hovering_pub_count = 0;
      }
    }
    else
    {
      if (has_default)
      {
        ros::Time time_now = ros::Time::now();
        double t_cur = (time_now - time_hovering_start).toSec();
        start_pt = default_traj.evaluateDeBoorT(t_cur);
        start_vel= default_vel.evaluateDeBoorT(t_cur);
        start_acc = default_acc.evaluateDeBoor(t_cur);
        has_default = false;
      }
      changeState(GEN_NEW_TRAJ, "STATE");
    }
    break;
  }

  case GEN_NEW_TRAJ: {
    bool success = trajGeneration();
    if (success)
      changeState(EXEC_TRAJ, "STATE");
    
    else
      changeState(GEN_NEW_TRAJ, "STATE");
    break;
  }

  case EXEC_TRAJ: {
    ros::Time time_now = ros::Time::now();
    double t_cur = (time_now - time_traj_start).toSec();
    double t_replan = ros::Duration(2, 0).toSec();
    t_cur = min(time_duration, t_cur);

    if (t_cur > time_duration - 1e-4) {
      has_target = false;
      changeState(WAIT_TARGET, "STATE");
      return;
    } else if ((target_pt - odom_pt).norm() < no_replan_thresh) {
      return;
    } else if ((start_pt - odom_pt).norm() < replan_thresh) {
      return;
    } else if (t_cur < t_replan) {
      return;
    } else {
      printf("replan due to replan time out\n");
      changeState(REPLAN_TRAJ, "STATE");
    }
    break;
  }
  case REPLAN_TRAJ: {
    ros::Time time_now = ros::Time::now();
    double t_cur = (time_now - time_traj_start).toSec();
    double t_delta = ros::Duration(0, 0).toSec();
    t_cur = t_delta + t_cur;
    //start_pt = getPos(t_cur);
    //start_vel = getVel(t_cur);
    start_pt = cur_bspline_traj.evaluateDeBoorT(t_cur);
    if (if_2d_search == 1) start_pt(2) = h_reference;
    start_vel = cur_bspline_vel.evaluateDeBoorT(t_cur);
    start_acc = cur_bspline_acc.evaluateDeBoorT(t_cur);
    //std::cout<< "current time: " << t_cur << std::endl << start_pt.transpose() << std::endl;
    //start_vel = odom_vel;
    bool success = trajGeneration();
    if (success)
      changeState(EXEC_TRAJ, "STATE");
    else
      changeState(GEN_NEW_TRAJ, "STATE");
    break;
  }
  }
}

void rcvWaypointsCallBack(const nav_msgs::Path &wp) {
  if (wp.poses[0].pose.position.z < 0.0)
    return;
  //target_pt << wp.poses[0].pose.position.x, wp.poses[0].pose.position.y,
  //    wp.poses[0].pose.position.z;
  target_pt << 255.0, 0.0, 0.000;
  if (if_2d_search == 1) target_pt(2) = h_reference;
  ROS_INFO("[node] receive the planning target");
  start_pt = odom_pt;
  if (if_2d_search == 1) start_pt(2) = h_reference;
  start_vel = odom_vel;
  has_target = true;

  if (exec_state == WAIT_TARGET)
  {
    changeState(GEN_NEW_TRAJ, "STATE");
    if (has_default)
    {
      ros::Time time_now = ros::Time::now();
      double t_cur = (time_now - time_hovering_start).toSec();
      start_pt = default_traj.evaluateDeBoorT(t_cur);
      start_vel= default_vel.evaluateDeBoorT(t_cur);
      start_acc = default_acc.evaluateDeBoorT(t_cur);
      has_default = false;
    }
  }
  else if (exec_state == EXEC_TRAJ){
    //printf("replan due to rcv waypoint\n");
    changeState(REPLAN_TRAJ, "STATE");}
}

void bspline_traj_pub(UniformBspline* uniform_bspline_ptr)
{
  quadrotor_msgs::Bspline bspline;
  bspline.order = 3;
  bspline.start_time = ros::Time::now();
  bspline.traj_id = 1;
  Eigen::MatrixXd pos_pts = uniform_bspline_ptr->getControlPoint();
  
  bspline.pos_pts.reserve(pos_pts.cols());
  printf("Trajectory_generator: publishing bspline traj\n");
  //std::cout << pos_pts << std::endl;
  for (int i = 0; i < pos_pts.cols(); ++i)
  {
    geometry_msgs::Point pt;
    pt.x = pos_pts(0, i);
    pt.y = pos_pts(1, i);
    pt.z = pos_pts(2, i);
    bspline.pos_pts.push_back(pt);
  }
  Eigen::VectorXd knots = uniform_bspline_ptr->getKnot();
  if (knots.rows() <= 7)
    return;
  bspline.knots.reserve(knots.rows());
  for (int i = 0; i < knots.rows(); ++i)
  {
    bspline.knots.push_back(knots(i));
  }

  bspline_pub.publish(bspline);
}

void rcvPointCloudCallBack(const sensor_msgs::PointCloud2 &pointcloud_map) {

  pcl::PointCloud<pcl::PointXYZ> cloud;
  pcl::PointCloud<pcl::PointXYZ> cloud_vis;
  sensor_msgs::PointCloud2 map_vis;

  pcl::fromROSMsg(pointcloud_map, cloud);

  if ((int)cloud.points.size() == 0)
    return;

  pcl::PointXYZ pt;
  //printf("receive point cloud %d\n", (int)cloud.points.size());
  for (int idx = 0; idx < (int)cloud.points.size(); idx++) {
    pt = cloud.points[idx];
    // set obstalces into grid map for path planning
    _astar_path_finder->setObs(pt.x, pt.y, pt.z);
  }
}

// trajectory generation: front-end + back-end
// front-end : A* search method
// back-end  : Minimum snap trajectory generation
bool trajGeneration() {
  /**
   *
   * STEP 1:  search A* path and initialize bspline traj 
   *
   * **/
  ros::Time t1, t2;
  _astar_path_finder->resetUsedGrids();
  t1 = ros::Time::now();
  _astar_path_finder->AstarGraphSearch(start_pt, target_pt, if_2d_search);
  auto grid_path = _astar_path_finder->getPath();
  grid_path[0] = target_pt;
  grid_path[grid_path.size()-1] = start_pt;
  BsplineOpt bspline_opt;
  bspline_opt.set_param(nh_ptr, _astar_path_finder);

  vector<Vector3d> start_target_derivative;
  start_target_derivative.push_back(start_vel);
  if (has_default) start_target_derivative.push_back(target_vel);
  else 
  {
    Vector3d vel_direction(0, 0, 0);
    vel_direction(0) = (grid_path[0]-grid_path[1])(0);
    vel_direction(1) = (grid_path[0]-grid_path[1])(1);
    start_target_derivative.push_back(vel_direction.normalized()*hovering_speed);
    target_vel = start_target_derivative[1];
  }
  start_target_derivative.push_back(start_acc);
  if (has_default) start_target_derivative.push_back(target_acc);
  else
  {
    start_target_derivative.push_back(hovering_speed*hovering_speed/hovering_radius*zRotate(start_target_derivative[1].normalized(), -dPI/2));
    target_acc = start_target_derivative[3];
  }

  bool success_flag = bspline_opt.set_bspline(grid_path, start_target_derivative);
  cur_bspline_traj = bspline_opt.get_bspline();
  visTrajectory(false, VIS_INIT);
  /**
   *
   * STEP 2:  optimize traj
   *
   * **/
  t2 = ros::Time::now();

  double init_time = (t2-t1).toSec();
  t1 = t2;
  if (!success_flag) return false;
  bspline_opt.optStage();
  t2 = ros::Time::now();
  cur_bspline_traj = bspline_opt.get_bspline();
  visTrajectory(false, VIS_OPT);
  double opt_time = (t2-t1).toSec();
  
  /**
   *
   * STEP 3:  Trajectory optimization
   *
   * **/
  t1 = t2;
  bspline_opt.adjStage();
  t2 = ros::Time::now();

  double adj_time = (t2-t1).toSec();
  printf("Init Opt Adj Time: %f, %f, %f\n", init_time, opt_time, adj_time);
  UniformBspline bspline(bspline_opt.get_bspline());
  //std::cout<<bspline.get_control_points()<<std::endl;
  time_traj_start = ros::Time::now();
  time_duration = bspline.getTimeSum();
  bspline_traj_pub(&bspline);
  cur_bspline_traj = bspline;
  cur_bspline_vel = cur_bspline_traj.getDerivative();
  cur_bspline_acc = cur_bspline_vel.getDerivative();
  visTrajectory(false, VIS_TRAJ);
  if (cur_bspline_traj.get_control_points().cols() > 2)
  {
    defaultTrajGeneration();
    time_duration_default = default_traj.getTimeSum();
    default_vel = default_traj.getDerivative();
    default_acc = default_vel.getDerivative();
    has_default = true;
  }

  return true;
}

/*
  generate cyclical waypoints for hovering above goal
*/
bool defaultTrajGeneration()
{
  MatrixXd cur_control_pt = cur_bspline_traj.getControlPoint();
  if (cur_control_pt.cols()<3)
    return false;
  Vector3d p0;
  p0 = -target_vel.normalized();
  Vector3d hovering_center;
  hovering_center = target_pt + hovering_radius*zRotate(p0, dPI/2);
  vector<Vector3d> points;
  double rotate_angles[7] = //{-2*dPI/3, dPI, 2*dPI/3, dPI/3, 0, -dPI/3, -2*dPI/3};
  {-5*dPI/6, 5*dPI/6, dPI/2, dPI/6, -dPI/6};
  points.push_back(target_pt);
  for (int i=0;i<5;i++)
  {
    points.push_back(hovering_center+hovering_radius * zRotate(p0, rotate_angles[i]));
  }
  points.push_back(target_pt);
  vector<Vector3d> start_target_derivative;
  start_target_derivative.push_back(-p0*hovering_speed);
  start_target_derivative.push_back(-p0*hovering_speed);
  start_target_derivative.push_back(hovering_speed*hovering_speed/hovering_radius*zRotate(p0, dPI/2));
  start_target_derivative.push_back(start_target_derivative[2]);

  MatrixXd control_pt;
  default_traj.parameterizeToBspline(hovering_radius*2*dPI/6/hovering_speed, points, start_target_derivative, control_pt);
  default_traj.setUniformBspline(control_pt, 3, 0.8*hovering_radius/hovering_speed);
  return true;
}

Vector3d zRotate(Vector3d origin_vec, double angle)
{
  Matrix3d rot_mat;
  rot_mat <<  cos(angle), sin(angle), 0,
              -sin(angle),cos(angle), 0,
              0         , 0         , 1;
  return rot_mat*origin_vec;
}

void trajOptimization(Eigen::MatrixXd path) {
  // if( !has_odom ) return;
  MatrixXd vel = MatrixXd::Zero(2, 3);
  MatrixXd acc = MatrixXd::Zero(2, 3);

  vel.row(0) = start_vel;

  /**
   *
   * STEP 3.1:  finish the timeAllocation() using resonable allocation
   *
   * **/
  _polyTime = timeAllocation(path);

  /**
   *
   * STEP 3.2:  generate a minimum-jerk piecewise monomial polynomial-based
   * trajectory
   *
   * **/
  _polyCoeff =
      _trajGene->PolyQPGeneration(_dev_order, path, vel, acc, _polyTime);

  // check if the trajectory is safe, if not, do reoptimize
  int unsafe_segment;

  /**
   *
   * STEP 3.3:  finish the safeCheck()
   *
   * **/
  unsafe_segment = _astar_path_finder->safeCheck(_polyCoeff, _polyTime);

  MatrixXd repath = path;
  int count = 0;
  while (unsafe_segment != -1) {
    /**
     *
     * STEP 3.4:  reoptimize
     * the method comes from: "Polynomial Trajectory
     * Planning for Aggressive Quadrotor Flight in Dense Indoor Environment"
     * part 3.5
     *
     * **/
  }
  // visulize path and trajectory
  visPath(repath);
  //visTrajectory(_polyCoeff, _polyTime);
}

void trajPublish(MatrixXd polyCoeff, VectorXd time) {
  if (polyCoeff.size() == 0 || time.size() == 0) {
    ROS_WARN("[trajectory_generator_waypoint] empty trajectory, nothing to "
             "publish.");
    return;
  }

  unsigned int poly_number;

  static int count =
      1; // The first trajectory_id must be greater than 0. zxzxzxzx

  quadrotor_msgs::PolynomialTrajectory traj_msg;

  traj_msg.header.seq = count;
  traj_msg.header.stamp = ros::Time::now();
  traj_msg.header.frame_id = std::string("world");
  traj_msg.trajectory_id = count;
  traj_msg.action = quadrotor_msgs::PolynomialTrajectory::ACTION_ADD;

  traj_msg.num_order = 2 * _dev_order - 1; // the order of polynomial
  traj_msg.num_segment = time.size();

  Vector3d initialVel, finalVel;
  initialVel = _trajGene->getVelPoly(_polyCoeff, 0, 0);
  finalVel = _trajGene->getVelPoly(_polyCoeff, traj_msg.num_segment - 1,
                                   _polyTime(traj_msg.num_segment - 1));
  traj_msg.start_yaw = atan2(initialVel(1), initialVel(0));
  traj_msg.final_yaw = atan2(finalVel(1), finalVel(0));

  poly_number = traj_msg.num_order + 1;

  for (unsigned int i = 0; i < traj_msg.num_segment; i++) {
    for (unsigned int j = 0; j < poly_number; j++) {
      traj_msg.coef_x.push_back(polyCoeff(i, j) * pow(time(i), j));
      traj_msg.coef_y.push_back(polyCoeff(i, poly_number + j) *
                                pow(time(i), j));
      traj_msg.coef_z.push_back(polyCoeff(i, 2 * poly_number + j) *
                                pow(time(i), j));
    }
    traj_msg.time.push_back(time(i));
    traj_msg.order.push_back(traj_msg.num_order);
  }  // cout << "p_order:" << poly_number << endl;
  // cout << "traj_msg.num_order:" << traj_msg.num_order << endl;
  // cout << "traj_msg.num_segment:" << traj_msg.num_segment << endl;
  traj_msg.mag_coeff = 1;

  count++;
  ROS_WARN("[traj..gen...node] traj_msg publish");
  _traj_pub.publish(traj_msg);
}

VectorXd timeAllocation(MatrixXd Path) {
  VectorXd time(Path.rows() - 1);

  return time;
}

void visTrajectory(bool if_default, VIS_TYPE vis_type) {
  visualization_msgs::Marker _traj_vis;

  _traj_vis.header.stamp = ros::Time::now();
  _traj_vis.header.frame_id = "world";

  _traj_vis.ns = "traj_node/trajectory";
  //_traj_vis.id = vis_traj_id;
  _traj_vis.id = 0;
  vis_traj_id++;
  _traj_vis.type = visualization_msgs::Marker::LINE_STRIP;
  _traj_vis.action = visualization_msgs::Marker::ADD;
  
  _traj_vis.pose.orientation.x = 0.0;
  _traj_vis.pose.orientation.y = 0.0;
  _traj_vis.pose.orientation.z = 0.0;
  _traj_vis.pose.orientation.w = 1.0;

  if (vis_type == VIS_TRAJ)
  {
    _traj_vis.scale.x = _vis_traj_width;
    _traj_vis.scale.y = _vis_traj_width;
    _traj_vis.scale.z = _vis_traj_width;
    
    _traj_vis.color.a = 0.9;
    _traj_vis.color.r = 0.5373;
    _traj_vis.color.g = 0.5137;
    _traj_vis.color.b = 0.7490;
  } else if (vis_type == VIS_OPT)
  {
    _traj_vis.scale.x = _vis_traj_width*0.5;
    _traj_vis.scale.y = _vis_traj_width*0.5;
    _traj_vis.scale.z = _vis_traj_width*0.5;
    
    _traj_vis.color.a = 1.0;
    _traj_vis.color.r = 0.0;
    _traj_vis.color.g = 0.5;
    _traj_vis.color.b = 1.0;
  } else if (vis_type == VIS_INIT)
  {
    _traj_vis.scale.x = _vis_traj_width*0.5;
    _traj_vis.scale.y = _vis_traj_width*0.5;
    _traj_vis.scale.z = _vis_traj_width*0.5;
    _traj_vis.color.a = 0.9;
    _traj_vis.color.r = 0.3294;
    _traj_vis.color.g = 0.7020;
    _traj_vis.color.b = 0.2706;
  }

  //_traj_vis.lifetime = ros::Duration(0.0);

  _traj_vis.points.clear();
  Vector3d pos;
  geometry_msgs::Point pt;

  double tm, tmp;
  if (if_default)
    default_traj.getTimeSpan(tm, tmp);
  else
    cur_bspline_traj.getTimeSpan(tm, tmp);

  for (double t = tm; t <= tmp; t += 0.1) {
    if (if_default)
      pos = default_traj.evaluateDeBoorT(t);
    else
      pos = cur_bspline_traj.evaluateDeBoorT(t);
    //cout << pos << endl;
    pt.x = pos(0);
    pt.y = pos(1);
    pt.z = pos(2);
    _traj_vis.points.push_back(pt);
  }
  if (vis_type == VIS_TRAJ)
    _traj_vis_pub.publish(_traj_vis);
  else if (vis_type == VIS_OPT)
    _traj_opt_vis_pub.publish(_traj_vis);
  else if (vis_type == VIS_INIT)
    _traj_init_vis_pub.publish(_traj_vis);
}

void visPath(MatrixXd nodes) {
  visualization_msgs::Marker points;

  int id = 0;
  points.id = id;
  points.type = visualization_msgs::Marker::SPHERE_LIST;
  points.header.frame_id = "world";
  points.header.stamp = ros::Time::now();
  points.ns = "traj_node/path";
  points.action = visualization_msgs::Marker::ADD;
  points.pose.orientation.w = 1.0;
  points.pose.orientation.x = 0.0;
  points.pose.orientation.y = 0.0;
  points.pose.orientation.z = 0.0;
  points.scale.x = 0.2;
  points.scale.y = 0.2;
  points.scale.z = 0.2;
  points.color.a = 1.0;
  points.color.r = 0.0;
  points.color.g = 0.0;
  points.color.b = 1.0;

  geometry_msgs::Point p;
  for (int i = 0; i < int(nodes.rows()); i++) {
    p.x = nodes(i, 0);
    p.y = nodes(i, 1);
    p.z = nodes(i, 2);

    points.points.push_back(p);
  }
  _path_vis_pub.publish(points);
}

Vector3d getPos(double t_cur) {
  double time = 0;
  Vector3d pos = Vector3d::Zero();
  for (int i = 0; i < _polyTime.size(); i++) {
    for (double t = 0.0; t < _polyTime(i); t += 0.01) {
      time = time + 0.01;
      if (time > t_cur) {
        pos = _trajGene->getPosPoly(_polyCoeff, i, t);
        return pos;
      }
    }
  }
  return pos;
}

Vector3d getVel(double t_cur) {
  double time = 0;
  Vector3d Vel = Vector3d::Zero();
  for (int i = 0; i < _polyTime.size(); i++) {
    for (double t = 0.0; t < _polyTime(i); t += 0.01) {
      time = time + 0.01;
      if (time > t_cur) {
        Vel = _trajGene->getVelPoly(_polyCoeff, i, t);
        return Vel;
      }
    }
  }
  return Vel;
}

int main(int argc, char **argv) {
  ros::init(argc, argv, "traj_node");
  ros::NodeHandle nh("~");
  nh_ptr = &nh;

  nh.param("planning/vel", _Vel, 1.0);
  nh.param("planning/acc", _Acc, 1.0);
  nh.param("planning/dev_order", _dev_order, 3);
  nh.param("planning/min_order", _min_order, 3);
  nh.param("vis/vis_traj_width", _vis_traj_width, 0.15);
  nh.param("map/resolution", _resolution, 0.2);
  nh.param("map/x_size", _x_size, 50.0);
  nh.param("map/y_size", _y_size, 50.0);
  nh.param("map/z_size", _z_size, 5.0);
  nh.param("path/resolution", _path_resolution, 0.05);
  nh.param("replanning/thresh_replan", replan_thresh, -1.0);
  nh.param("replanning/thresh_no_replan", no_replan_thresh, -1.0);
  nh.param("hovering_speed", hovering_speed, 3.0);
  nh.param("hovering_radius", hovering_radius, 3.0);
  nh.param("planning/if_2d_search", if_2d_search, 1);
  nh.param("planning/h_reference", h_reference, 15.0);
  num_cp = 40;

  _poly_num1D = 2 * _dev_order;

  _exec_timer = nh.createTimer(ros::Duration(0.01), execCallback);

  _odom_sub = nh.subscribe("odom", 10, rcvOdomCallback);
  _map_sub = nh.subscribe("/pcl_render_node/local_pointcloud", 1, rcvPointCloudCallBack);
  _pts_sub = nh.subscribe("waypoints", 1, rcvWaypointsCallBack);

  _traj_pub =
      nh.advertise<quadrotor_msgs::PolynomialTrajectory>("trajectory", 50);
  _default_traj_pub = nh.advertise<visualization_msgs::Marker>("default_path", 1);
  _traj_init_vis_pub = nh.advertise<visualization_msgs::Marker>("vis_init_trajectory", 1);
  _traj_opt_vis_pub = nh.advertise<visualization_msgs::Marker>("vis_opt_trajectory", 1);
  _traj_vis_pub = nh.advertise<visualization_msgs::Marker>("vis_trajectory", 1);
  _path_vis_pub = nh.advertise<visualization_msgs::Marker>("vis_path", 1);
  bspline_pub = nh.advertise<quadrotor_msgs::Bspline>("bspline_trajectory", 10);
  // set the obstacle map
  _map_lower << -_x_size / 2.0, -_y_size / 2.0, 0.0;
  _map_upper << +_x_size / 2.0, +_y_size / 2.0, _z_size;
  _inv_resolution = 1.0 / _resolution;
  _max_x_id = (int)(_x_size * _inv_resolution);
  _max_y_id = (int)(_y_size * _inv_resolution);
  _max_z_id = (int)(_z_size * _inv_resolution);

  _astar_path_finder = new AstarPathFinder();
  //TrajectoryServer server(nh);
  _astar_path_finder->initGridMap(_resolution, _map_lower, _map_upper,
                                  _max_x_id, _max_y_id, _max_z_id);

  ros::Rate rate(100);
  bool status = ros::ok();
  while (status) {
    ros::spinOnce();
    status = ros::ok();
    rate.sleep();
  }

  return 0;
}