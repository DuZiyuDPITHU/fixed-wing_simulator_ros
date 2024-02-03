#include "uniform_bspline.h"
#include "nav_msgs/Odometry.h"
#include "quadrotor_msgs/Bspline.h"
#include "quadrotor_msgs/PositionCommand.h"
#include <mav_msgs/Commands.h>
#include "std_msgs/Empty.h"
#include "visualization_msgs/Marker.h"
#include <ros/ros.h>
#include <cmath>

ros::Publisher pos_cmd_pub;
ros::Publisher fw_pos_cmd_pub;
ros::Subscriber odom_sub;

quadrotor_msgs::PositionCommand cmd;
mav_msgs::Commands fw_cmd;
double pos_gain[3] = {0, 0, 0};
double vel_gain[3] = {0, 0, 0};


bool receive_traj_ = false;
vector<UniformBspline> traj_;
double traj_duration_;
ros::Time start_time_;
int traj_id_;

// yaw control
double last_yaw_, last_yaw_dot_;
double time_forward_;
Eigen::Vector3d odom_pt;
Eigen::Vector3d odom_vel;

/*
class PDcontrol
{
  double kp;
  double kd;
  double sigma = 0.01;
  double y_prev = 0.0;
  double eprev_ = 0.0;
  double y_dot = 0.0;
};

PDcontrol Va_controller, chi_controller;*/ 


void bsplineCallback(quadrotor_msgs::BsplineConstPtr msg)
{
  // parse pos traj
  Eigen::MatrixXd pos_pts(3, msg->pos_pts.size());

  Eigen::VectorXd knots(msg->knots.size());
  printf("points size: %d, knots size: %d\n", msg->pos_pts.size(), msg->knots.size());
  for (size_t i = 0; i < msg->knots.size(); ++i)
  {
    knots(i) = msg->knots[i];
  }
  
  for (size_t i = 0; i < msg->pos_pts.size(); ++i)
  {
    pos_pts(0, i) = msg->pos_pts[i].x;
    pos_pts(1, i) = msg->pos_pts[i].y;
    pos_pts(2, i) = msg->pos_pts[i].z;
  }

  UniformBspline pos_traj(pos_pts, msg->order, 0.1);
  pos_traj.setKnot(knots);

  // parse yaw traj

  // Eigen::MatrixXd yaw_pts(msg->yaw_pts.size(), 1);
  // for (int i = 0; i < msg->yaw_pts.size(); ++i) {
  //   yaw_pts(i, 0) = msg->yaw_pts[i];
  // }

  //UniformBspline yaw_traj(yaw_pts, msg->order, msg->yaw_dt);

  start_time_ = msg->start_time;
  traj_id_ = msg->traj_id;

  traj_.clear();
  traj_.push_back(pos_traj);
  traj_.push_back(traj_[0].getDerivative());
  traj_.push_back(traj_[1].getDerivative());

  traj_duration_ = traj_[0].getTimeSum();

  receive_traj_ = true;
  printf("received bspline traj for %f\n", traj_duration_);
}

double calculate_phi_ff(double t_cur)
{
  // clockwise rotation is declared as positive
  double phi_ff;
  double V_a = odom_vel.norm();
  Eigen::Vector3d dir = t_cur + time_forward_ <= traj_duration_ ? traj_[0].evaluateDeBoorT(t_cur + time_forward_) - odom_pt : traj_[0].evaluateDeBoorT(traj_duration_) - odom_pt;
  double L1 = dir.norm();
  double c_ita = dir.dot(odom_vel)/(V_a*L1);
  double s_ita = sqrt(1-c_ita);
  phi_ff = atan((2*V_a*V_a*s_ita)/(L1*9.8));

  Eigen::Vector3d L1_dir = dir.dot(odom_vel)>0? dir: -1*dir;
  if (odom_vel.cross(L1_dir)(2)>0)
    phi_ff *= -1;
  return phi_ff;
  /*
  double V_g = vel.norm();
  Eigen::Vector3d v_cross_a = vel.cross(acc);
  double kappa = v_cross_a.norm()/(V_g*V_g*V_g);
  double phi_ff_abs = abs(atan(kappa*V_g*V_g/9.8));
  if (v_cross_a(2)>0)
  {
    phi_ff_abs *= -1;
  }
  return phi_ff_abs;
  */
}

double calculate_Va(double t_cur)
{
  double tar_Va;
  Eigen::Vector3d dir = t_cur + time_forward_ <= traj_duration_ ? traj_[0].evaluateDeBoorT(t_cur + time_forward_) - odom_pt : traj_[0].evaluateDeBoorT(traj_duration_) - odom_pt;
  if (odom_vel.norm()==0) tar_Va = dir.dot(odom_vel.normalized())/time_forward_;
  else tar_Va = dir.dot(odom_vel.normalized())/time_forward_;
  if (tar_Va<0) tar_Va = 0;
  return tar_Va;
}

std::pair<double, double> calculate_yaw(double t_cur, Eigen::Vector3d &pos, ros::Time &time_now, ros::Time &time_last)
{
  constexpr double PI = 3.1415926;
  constexpr double YAW_DOT_MAX_PER_SEC = PI;
  // constexpr double YAW_DOT_DOT_MAX_PER_SEC = PI;
  std::pair<double, double> yaw_yawdot(0, 0);
  double yaw = 0;
  double yawdot = 0;

  Eigen::Vector3d dir = t_cur + time_forward_ <= traj_duration_ ? pos - traj_[0].evaluateDeBoorT(t_cur + time_forward_) : pos - traj_[0].evaluateDeBoorT(traj_duration_);
  //std::cout << "pos: " << pos << std::endl << "tar" << traj_[0].evaluateDeBoorT(t_cur + time_forward_) << std::endl;
  double yaw_temp = dir.norm() > 0.1 ? atan2(dir(1), dir(0)) : last_yaw_;
  
  double max_yaw_change = YAW_DOT_MAX_PER_SEC * (time_now - time_last).toSec();
  
  if (yaw_temp - last_yaw_ > PI)
  {
    if (yaw_temp - last_yaw_ - 2 * PI < -max_yaw_change)
    {
      yaw = last_yaw_ - max_yaw_change;
      if (yaw < -PI)
        yaw += 2 * PI;

      yawdot = -YAW_DOT_MAX_PER_SEC;
    }
    else
    {
      yaw = yaw_temp;
      if (yaw - last_yaw_ > PI)
        yawdot = -YAW_DOT_MAX_PER_SEC;
      else
        yawdot = (yaw_temp - last_yaw_) / (time_now - time_last).toSec();
    }
  }
  else if (yaw_temp - last_yaw_ < -PI)
  {
    if (yaw_temp - last_yaw_ + 2 * PI > max_yaw_change)
    {
      yaw = last_yaw_ + max_yaw_change;
      if (yaw > PI)
        yaw -= 2 * PI;

      yawdot = YAW_DOT_MAX_PER_SEC;
    }
    else
    {
      yaw = yaw_temp;
      if (yaw - last_yaw_ < -PI)
        yawdot = YAW_DOT_MAX_PER_SEC;
      else
        yawdot = (yaw_temp - last_yaw_) / (time_now - time_last).toSec();
    }
  }
  else
  {
    if (yaw_temp - last_yaw_ < -max_yaw_change)
    {
      yaw = last_yaw_ - max_yaw_change;
      if (yaw < -PI)
        yaw += 2 * PI;

      yawdot = -YAW_DOT_MAX_PER_SEC;
    }
    else if (yaw_temp - last_yaw_ > max_yaw_change)
    {
      yaw = last_yaw_ + max_yaw_change;
      if (yaw > PI)
        yaw -= 2 * PI;

      yawdot = YAW_DOT_MAX_PER_SEC;
    }
    else
    {
      yaw = yaw_temp;
      if (yaw - last_yaw_ > PI)
        yawdot = -YAW_DOT_MAX_PER_SEC;
      else if (yaw - last_yaw_ < -PI)
        yawdot = YAW_DOT_MAX_PER_SEC;
      else
        yawdot = (yaw_temp - last_yaw_) / (time_now - time_last).toSec();
    }
  }
  
  if (fabs(yaw - last_yaw_) <= max_yaw_change)
    yaw = 0.5 * last_yaw_ + 0.5 * yaw; // nieve LPF
  yawdot = 0.5 * last_yaw_dot_ + 0.5 * yawdot;
  last_yaw_ = yaw;
  last_yaw_dot_ = yawdot;

  yaw_yawdot.first = yaw;
  //yaw_yawdot.first = yaw_temp;
  yaw_yawdot.second = yawdot;

  return yaw_yawdot;
}

void cmdCallback(const ros::TimerEvent &e)
{
  /* no publishing before receive traj_ */
  if (!receive_traj_)
    return;

  ros::Time time_now = ros::Time::now();
  double t_cur = (time_now - start_time_).toSec();
  //printf("time now: %f, t_cur: %f, traj duration: %f\n", time_now.toSec(), t_cur, traj_duration_);

  Eigen::Vector3d pos(Eigen::Vector3d::Zero()), vel(Eigen::Vector3d::Zero()), acc(Eigen::Vector3d::Zero()), pos_f;
  std::pair<double, double> yaw_yawdot(0, 0);
  std::pair<double, double> fw_yaw_yawdot(0, 0);
  double phi_ff;
  double cmd_Va;

  static ros::Time time_last = ros::Time::now();
  if (t_cur < traj_duration_ && t_cur >= 0.0)
  {
    pos = traj_[0].evaluateDeBoorT(t_cur);
    vel = traj_[1].evaluateDeBoorT(t_cur);
    acc = traj_[2].evaluateDeBoorT(t_cur);

    /*** calculate yaw ***/
    //yaw_yawdot = calculate_yaw(t_cur, pos, time_now, time_last);
    fw_yaw_yawdot = calculate_yaw(t_cur, odom_pt, time_now, time_last);
    /*** calculate yaw ***/

    /*** calculate phi_ff and Va ***/
    phi_ff = calculate_phi_ff(t_cur);
    cmd_Va = calculate_Va(t_cur);
    /*** calculate phi_ff and Va ***/

    double tf = min(traj_duration_, t_cur + 2.0);
    pos_f = traj_[0].evaluateDeBoorT(tf);
  }
  else if (t_cur >= traj_duration_)
  {
    /* hover when finish traj_ */
    pos = traj_[0].evaluateDeBoorT(traj_duration_);
    vel.setZero();
    acc.setZero();

    //yaw_yawdot.first = last_yaw_;
    //yaw_yawdot.second = 0;

    pos_f = pos;
  }
  else
  {
    cout << "[Traj server]: invalid time." << endl;
  }
  time_last = time_now;

  cmd.header.stamp = time_now;
  cmd.header.frame_id = "world";
  cmd.trajectory_flag = quadrotor_msgs::PositionCommand::TRAJECTORY_STATUS_READY;
  cmd.trajectory_id = traj_id_;

  cmd.position.x = pos(0);
  cmd.position.y = pos(1);
  cmd.position.z = pos(2);
/*
  cmd.velocity.x = vel(0);
  cmd.velocity.y = vel(1);
  cmd.velocity.z = vel(2);

  cmd.acceleration.x = acc(0);
  cmd.acceleration.y = acc(1);
  cmd.acceleration.z = acc(2);

  cmd.yaw = yaw_yawdot.first;
  cmd.yaw_dot = yaw_yawdot.second;
  last_yaw_ = cmd.yaw;
*/
  //pos_cmd_pub.publish(cmd);

  fw_cmd.header.stamp = time_now;
  fw_cmd.header.frame_id = "world";

  fw_cmd.Va_cmd = cmd_Va;

  double fw_yaw_cmd = (fw_yaw_yawdot.first)/3.14159265*(180)+180;
  if (fw_yaw_cmd >=180) fw_yaw_cmd-=360;
  else if (fw_yaw_cmd <= -180) fw_yaw_cmd += 360;
  fw_cmd.chi_cmd = fw_yaw_cmd;
  fw_cmd.h_cmd = pos(2);
  fw_cmd.phi_ff = phi_ff;
  //printf("fw_yaw_cmd: %f\n", fw_yaw_cmd);
  fw_pos_cmd_pub.publish(fw_cmd);
//printf("chi: %f, phi_ff: %f, V_a: %f\n", yaw_yawdot.first/3.14159265*(180), phi_ff, vel.norm());
}

void odomCallback(const nav_msgs::Odometry::ConstPtr& odom)
{
  odom_pt(0) = odom->pose.pose.position.x;
  odom_pt(1) = odom->pose.pose.position.y;
  odom_pt(2) = odom->pose.pose.position.z;

  odom_vel(0) = odom->twist.twist.linear.x;
  odom_vel(1) = odom->twist.twist.linear.y;
  odom_vel(2) = odom->twist.twist.linear.z;
  //printf("V_a: %f\n", odom_vel.norm());
}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "traj_server");
  ros::NodeHandle node;
  ros::NodeHandle nh("~");
  printf("start traj server node\n");
  ros::Subscriber bspline_sub = node.subscribe("trajectory_generator_node/bspline_trajectory", 10, bsplineCallback);
  odom_sub = node.subscribe("/visual_slam/odom", 10, odomCallback);

  pos_cmd_pub = node.advertise<quadrotor_msgs::PositionCommand>("pos_cmd", 20);
  fw_pos_cmd_pub = node.advertise<mav_msgs::Commands>("fw_position_cmd", 20);

  ros::Timer cmd_timer = node.createTimer(ros::Duration(0.05), cmdCallback);

  /* control parameter */
  cmd.kx[0] = pos_gain[0];
  cmd.kx[1] = pos_gain[1];
  cmd.kx[2] = pos_gain[2];

  cmd.kv[0] = vel_gain[0];
  cmd.kv[1] = vel_gain[1];
  cmd.kv[2] = vel_gain[2];

  nh.param("traj_server/time_forward", time_forward_, 1.0);
  last_yaw_ = 0.0;
  last_yaw_dot_ = 0.0;

  ros::Duration(1.0).sleep();

  ROS_WARN("[Traj server]: ready.");

  ros::spin();

  return 0;
}