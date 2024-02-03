#include "estimator/ekf.h"
#include <nav_msgs/Odometry.h>
#include <cmath>

namespace est
{
  EKF::EKF(): nh_(ros::NodeHandle()), nh_p("~")
  {
    state_sub = nh_.subscribe("true_states", 1, &EKF::stateCallback, this);
    estState_pub = nh_.advertise<mav_msgs::State>("estimated_states", 1);
    odom_pub = nh_.advertise<nav_msgs::Odometry>("/visual_slam/odom", 100);
    nh_.param("odom_rate", odom_rate, 20.0);
    last_odom_pub_time = ros::Time::now();
    //nh_.param("quadrotor_name", name, std::string("quadrotor"));
    name = "quadrotor";
  }

  EKF::~EKF(){}

  void EKF::stateCallback(const mav_msgs::StateConstPtr &msg)
  {
    mav_msgs::State estState;
    estState.header.stamp = ros::Time::now();
    estState.pn = msg->pn;
    estState.pe = msg->pe;
    estState.h = msg->h;
    estState.vn = msg->vn;
    estState.ve = msg->ve;
    estState.vh = msg->vh;
    estState.phi = msg->phi;
    estState.theta = msg->theta;
    estState.psi = msg->psi;
    estState.Va = msg->Va;
    estState.alpha = msg->alpha;
    estState.beta = msg->beta;
    estState.p = msg->p;
    estState.q = msg->q;
    estState.r = msg->r;
    estState.Vg = msg->Vg;
    estState.gamma = msg->gamma;
    estState.chi = msg->chi;
    estState.wn = msg->wn;
    estState.we = msg->we;
    estState.bx = 0.0;
    estState.by = 0.0;
    estState.bz = 0.0;

    estState_pub.publish(estState);

    ros::Time time_now = ros::Time::now();
    if ((time_now-last_odom_pub_time).toSec()>=(1/odom_rate))
    {
      nav_msgs::Odometry odom_msg;
      odom_msg.header.frame_id = "/simulator";
      odom_msg.child_frame_id  = "/" + name;
      odom_msg.pose.pose.position.y = msg->pe*1;
      odom_msg.pose.pose.position.x = msg->pn;
      odom_msg.pose.pose.position.z = msg->h*1;

      odom_msg.pose.pose.orientation.w = cos(msg->psi/2)*cos(msg->theta/2)*cos(msg->phi/2*(-1))+sin(msg->psi/2)*sin(msg->theta/2)*sin(msg->phi/2*(-1));
      odom_msg.pose.pose.orientation.x = cos(msg->psi/2)*cos(msg->theta/2)*sin(msg->phi/2*(-1))-sin(msg->psi/2)*sin(msg->theta/2)*cos(msg->phi/2*(-1));
      odom_msg.pose.pose.orientation.y = cos(msg->psi/2)*sin(msg->theta/2)*cos(msg->phi/2*(-1))+sin(msg->psi/2)*cos(msg->theta/2)*sin(msg->phi/2*(-1));
      odom_msg.pose.pose.orientation.z = sin(msg->psi/2)*cos(msg->theta/2)*cos(msg->phi/2*(-1))-cos(msg->psi/2)*sin(msg->theta/2)*sin(msg->phi/2*(-1));

      odom_msg.twist.twist.linear.y = msg->ve*1;
      odom_msg.twist.twist.linear.x = msg->vn;
      odom_msg.twist.twist.linear.z = msg->vh*1;

      odom_msg.twist.twist.angular.x = msg->p;
      odom_msg.twist.twist.angular.y = 1*msg->q;
      odom_msg.twist.twist.angular.z = msg->r;

      odom_pub.publish(odom_msg);
    }
  }
}
