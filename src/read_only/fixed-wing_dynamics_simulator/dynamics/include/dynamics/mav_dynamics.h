#ifndef MAV_DYNAMICS
#define MAV_DYNAMICS

#include <ros/ros.h>
#include <mav_msgs/State.h>
#include <mav_msgs/Wind.h>
#include <mav_msgs/ControlInputs.h>
#include <Eigen/Core>
#include <cmath>
#include "tools/rotations.h"

typedef Eigen::Matrix<double, 13, 1> StateVec;
namespace dyn
{

class Dynamics
{
public:
  enum
  {
    POS = 0,
    VEL = 3,
    ATT = 6,
    OMEGA = 10,
    F = 0,
    M = 3
  };
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Dynamics();
  ~Dynamics();
  void run();

private:
  bool if_started;
  ros::NodeHandle nh_;
  ros::NodeHandle nh_p_;

  //Subscribers
  ros::Subscriber wind_sub;
  ros::Subscriber inputs_sub;

  //Publishers
  ros::Publisher state_pub;

  //callbacks
  void windCallback(const mav_msgs::WindConstPtr &msg);
  void inputCallback(const mav_msgs::ControlInputsConstPtr &msg);

  //Other functions
  void propogateDynamics(bool if_start);
  void updateVelocityData(const Eigen::Vector3d& gust);
  void calcGammaAndChi();
  void calculateForcesAndMoments();
  void calculateLongitudinalForces(double de);
  void calculateLateralForces(double da, double dr);
  void calculateThrustForce(double dt);
  StateVec derivatives(const StateVec& x);
  void loadParams();

  double Ts_, Va_, alpha_, beta_, chi_, flight_path_;
  StateVec x_;
  Eigen::Vector3d wind_, windg_, wind_ss;
  Eigen::Matrix<double, 6, 1> forces_;

  const double PI{3.14159265};
  double t0;
  double tprev;

  //Parameters from the yaml file
  double mass, gamma, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8;
  double Jy, g_;
  double delta_r, delta_a, delta_t, delta_e;

  double M_, rho, S_wing, c, alpha0, e, AR, b;
  double CL_0, CL_alpha, CL_q, CL_de;
  double CD_p, CD_alpha, CD_q, CD_de, CD_0;
  double Cm_0, Cm_alpha, Cm_q, Cm_de;
  double CY_0, CY_beta, CY_p, CY_r, CY_da, CY_dr;
  double Cell_0, Cell_beta, Cell_p, Cell_r, Cell_da, Cell_dr;
  double Cn_0, Cn_beta, Cn_p, Cn_r, Cn_da, Cn_dr;
  double D_prop, C_Q0, C_Q1, KQ, R_motor, C_Q2, i0;
  double C_T2, C_T1, C_T0, V_max;

};
}

#endif
