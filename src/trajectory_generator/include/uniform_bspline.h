#ifndef _UNIFORM_BSPLINE_H_
#define _UNIFORM_BSPLINE_H_

#include <Eigen/Eigen>
#include <algorithm>
#include <iostream>
#include <ros/ros.h>
#include "Astar_searcher.h"

using namespace std;

// An implementation of non-uniform B-spline with different dimensions
// It also represents uniform B-spline which is a special case of non-uniform
class UniformBspline
{
private:
  // control points for B-spline with different dimensions.
  // Each row represents one single control point
  // The dimension is determined by column number
  // e.g. B-spline with N points in 3D space -> Nx3 matrix
  Eigen::MatrixXd control_points_;

  int p_, n_, m_;     // p degree, n+1 control points, m = n+p+1
  bool flag;
  Eigen::VectorXd u_; // knots vector
  double interval_;   // knot span \delta t

  Eigen::MatrixXd getDerivativeControlPoints();

  double limit_vel_, limit_acc_, limit_ratio_, feasibility_tolerance_; // physical limits and time adjustment ratio

public:
  UniformBspline() {flag = false;}
  UniformBspline(const Eigen::MatrixXd &points, const int &order, const double &interval);
  ~UniformBspline();

  Eigen::MatrixXd get_control_points(void) { return control_points_; }

  // initialize as an uniform B-spline
  void setUniformBspline(const Eigen::MatrixXd &points, const int &order, const double &interval);

  // get / set basic bspline info

  void setKnot(const Eigen::VectorXd &knot);
  Eigen::VectorXd getKnot();
  bool getFlag() {return flag;}
  Eigen::MatrixXd getControlPoint();
  double getInterval();
  bool getTimeSpan(double &um, double &um_p);

  // compute position / derivative

  Eigen::VectorXd evaluateDeBoor(const double &u);                                               // use u \in [up, u_mp]
  inline Eigen::VectorXd evaluateDeBoorT(const double &t) { return evaluateDeBoor(t + u_(p_)); } // use t \in [0, duration]
  inline int getControlPtSize(){return control_points_.cols();}
  UniformBspline getDerivative();

  // 3D B-spline interpolation of points in point_set, with boundary vel&acc
  // constraints
  // input : (K+2) points with boundary vel/acc; ts
  // output: (K+6) control_pts
  static void parameterizeToBspline(const double &ts, const vector<Eigen::Vector3d> &point_set,
                                    const vector<Eigen::Vector3d> &start_end_derivative,
                                    Eigen::MatrixXd &ctrl_pts);

  /* check feasibility, adjust time */

  void setPhysicalLimits(const double &vel, const double &acc, const double &tolerance);
  bool checkFeasibility(double &ratio, bool show = false);
  void lengthenTime(const double &ratio);

  /* for performance evaluation */

  double getTimeSum();
  double getLength(const double &res = 0.01);
  double getJerk();
  void getMeanAndMaxVel(double &mean_v, double &max_v);
  void getMeanAndMaxAcc(double &mean_a, double &max_a);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class BsplineOpt
{
private:
    enum FORCE_STOP_OPTIMIZE_TYPE
    {
      DONT_STOP,
      STOP_FOR_REBOUND,
      STOP_FOR_ERROR
    } force_stop_type_;
    int iter_num_;
    int variable_num_;
    double min_cost_;

    double lambda1_;              // smoothness weight
    double lambda2_;              // collision weight
    double lambda3_;              // feasibility weight
    double lambda4_;              // curvature weight
    double lambda5_;              // fitness weight
    double dist0_;                // safe distance
    double max_vel_, min_vel_, max_acc_;    // dynamic limits
    int order_;                   // spline order
    double cp_dist_;              // control points distance
    Eigen::Vector3d end_pt_;
    Eigen::Vector3d start_pt_;
    Eigen::MatrixXd control_pts;

    UniformBspline bspline;
    AstarPathFinder* path_finder;

    // q control points mat, double cost, mat gradient
    void calcSmoothnessCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient);
    void calcFeasibilityCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient);
    void calcCollisionCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient);
    void calcFitnessCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient);
    void calcCurvatureCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient);

public:
    BsplineOpt() {};
    ~BsplineOpt() {};
    void set_param(ros::NodeHandle* nh, AstarPathFinder* new_path_finder);
    void set_bspline(std::vector<Eigen::Vector3d> A_Star_Path, std::vector<Eigen::Vector3d> start_target_derivetive);
    UniformBspline get_bspline();
    void combineOptCost(const double *x, double *grad, double &f_combine, const int n);
    void combineAdjCost(const double *x, double *grad, double &f_combine, const int n);
    static int earlyExit(void *func_data, const double *x, const double *g, const double fx, const double xnorm, const double gnorm, const double step, int n, int k, int ls);
    static double costFunctionOpt(void *func_data, const double *x, double *grad, const int n);
    bool optStage();
};
#endif