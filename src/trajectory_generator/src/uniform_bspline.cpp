#include "uniform_bspline.h"
#include "lbfgs.hpp"
#include <ros/ros.h>

UniformBspline::UniformBspline(const Eigen::MatrixXd &points, const int &order,
                                const double &interval)
{
  setUniformBspline(points, order, interval);
}

UniformBspline::~UniformBspline() {}

void UniformBspline::setUniformBspline(const Eigen::MatrixXd &points, const int &order,
                                        const double &interval)
{
  flag = true;
  control_points_ = points;
  p_ = order;
  interval_ = interval;

  n_ = points.cols() - 1;
  m_ = n_ + p_ + 1;

  u_ = Eigen::VectorXd::Zero(m_ + 1);
  for (int i = 0; i <= m_; ++i)
  {

    if (i <= p_)
    {
      u_(i) = double(-p_ + i) * interval_;
    }
    else if (i > p_ && i <= m_ - p_)
    {
      u_(i) = u_(i - 1) + interval_;
    }
    else if (i > m_ - p_)
    {
      u_(i) = u_(i - 1) + interval_;
    }
  }
  //std::cout<< "u_: " << u_ <<std::endl;
}

void UniformBspline::setKnot(const Eigen::VectorXd &knot) { this->u_ = knot; }

Eigen::VectorXd UniformBspline::getKnot() { return this->u_; }

bool UniformBspline::getTimeSpan(double &um, double &um_p)
{
  if (p_ > u_.rows() || m_ - p_ > u_.rows())
    return false;

  um = u_(p_);
  um_p = u_(m_ - p_);
  printf("time span: %f\n", um_p-um);
  return true;
}

Eigen::MatrixXd UniformBspline::getControlPoint() { return control_points_; }

Eigen::VectorXd UniformBspline::evaluateDeBoor(const double &u)
{

  double ub = min(max(u_(p_), u), u_(m_ - p_));

  // determine which [ui,ui+1] lay in
  int k = p_;
  while (true)
  {
    if (u_(k + 1) >= ub)
      break;
    ++k;
  }

  /* deBoor's alg */
  vector<Eigen::VectorXd> d;
  for (int i = 0; i <= p_; ++i)
  {
    d.push_back(control_points_.col(k - p_ + i));
    // cout << d[i].transpose() << endl;
  }

  for (int r = 1; r <= p_; ++r)
  {
    for (int i = p_; i >= r; --i)
    {
      double alpha = (ub - u_[i + k - p_]) / (u_[i + 1 + k - r] - u_[i + k - p_]);
      // cout << "alpha: " << alpha << endl;
      d[i] = (1 - alpha) * d[i - 1] + alpha * d[i];
    }
  }

  return d[p_];
}

// Eigen::VectorXd UniformBspline::evaluateDeBoorT(const double& t) {
//   return evaluateDeBoor(t + u_(p_));
// }

Eigen::MatrixXd UniformBspline::getDerivativeControlPoints()
{
  // The derivative of a b-spline is also a b-spline, its order become p_-1
  // control point Qi = p_*(Pi+1-Pi)/(ui+p_+1-ui+1)
  Eigen::MatrixXd ctp(control_points_.rows(), control_points_.cols() - 1);
  for (int i = 0; i < ctp.cols(); ++i)
  {
    ctp.col(i) =
        p_ * (control_points_.col(i + 1) - control_points_.col(i)) / (u_(i + p_ + 1) - u_(i + 1));
  }
  return ctp;
}

UniformBspline UniformBspline::getDerivative()
{
  Eigen::MatrixXd ctp = getDerivativeControlPoints();
  UniformBspline derivative(ctp, p_ - 1, interval_);

  /* cut the first and last knot */
  Eigen::VectorXd knot(u_.rows() - 2);
  knot = u_.segment(1, u_.rows() - 2);
  derivative.setKnot(knot);

  return derivative;
}

double UniformBspline::getInterval() { return interval_; }

void UniformBspline::setPhysicalLimits(const double &vel, const double &acc, const double &tolerance)
{
  limit_vel_ = vel;
  limit_acc_ = acc;
  limit_ratio_ = 1.1;
  feasibility_tolerance_ = tolerance;
}

bool UniformBspline::checkFeasibility(double &ratio, bool show)
{
  bool fea = true;

  Eigen::MatrixXd P = control_points_;
  int dimension = control_points_.rows();

  /* check vel feasibility and insert points */
  double max_vel = -1.0;
  double enlarged_vel_lim = limit_vel_ * (1.0 + feasibility_tolerance_) + 1e-4;
  for (int i = 0; i < P.cols() - 1; ++i)
  {
    Eigen::VectorXd vel = p_ * (P.col(i + 1) - P.col(i)) / (u_(i + p_ + 1) - u_(i + 1));

    if (fabs(vel(0)) > enlarged_vel_lim || fabs(vel(1)) > enlarged_vel_lim ||
        fabs(vel(2)) > enlarged_vel_lim)
    {

      if (show)
        cout << "[Check]: Infeasible vel " << i << " :" << vel.transpose() << endl;
      fea = false;

      for (int j = 0; j < dimension; ++j)
      {
        max_vel = max(max_vel, fabs(vel(j)));
      }
    }
  }

  /* acc feasibility */
  double max_acc = -1.0;
  double enlarged_acc_lim = limit_acc_ * (1.0 + feasibility_tolerance_) + 1e-4;
  for (int i = 0; i < P.cols() - 2; ++i)
  {

    Eigen::VectorXd acc = p_ * (p_ - 1) *
                          ((P.col(i + 2) - P.col(i + 1)) / (u_(i + p_ + 2) - u_(i + 2)) -
                            (P.col(i + 1) - P.col(i)) / (u_(i + p_ + 1) - u_(i + 1))) /
                          (u_(i + p_ + 1) - u_(i + 2));

    if (fabs(acc(0)) > enlarged_acc_lim || fabs(acc(1)) > enlarged_acc_lim ||
        fabs(acc(2)) > enlarged_acc_lim)
    {

      if (show)
        cout << "[Check]: Infeasible acc " << i << " :" << acc.transpose() << endl;
      fea = false;

      for (int j = 0; j < dimension; ++j)
      {
        max_acc = max(max_acc, fabs(acc(j)));
      }
    }
  }

  ratio = max(max_vel / limit_vel_, sqrt(fabs(max_acc) / limit_acc_));

  return fea;
}

void UniformBspline::lengthenTime(const double &ratio)
{
  int num1 = 5;
  int num2 = getKnot().rows() - 1 - 5;

  double delta_t = (ratio - 1.0) * (u_(num2) - u_(num1));
  double t_inc = delta_t / double(num2 - num1);
  for (int i = num1 + 1; i <= num2; ++i)
    u_(i) += double(i - num1) * t_inc;
  for (int i = num2 + 1; i < u_.rows(); ++i)
    u_(i) += delta_t;
}

// void UniformBspline::recomputeInit() {}

void UniformBspline::parameterizeToBspline(const double &ts, const vector<Eigen::Vector3d> &point_set,
                                            const vector<Eigen::Vector3d> &start_end_derivative,
                                            Eigen::MatrixXd &ctrl_pts)
{
  if (ts <= 0)
  {
    cout << "[B-spline]:time step error." << endl;
    return;
  }

  if (point_set.size() <= 3)
  {
    cout << "[B-spline]:point set have only " << point_set.size() << " points." << endl;
    return;
  }

  if (start_end_derivative.size() != 4)
  {
    cout << "[B-spline]:derivatives error." << endl;
  }

  int K = point_set.size();

  // write A
  Eigen::Vector3d prow(3), vrow(3), arow(3);
  prow << 1, 4, 1;
  vrow << -1, 0, 1;
  arow << 1, -2, 1;

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(K + 4, K + 2);

  for (int i = 0; i < K; ++i)
    A.block(i, i, 1, 3) = (1 / 6.0) * prow.transpose();

  A.block(K, 0, 1, 3) = (1 / 2.0 / ts) * vrow.transpose();
  A.block(K + 1, K - 1, 1, 3) = (1 / 2.0 / ts) * vrow.transpose();

  A.block(K + 2, 0, 1, 3) = (1 / ts / ts) * arow.transpose();
  A.block(K + 3, K - 1, 1, 3) = (1 / ts / ts) * arow.transpose();

  //cout << "A" << endl << A << endl << endl;

  // write b
  Eigen::VectorXd bx(K + 4), by(K + 4), bz(K + 4);
  for (int i = 0; i < K; ++i)
  {
    bx(i) = point_set[i](0);
    by(i) = point_set[i](1);
    bz(i) = point_set[i](2);
  }

  for (int i = 0; i < 4; ++i)
  {
    bx(K + i) = start_end_derivative[i](0);
    by(K + i) = start_end_derivative[i](1);
    bz(K + i) = start_end_derivative[i](2);
  }

  // solve Ax = b
  Eigen::VectorXd px = A.colPivHouseholderQr().solve(bx);
  Eigen::VectorXd py = A.colPivHouseholderQr().solve(by);
  Eigen::VectorXd pz = A.colPivHouseholderQr().solve(bz);

  // convert to control pts
  ctrl_pts.resize(3, K + 2);
  ctrl_pts.row(0) = px.transpose();
  ctrl_pts.row(1) = py.transpose();
  ctrl_pts.row(2) = pz.transpose();

  // cout << "[B-spline]: parameterization ok." << endl;
}

double UniformBspline::getTimeSum()
{
  double tm, tmp;
  if (getTimeSpan(tm, tmp))
    return tmp - tm;
  else
    return -1.0;
}

double UniformBspline::getLength(const double &res)
{
  double length = 0.0;
  double dur = getTimeSum();
  Eigen::VectorXd p_l = evaluateDeBoorT(0.0), p_n;
  for (double t = res; t <= dur + 1e-4; t += res)
  {
    p_n = evaluateDeBoorT(t);
    length += (p_n - p_l).norm();
    p_l = p_n;
  }
  return length;
}

double UniformBspline::getJerk()
{
  UniformBspline jerk_traj = getDerivative().getDerivative().getDerivative();

  Eigen::VectorXd times = jerk_traj.getKnot();
  Eigen::MatrixXd ctrl_pts = jerk_traj.getControlPoint();
  int dimension = ctrl_pts.rows();

  double jerk = 0.0;
  for (int i = 0; i < ctrl_pts.cols(); ++i)
  {
    for (int j = 0; j < dimension; ++j)
    {
      jerk += (times(i + 1) - times(i)) * ctrl_pts(j, i) * ctrl_pts(j, i);
    }
  }

  return jerk;
}

void UniformBspline::getMeanAndMaxVel(double &mean_v, double &max_v)
{
  UniformBspline vel = getDerivative();
  double tm, tmp;
  vel.getTimeSpan(tm, tmp);

  double max_vel = -1.0, mean_vel = 0.0;
  int num = 0;
  for (double t = tm; t <= tmp; t += 0.01)
  {
    Eigen::VectorXd vxd = vel.evaluateDeBoor(t);
    double vn = vxd.norm();

    mean_vel += vn;
    ++num;
    if (vn > max_vel)
    {
      max_vel = vn;
    }
  }

  mean_vel = mean_vel / double(num);
  mean_v = mean_vel;
  max_v = max_vel;
}

void UniformBspline::getMeanAndMaxAcc(double &mean_a, double &max_a)
{
  UniformBspline acc = getDerivative().getDerivative();
  double tm, tmp;
  acc.getTimeSpan(tm, tmp);

  double max_acc = -1.0, mean_acc = 0.0;
  int num = 0;
  for (double t = tm; t <= tmp; t += 0.01)
  {
    Eigen::VectorXd axd = acc.evaluateDeBoor(t);
    double an = axd.norm();

    mean_acc += an;
    ++num;
    if (an > max_acc)
    {
      max_acc = an;
    }
  }

  mean_acc = mean_acc / double(num);
  mean_a = mean_acc;
  max_a = max_acc;
}

void BsplineOpt::set_param(ros::NodeHandle* nh, AstarPathFinder* new_path_finder)
{
  printf("setting bspline param\n");
  nh->param("optimization/lambda_smooth", lambda1_, -1.0);
  nh->param("optimization/lambda_collision", lambda2_, -1.0);
  nh->param("optimization/lambda_feasibility", lambda3_, -1.0);
  nh->param("optimization/lambda_curvature", lambda4_, -1.0);
  nh->param("optimization/lambda_fitness", lambda5_, -1.0);
  
  nh->param("optimization/dist0", dist0_, -1.0);
  nh->param("optimization/max_vel", max_vel_, -1.0);
  nh->param("optimization/max_acc", max_acc_, -1.0);
  nh->param("optimization/min_vel", min_vel, -1.0);
  nh->param("optimization/order_", order_, 3);
  nh->param("optimization/control_points_distance", cp_dist_, 5.0);
  path_finder = new_path_finder;
}

void BsplineOpt::set_bspline(std::vector<Eigen::Vector3d> A_Star_Path, std::vector<Eigen::Vector3d> start_target_derivative)
{
  //printf("setting bspline\n");
  start_pt_ = A_Star_Path.front();
  end_pt_ = A_Star_Path.back();
  double dist_count = 0;
  //Eigen::MatrixXd points = start_pt_;
  std::vector<Eigen::Vector3d> points;
  Eigen::Vector3d fore_point;
  for (int i=1; i<A_Star_Path.size(); i++)
  {
    dist_count += (A_Star_Path[i] - A_Star_Path[i-1]).norm();
    if (dist_count > cp_dist_)
    {
      fore_point = start_pt_ - (A_Star_Path[i]-start_pt_);
      break;
    }
  }
  points.push_back(fore_point);
  points.push_back(start_pt_);
  
  for (int i=1; i<A_Star_Path.size(); i++)
  {
    dist_count += (A_Star_Path[i] - A_Star_Path[i-1]).norm();
    if (dist_count > cp_dist_)
    {
      points.push_back(A_Star_Path[i]);
      dist_count = 0;
    }
  }
  if (points.back() != end_pt_)
  {
    points.push_back(end_pt_);
  }
  //printf("setting control points\n");
  Eigen::MatrixXd control_points(int(points.size()), 3);
  int count = 0;
  vector<Eigen::Vector3d> points_inv;
  for (int i=points.size()-1; i>=0; i--)
  {
    /*
    for (int j=0;j<order_;j++)
    {
      //control_points.row(i) << points[i];
      control_points(count, j) = double(points[i](j));
    }*/
    points_inv.push_back(points[i]);
    count ++;
  }
  //cout<<control_points.transpose() << endl;
  double ts = cp_dist_ / max_vel_ * 5;
  bspline.parameterizeToBspline(ts,points_inv, start_target_derivative, control_points);
  bspline = UniformBspline(control_points, order_, ts);
  printf("gen\n");
}

UniformBspline BsplineOpt::get_bspline()
{
  return this->bspline;
}

void BsplineOpt::calcSmoothnessCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient)
{
  cost = 0.0;
  Eigen::Vector3d acc, temp_acc;

  for (int i = 0; i < q.cols() - 2; i++)
  {
    /* evaluate acc */
    acc = q.col(i + 2) - 2 * q.col(i + 1) + q.col(i);
    cost += acc.squaredNorm();
    temp_acc = 2.0 * acc;
    /* acc gradient */
    gradient.col(i + 0) += temp_acc;
    gradient.col(i + 1) += -2.0 * temp_acc;
    gradient.col(i + 2) += temp_acc;
  }
}

void BsplineOpt::calcFeasibilityCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient)
{
  cost = 0.0;
  double ts, ts_inv2;
  ts = bspline.getInterval();
  ts_inv2 = 1 / ts / ts;

  /* velocity feasibility */
  for (int i = 0; i < q.cols() - 1; i++)
  {
    Eigen::Vector3d vi = (q.col(i + 1) - q.col(i)) / ts;
    for (int j = 0; j < 3; j++)
    {
      if (vi(j) > max_vel_)
      {
        cost += pow(vi(j) - max_vel_, 2) * ts_inv2; // multiply ts_inv3 to make vel and acc has similar magnitude

        gradient(j, i + 0) += -2 * (vi(j) - max_vel_) / ts * ts_inv2;
        gradient(j, i + 1) += 2 * (vi(j) - max_vel_) / ts * ts_inv2;
      }
      else if (vi(j) < -max_vel_)
      {
        cost += pow(vi(j) + max_vel_, 2) * ts_inv2;

        gradient(j, i + 0) += -2 * (vi(j) + max_vel_) / ts * ts_inv2;
        gradient(j, i + 1) += 2 * (vi(j) + max_vel_) / ts * ts_inv2;
      }
      else if (vi(j) < min_vel_ && vi(j) > -min_vel_)
      {
        if (vi(j) < min_vel_ && vi(j) >= 0)
        {
          cost += pow(vi(j) - min_vel_, 2) * ts_inv2;

          gradient(j, i + 0) += -2 * (min_vel_ - vi(j)) / ts * ts_inv2;
          gradient(j, i + 1) += 2 * (min_vel_ - vi(j)) / ts * ts_inv2;
        }
        else
        {
          cost += pow(vi(j) + min_vel_, 2) * ts_inv2;

          gradient(j, i + 0) += 2 * (min_vel_ + vi(j)) / ts * ts_inv2;
          gradient(j, i + 1) += -2 * (min_vel_ + vi(j)) / ts * ts_inv2;
        }
      }
    }
  }
  /* acceleration feasibility */
  for (int i = 0; i < q.cols() - 2; i++)
  {
    Eigen::Vector3d ai = (q.col(i + 2) - 2 * q.col(i + 1) + q.col(i)) * ts_inv2;

    //cout << "temp_a * ai=" ;
    for (int j = 0; j < 3; j++)
    {
      if (ai(j) > max_acc_)
      {
        // cout << "fuck ACC" << endl;
        // cout << ai(j) << endl;
        cost += pow(ai(j) - max_acc_, 2);

        gradient(j, i + 0) += 2 * (ai(j) - max_acc_) * ts_inv2;
        gradient(j, i + 1) += -4 * (ai(j) - max_acc_) * ts_inv2;
        gradient(j, i + 2) += 2 * (ai(j) - max_acc_) * ts_inv2;
      }
      else if (ai(j) < -max_acc_)
      {
        cost += pow(ai(j) + max_acc_, 2);

        gradient(j, i + 0) += 2 * (ai(j) + max_acc_) * ts_inv2;
        gradient(j, i + 1) += -4 * (ai(j) + max_acc_) * ts_inv2;
        gradient(j, i + 2) += 2 * (ai(j) + max_acc_) * ts_inv2;
      }
      else
      {
        /* code */
      }
    }
  }
}

void BsplineOpt::calcCollisionCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient)
{
  cost = 0.0;
  int end_idx = q.cols() - order_;
  double demarcation = dist0_/2;
  double a = 3 * demarcation, b = -3 * pow(demarcation, 2), c = pow(demarcation, 3);
  for (auto i = order_; i < end_idx; ++i)
  {
    Eigen::Vector3d pt = q.col(i);
    Eigen::Vector3d gradient(0, 0, 0);
    double edt_result;
    path_finder->getEDTValueGradient(pt, edt_result, gradient);
    double dist_err = dist0_ - edt_result;

    if (dist_err < 0){}
    else if (dist_err < demarcation)
    {
      cost += pow(dist_err, 3);
      gradient.col(i) += -3.0 * dist_err * dist_err * gradient;
    }
    else
    {
      cost += a * dist_err * dist_err + b * dist_err + c;
      gradient.col(i) += -(2.0 * a * dist_err + b) * gradient;
    }
  }
}

void BsplineOpt::calcFitnessCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient)
{

}

void BsplineOpt::calcCurvatureCost(const Eigen::MatrixXd &q, double &cost, Eigen::MatrixXd &gradient)
{

}

void BsplineOpt::combineOptCost()
{
  
}

void BsplineOpt::combineAdjCost()
{

}