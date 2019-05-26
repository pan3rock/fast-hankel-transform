#ifndef FHT_APPROX_HANKEL_TRANSFORM_H_
#define FHT_APPROX_HANKEL_TRANSFORM_H_

#include <Eigen/Dense>

class ApproxHankelTransform {
public:
  ApproxHankelTransform(int num_sample);
  double calculate(double r);
  Eigen::VectorXd sampling(double lb, double ub);
  void set_feval(const Eigen::Ref<const Eigen::VectorXd> &feval);

private:
  double basis(double x, int k);
  double prime_basis(double x, int k);
  double p1j0_p2j1(const Eigen::Ref<const Eigen::VectorXd> &coef, double r,
                   double x);
  const int num_sample_;
  double lb_, ub_;
  Eigen::VectorXd x_, f_;
};

#endif
