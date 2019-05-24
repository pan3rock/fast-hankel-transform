#ifndef FASTHANKELTRANSFORM_FASTHANKELTRANSFORM_H_
#define FASTHANKELTRANSFORM_FASTHANKELTRANSFORM_H_

#include <Eigen/Dense>
#include <complex>

class FastHankelTransform {
public:
  FastHankelTransform(int num_sample);
  ~FastHankelTransform();
  Eigen::VectorXd sampling();
  void set_feval(const Eigen::Ref<const Eigen::VectorXd> &feval);
  Eigen::VectorXd calculate();

  Eigen::VectorXd get_phi();

private:
  double evaluate_alpha();
  double evaluate_k0(double alpha);
  void evaluate_phi();
  void evaluate_j1();
  const int num_sample_;
  double alpha_;
  double k0_;

  Eigen::VectorXd x_;
  Eigen::VectorXd f_;
  std::complex<double> *phi_;
  std::complex<double> *j1_;

  bool x_updated_ = false;
  bool f_updated_ = false;

  const double tol_ = 1.0e-10;
};

#endif