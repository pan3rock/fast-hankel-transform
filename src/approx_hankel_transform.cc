#include "approx_hankel_transform.hpp"

#include <Eigen/Dense>
#include <boost/math/special_functions/bessel.hpp>
#include <cmath>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <iostream>

using namespace Eigen;
using std::pow;

ApproxHankelTransform::ApproxHankelTransform(int num_sample)
    : num_sample_(num_sample), x_(VectorXd::Zero(num_sample_)),
      f_(VectorXd::Zero(num_sample_ * 2)) {}

VectorXd ApproxHankelTransform::sampling(double lb, double ub) {
  lb_ = lb;
  ub_ = ub;

  for (auto i = 0; i < num_sample_; ++i) {
    x_(i) = lb_ + (ub_ - lb_) * i / (num_sample_ - 1);
  }
  return x_;
}

void ApproxHankelTransform::set_feval(const Ref<const VectorXd> &feval) {
  f_.head(num_sample_) = feval;
}

double ApproxHankelTransform::calculate(double r) {
  MatrixXd mat(num_sample_ * 2, num_sample_ * 2);
  for (auto i = 0; i < num_sample_; ++i) {
    for (auto j = 0; j < num_sample_; ++j) {
      mat(j, i) = prime_basis(x_(j), i);
    }

    for (auto j = 0; j < num_sample_; ++j) {
      mat(j + num_sample_, i) = -r * basis(x_(j), i);
    }
  }

  for (auto i = 0; i < num_sample_; ++i) {
    for (auto j = 0; j < num_sample_; ++j) {
      mat(j, i + num_sample_) = r * basis(x_(j), i);
    }

    for (auto j = 0; j < num_sample_; ++j) {
      mat(j + num_sample_, i + num_sample_) =
          prime_basis(x_(j), i) - 1.0 / x_(j) * basis(x_(j), i);
    }
  }
  VectorXd coef = mat.fullPivLu().solve(f_);
  double ret = p1j0_p2j1(coef, r, ub_) - p1j0_p2j1(coef, r, lb_);
  return ret;
}

// \int_a^b [f(x)j0(rx) + g(x)j1(rx)]dx
double ApproxHankelTransform::p1j0_p2j1(const Ref<const VectorXd> &coef,
                                        double r, double x) {
  double p1 = 0.0;
  double p2 = 0.0;
  for (auto i = 0; i < num_sample_; ++i) {
    p1 += coef(i) * basis(x, i);
    p2 += coef(i + num_sample_) * basis(x, i);
  }
  double ret = p1 * boost::math::cyl_bessel_j(0, r * x) +
               p2 * boost::math::cyl_bessel_j(1, r * x);
  return ret;
}

double ApproxHankelTransform::basis(double x, int k) {
  const double d = (lb_ + ub_) / 2.0 + 1.0e-8;
  double ret = pow(x - d, k);
  return ret;
}

double ApproxHankelTransform::prime_basis(double x, int k) {
  const double d = (lb_ + ub_) / 2.0 + 1.0e-8;
  if (k == 0)
    return 0.0;
  double ret = k * pow(x - d, k - 1);
  return ret;
}