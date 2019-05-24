#include "fast_hankel_transform.hpp"

#include <Eigen/Dense>
#include <boost/math/special_functions/bessel.hpp>
#include <cmath>
#include <complex>
#include <fftw3.h>

const double PI = 3.14159265358979323846;

using namespace Eigen;
using std::abs;
using std::exp;
using std::log;
using std::pow;

FastHankelTransform::FastHankelTransform(int num_sample)
    : num_sample_(num_sample), x_(VectorXd::Zero(num_sample_)),
      f_(VectorXd::Zero(num_sample_ + 1)),
      // phi_(VectorXd::Zero(num_sample_ * 2)),
      phi_(new std::complex<double>[num_sample_ * 2]),
      j1_(new std::complex<double>[num_sample_ * 2]) {
  alpha_ = evaluate_alpha();
  k0_ = evaluate_k0(alpha_);
}

FastHankelTransform::~FastHankelTransform() {
  delete[] phi_;
  delete[] j1_;
}

double FastHankelTransform::evaluate_alpha() {
  auto func = [&](auto a) { return -log(1.0 - exp(a)) / (num_sample_ - 1); };

  double alpha = 1.0;
  while (true) {
    double alpha_new = func(alpha);
    if (abs(alpha_new - alpha) < tol_) {
      alpha = alpha_new;
      break;
    } else {
      alpha = alpha_new;
    }
  }
  return alpha;
}

double FastHankelTransform::evaluate_k0(double alpha) {
  double k0 = (2.0 * exp(alpha) + exp(2.0 * alpha)) /
              (pow(1 + exp(alpha), 2) * (1 - exp(-2.0 * alpha)));
  return k0;
}

VectorXd FastHankelTransform::sampling() {
  x_(0) = (1.0 + exp(alpha_)) * exp(-alpha_ * num_sample_) / 2.0;
  for (auto i = 1; i < num_sample_; ++i) {
    x_(i) = x_(0) * exp(alpha_ * i);
  }
  return x_;
}

void FastHankelTransform::get_feval(const Ref<const VectorXd> &feval) {
  f_.head(num_sample_) = feval;
}

void FastHankelTransform::evaluate_phi() {
  phi_[0] = k0_ * (f_(0) - f_(1)) * exp(alpha_ * (1 - num_sample_));
  for (auto i = 1; i < num_sample_; ++i) {
    phi_[i] = (f_(i) - f_(i + 1)) * exp(alpha_ * (i + 1 - num_sample_));
  }
}

void FastHankelTransform::evaluate_j1() {
  for (auto i = 0; i < num_sample_ * 2; ++i) {
    double x = 2.0 * PI * x_(0) * exp(alpha_ * (i + 1 - num_sample_));
    j1_[i] = boost::math::sph_bessel(0, x);
  }
}

VectorXd FastHankelTransform::calculate() {
  int nsample = num_sample_ * 2;
  fftw_complex *fft_phi =
      (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nsample);
  fftw_complex *phi = reinterpret_cast<fftw_complex *>(phi_);

  fftw_plan p1 =
      fftw_plan_dft_1d(nsample, phi, fft_phi, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_complex *fft_j1 =
      (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nsample);
  fftw_complex *j1 = reinterpret_cast<fftw_complex *>(j1_);

  fftw_plan p2 =
      fftw_plan_dft_1d(nsample, j1, fft_j1, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftw_complex *in =
      (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nsample);

  for (auto i = 0; i < nsample; ++i) {
    in[i][0] = fft_phi[i][0] * fft_j1[i][0] - fft_phi[i][1] * fft_j1[i][1];
    in[i][1] = fft_phi[i][0] * fft_j1[i][1] + fft_phi[i][1] * fft_j1[i][0];
  }

  fftw_complex *out =
      (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nsample);

  fftw_plan p3 =
      fftw_plan_dft_1d(nsample, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  VectorXd ret(num_sample_);
  for (auto i = 0; i < num_sample_; ++i) {
    ret(i) = 1.0 / x_(i) * out[i][0];
  }

  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);
  fftw_destroy_plan(p3);
  fftw_free(fft_phi);
  fftw_free(fft_j1);
  fftw_free(in);
  fftw_free(out);

  return ret;
}