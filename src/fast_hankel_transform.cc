/*
 * File: fast_hankel_transform.cc
 * Created Date: 2019-09-11
 * Author: Lei Pan
 * Contact: <panlei7@gmail.com>
 *
 * Last Modified: Wednesday September 25th 2019 11:37:38 am
 *
 * MIT License
 *
 * Copyright (c) 2019 Lei Pan
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the 'Software'), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * -----
 * HISTORY:
 * Date      	 By	Comments
 * ----------	---
 * ----------------------------------------------------------
 */

#include "fast_hankel_transform.hpp"

#include <Eigen/Dense>
#include <boost/math/special_functions/bessel.hpp>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <iostream>

using namespace Eigen;
using std::abs;
using std::exp;
using std::log;
using std::pow;

FastHankelTransform::FastHankelTransform(int num_sample, double ux, double uy)
    : num_sample_(num_sample), ux_(ux), uy_(uy),
      x_(VectorXd::Zero(num_sample_)), f_(VectorXd::Zero(num_sample_ + 1)),
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
  auto func = [&](auto a) { return -log(1.0 - exp(-a)) / (num_sample_ - 1); };

  const int maxiter = 100;
  double alpha = 1.0;
  for (auto i = 0; i < maxiter; ++i) {
    alpha = func(alpha);
  }
  return alpha;
}

double FastHankelTransform::evaluate_k0(double alpha) {
  double k0 = (2.0 * exp(alpha) + exp(2.0 * alpha)) /
              (pow(1 + exp(alpha), 2) * (1 - exp(-2.0 * alpha)));
  return k0;
}

VectorXd FastHankelTransform::sampling() {
  x_updated_ = true;
  x_(0) = (1.0 + exp(alpha_)) * exp(-alpha_ * num_sample_) / 2.0;
  for (auto i = 1; i < num_sample_; ++i) {
    x_(i) = x_(0) * exp(alpha_ * i);
  }
  return x_;
}

void FastHankelTransform::set_feval(const Ref<const VectorXd> &feval) {
  assert(x_updated_ == true);
  f_updated_ = true;
  f_.head(num_sample_) = feval;
}

void FastHankelTransform::evaluate_phi() {
  phi_[0] = k0_ * (f_(0) - f_(1)) * exp(alpha_ * (1 - num_sample_));
  for (auto i = 1; i < num_sample_; ++i) {
    phi_[i] = (f_(i) - f_(i + 1)) * exp(alpha_ * (i + 1 - num_sample_));
  }
  for (auto i = num_sample_; i < num_sample_ * 2; ++i) {
    phi_[i] = 0.0;
  }
}

void FastHankelTransform::evaluate_j1() {
  for (auto i = 0; i < num_sample_ * 2; ++i) {
    double x = ux_ * uy_ * x_(0) * exp(alpha_ * (i + 1 - num_sample_));
    j1_[i] = boost::math::cyl_bessel_j(1, x);
  }
}

VectorXd FastHankelTransform::calculate() {
  assert(f_updated_ == true);
  evaluate_phi();
  evaluate_j1();

  int nsample = num_sample_ * 2;
  fftw_complex *fft_phi =
      (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nsample);
  fftw_complex *phi = reinterpret_cast<fftw_complex *>(phi_);
  fftw_plan p1 =
      fftw_plan_dft_1d(nsample, phi, fft_phi, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p1);

  fftw_complex *fft_j1 =
      (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nsample);
  fftw_complex *j1 = reinterpret_cast<fftw_complex *>(j1_);
  fftw_plan p2 =
      fftw_plan_dft_1d(nsample, j1, fft_j1, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p2);

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
  fftw_execute(p3);

  VectorXd ret(num_sample_);
  for (auto i = 0; i < num_sample_; ++i) {
    ret(i) = 2.0 / (x_(i) * pow(ux_, 2) * uy_) * out[i][0] / nsample;
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
