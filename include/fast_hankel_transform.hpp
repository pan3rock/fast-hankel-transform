/*
 * File: fast_hankel_transform.hpp
 * Created Date: 2019-09-11
 * Author: Lei Pan
 * Contact: <panlei7@gmail.com>
 *
 * Last Modified: Wednesday September 25th 2019 11:38:09 am
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

#ifndef FASTHANKELTRANSFORM_FASTHANKELTRANSFORM_H_
#define FASTHANKELTRANSFORM_FASTHANKELTRANSFORM_H_

#include <Eigen/Dense>
#include <complex>

class FastHankelTransform {
public:
  FastHankelTransform(int num_sample, double ux_, double uy_);
  ~FastHankelTransform();
  Eigen::VectorXd sampling();
  void set_feval(const Eigen::Ref<const Eigen::VectorXd> &feval);
  Eigen::VectorXd calculate();

private:
  double evaluate_alpha();
  double evaluate_k0(double alpha);
  void evaluate_phi();
  void evaluate_j1();
  const int num_sample_;
  const double ux_;
  const double uy_;
  double alpha_;
  double k0_;

  Eigen::VectorXd x_;
  Eigen::VectorXd f_;
  std::complex<double> *phi_;
  std::complex<double> *j1_;

  bool x_updated_ = false;
  bool f_updated_ = false;
};

#endif