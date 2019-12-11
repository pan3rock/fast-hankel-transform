'''
File: test_fht.py
Created Date: 2019-09-11
Author: Lei Pan
Contact: <panlei7@gmail.com>

Last Modified: Wednesday September 25th 2019 11:37:00 am

MIT License

Copyright (c) 2019 Lei Pan

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the 'Software'), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
-----
HISTORY:
Date      	 By	Comments
----------	---	----------------------------------------------------------
'''

from .fhtcxx import FastHankelTransform
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv


def test_alpha():
    n = 128
    def func(x): return -np.log(1 - np.exp(-x)) / (n - 1)
    plt.figure()
    t = np.linspace(-1000.0, 1000.0, 1.0e6)
    plt.plot(t, func(t))
    plt.show()
    x = 1.0
    for i in range(100):
        x_next = func(x)
        print("{:6d}{:12.5e}{:12.5e}".format(i, x, abs(x_next - x)))
        x = x_next


def func1(x):
    return np.sqrt(5.0 / 2.0 * np.pi) * x ** 2


def target(x):
    #     eta = 2.0 * np.pi * x * 10
    eta = x
    ret = np.sqrt(10.0 * np.pi) * eta ** (-4.0) * (2.0 * eta ** 2 * jv(0, eta)
                                                   + (eta ** 3 - 4.0 * eta) * jv(1, eta))
    return ret


def test_func1():
    n = int(1e4)
    ux = 1.0
    uy = 10.0
    fht = FastHankelTransform(n, ux, uy)
    x = fht.sampling()
    y = func1(x)
    x2 = x * ux
    y2 = func1(x * ux)
    fht.set_feval(y2)
    hk_eval = fht.calculate()
    plt.figure()
    plt.plot(y*uy, target(y * ux * uy), 'b', alpha=0.5)
    plt.plot(x*uy, hk_eval, 'r', alpha=0.5)
    plt.show()
