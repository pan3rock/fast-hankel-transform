from .fhtcxx import FastHankelTransform
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv


def test_alpha():
    n = 100
    def func(x): return -np.log(1 - np.exp(-x)) / (n - 1)
    x = 1.0
    for i in range(100):
        x_next = func(x)
        print("{:6d}{:12.5e}{:12.5e}".format(i, x, abs(x_next - x)))
        x = x_next


def func1(x):
    return np.sqrt(5.0 / 2.0 * np.pi) * x ** 2


def target(x):
    eta = 2.0 * np.pi * x
    ret = np.sqrt(10.0 * np.pi) * eta ** (-4.0) * (2.0 * eta ** 2 * jv(0, eta)
                                                   + (eta ** 3 - 4.0 * eta) * jv(1, eta))
    return ret


def test_func1():
    n = int(1e2)
    fht = FastHankelTransform(n)
    x = fht.sampling()
    print(x)
    y = func1(x)
    fht.set_feval(y)
    hk_eval = fht.calculate()
    plt.figure()
    plt.plot(y, target(y), 'b', alpha=0.5)
    plt.plot(y, hk_eval, 'r', alpha=0.5)
    plt.show()
