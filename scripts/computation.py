import math

import numpy as np
from numpy import sin, cos, power, sqrt
from scipy.special import ellipj, ellipk


class Constants:
    def __init__(self, theta0, g, length):
        self.theta0 = theta0
        self.g = g
        self.length = length


def agm(a, b):
    tolerance = 10 ** -10
    while abs(a - b) > tolerance:
        a, b = (a + b) / 2, math.sqrt(a * b)

    return a


def equation(consts: Constants, t, y0):
    theta, x = y0
    f = [x, -(consts.g / consts.length) * sin(theta)]
    return f


def exact_solve(consts: Constants, t):
    sn_n = ellipk(power(sin(consts.theta0 / 2), 2)) - sqrt(consts.g / consts.length) * t
    sn_m = power(sin(consts.theta0 / 2), 2)
    sn, cn, dn, ph = ellipj(sn_n, sn_m)
    return 2 * np.arcsin(sin(consts.theta0 / 2) * sn)


def fourier(consts: Constants, t):
    eps = (1 - sqrt(cos(consts.theta0 / 2))) / (2 + 2 * sqrt(cos(consts.theta0 / 2)))
    q = eps + 2 * power(eps, 5) + 15 * power(eps, 9) + 150 * power(eps, 13) + 1707 * power(eps, 17) + 20910 * power(
        eps, 21)
    result = 0
    x = agm(1, cos(consts.theta0 / 2)) / sqrt(consts.length / consts.g)
    for n in range(1, 4, 2):
        result += (power(-1, n // 2) * power(q, n / 2) * cos(n * x * t)) / (n * (1 + power(q, n)))
    result *= 8

    return result
