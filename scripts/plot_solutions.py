from functools import partial

import numpy as np
from matplotlib import pyplot as plt
from numpy import cos
from scipy.integrate import solve_ivp

from computation import exact_solve, equation, fourier, Constants


def plot_results(initial_angle, time, sol_rk45, sol_linear, sol_exact, sol_fourier):
    plt.plot(sol_rk45.t, sol_rk45.y[0], alpha=0.7)
    plt.plot(time, sol_linear)
    plt.plot(time, sol_exact)
    plt.plot(time, sol_fourier, alpha=0.5)

    plt.title('Pendulum Motion (Initial Angle = {} degrees)'.format(initial_angle))
    plt.xlabel('time (s)')
    plt.ylabel('angle (rad)')
    plt.grid(True)
    plt.legend(['RK45', 'linear', 'exact solution', 'fourier'], loc='lower right')
    plt.show()


def main():
    g = 10
    length = 1

    time = np.arange(0, 10.0, 0.025)

    initial_angle = 150
    theta0 = np.radians(initial_angle)
    x0 = np.radians(0.0)
    consts = Constants(theta0, g, length)

    sol_rk45 = solve_ivp(partial(equation, consts), [0, 10], [theta0, x0], t_eval=time, method='RK45')

    w = np.sqrt(g / length)
    sol_linear = [theta0 * cos(w * t) for t in time]
    sol_exact = [partial(exact_solve, consts)(t) for t in time]
    sol_fourier = [partial(fourier, consts)(t) for t in time]

    plot_results(initial_angle, time, sol_rk45, sol_linear, sol_exact, sol_fourier)


if __name__ == '__main__':
    main()
