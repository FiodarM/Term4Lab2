__author__ = 'fiodar'

import numpy as np


def divdiff(x, y):
    if len(x) < 1:
        print 'Error in function divdiff(): Invalid input'
        return
    if len(x) == 1:
        return y[0]

    return (divdiff(x[1:], y[1:]) - divdiff(x[:-1], y[:-1])) \
           / (x[-1] - x[0])


def eval_divdiffs(points):
    X, Y = points
    divdiffs = np.array([])
    for i, node in enumerate(X):
        divdiffs = np.append(divdiffs, divdiff(X[:i + 1], Y[:i + 1]))

    return divdiffs


def newton_p(points, x, divdiffs=None):
    X, Y = points

    if divdiffs is None:
        divdiffs = eval_divdiffs(points)

    f = divdiffs[-1]
    for i, node in enumerate(X[:0:-1]):
        f = divdiffs[-2 - i] + (x - X[-2 - i]) * f

    return f


def cubic_spline(points):
    X, Y = np.array(points)[0], np.array(points)[1]

    coefs = np.transpose(eval_abcd(X, Y))
    for i, item in enumerate(X[:-1]):
        print 'spline_{0} = {1:.3f} {2:.3f}x {3:.3f}x**2 {4:.3f}x**3'.format(i + 1, *coefs[i])

    delta = lambda x: [np.greater_equal(x, X[k]) *
                       np.less(x, X[k + 1]) + (node == X[-2]) * np.greater_equal(x, X[-1])
                       + (node == X[0]) * np.less(x, X[0]) for k, node in enumerate(X[:-1])]
    i = lambda x: np.nonzero(delta(x))[0]
    coef = lambda x: coefs[i(x)].transpose()
    spline = lambda x: coef(x)[0] + (x - X[i(x)]) * \
                                    (coef(x)[1] + (x - X[i(x)]) * (coef(x)[2] + (x - X[i(x)]) * coef(x)[3]))
    return spline


def eval_abcd(X, Y):

    h, dy = X[1:] - X[:-1], Y[1:] - Y[:-1]
    h = np.array(h)
    dy = np.array(dy)
    d_sub = h[:-1]
    d_main = 2. * (h[:-1] + h[1:])
    d_super = h[1:]
    data = (d_sub, d_main, d_super)
    f = 3. * (dy[1:] / h[1:] - dy[:-1] / h[:-1])
    a = Y[:-1]
    c = [0., 0.]
    c = np.insert(c, 1, tridiag_solve(data, f))
    b = (dy / h) - h * (c[1:] + 2. * c[:-1]) / 3.
    d = (c[1:] - c[:-1]) / (3. * h)

    return a, b, c[:-1], d


def tridiag_solve(diags, f):

    alpha, beta, gamma = diags
    x = np.zeros_like(f)
    k, l = np.zeros_like(f), np.zeros_like(f)
    k[0] = f[0] / beta[0]
    l[0] = - gamma[0] / beta[0]
    for i in range(1, len(f) - 1):
        k[i] = (f[i] - alpha[i]*k[i-1])/(alpha[i]*l[i-1] + beta[i])
        l[i] = - gamma[i] / (alpha[i]*l[i-1] + beta[i])

    x[-1] = (f[-1] - alpha[-1]*k[-2]) / (alpha[-1]*l[-2] + beta[-1])
    for i, item in enumerate(x[-2:0:-1]):
        x[-2-i] = (f[-2-i] - alpha[-2-i]*k[-3-i] - gamma[-2-i]*x[-1-i]) /\
                  (alpha[-2-i]*l[-3-i] + beta[-2-i])

    x[0] = k[0] + l[0] * x[1]

    return x