from math import cos, exp, pi, log

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams['figure.figsize'] = (12, 12)


# аналитические решения
def analitic(x, y, t, a, mu1, mu2):
    return cos(x) * cos(y) * exp((-(mu1 ** 2 + mu2 ** 2)) * a * t)


# Граничные условия.
def u_x0(y, t, a, mu1, mu2):
    return cos(mu2 * y) * exp(-(mu1 ** 2 + mu2 ** 2) * a * t)


def u_xpi(y, t, a, mu1, mu2):
    return (-1) ** (mu1) * cos(mu2 * y) * exp(-(mu1 ** 2 + mu2 ** 2) * a * t)


def u_y0(x, t, a, mu1, mu2):
    return cos(mu1 * x) * exp(-(mu1 ** 2 + mu2 ** 2) * a * t)


def u_ypi(x, t, a, mu1, mu2):
    return (-1) ** (mu2) * cos(mu1 * x) * exp(-(mu1 ** 2 + mu2 ** 2) * a * t)


# Начальное условие.
def u_xy(x, y, a, mu1, mu2):
    return cos(mu1 * x) * cos(mu2 * y)


# Аналитическое решение на сетке

def analitic_on_net(a, dx, dy, t_, lc, mu1, mu2):
    x_count = int(pi / (2 * dx)) + 1
    y_count = int(log(2) / dy) + 1
    layers = np.zeros((lc, y_count, x_count))

    for i in range(lc):
        for j in range(y_count):
            layers[i][j] = np.array([analitic(k * dx, j * dy, i * t_, a, mu1, mu2) for k in range(x_count)])
    return layers


# погрешность
def error(explicit, analitic):
    er = max((explicit - analitic) ** 2)
    return er


# Метод переменных направлений
def variable_direction_method(a, dx, dy, t_, lc, mu1, mu2):
    x_count = int(pi / (2 * dx)) + 1
    y_count = int(log(2) / dy) + 1

    layers = np.zeros((lc, y_count, x_count))
    layers[0] = np.array([[u_xy(i * dx, j * dy, a, mu1, mu2) for i in range(x_count)] for j in range(y_count)])

    for i in range(1, lc):

        aa = -a * t_ / (2 * dx ** 2)
        bb = 1 + a * t_ / (dx ** 2)
        cc = aa

        fract_layer = np.zeros((y_count, x_count))

        for k in range(1, y_count - 1):

            alpha = np.zeros(x_count - 1)
            beta = np.zeros(x_count - 1)
            beta[0] = u_x0(k * dy, t_ * (i - 0.5), a, mu1, mu2)

            for j in range(1, x_count - 1):
                alpha[j] = -aa / (bb + cc * alpha[j - 1])
                xi_jk = layers[i - 1][k][j] + a * t_ / 2 * (
                            layers[i - 1][k + 1][j] - 2 * layers[i - 1][k][j] + layers[i - 1][k - 1][j]) / (dy ** 2)
                beta[j] = (xi_jk - cc * beta[j - 1]) / (bb + cc * alpha[j - 1])

            fract_layer[k][x_count - 1] = u_xpi(k * dy, t_ * (i - 0.5), a, mu1, mu2)

            for j in range(x_count - 2, -1, -1):
                fract_layer[k][j] = fract_layer[k][j + 1] * alpha[j] + beta[j]
        fract_layer[0] = np.array([u_y0(k * dx, t_ * (i - 0.5), a, mu1, mu2) for k in range(x_count)])
        fract_layer[y_count - 1] = np.array([u_ypi(k * dx, t_ * (i - 0.5), a, mu1, mu2) for k in range(x_count)])

        aa = -a * t_ / (2 * dy ** 2)
        bb = 1 + a * t_ / (dy ** 2)
        cc = aa

        layer = np.zeros((y_count, x_count))

        for k in range(1, x_count - 1):

            alpha = np.zeros(y_count - 1)
            beta = np.zeros(y_count - 1)
            beta[0] = u_y0(k * dx, t_ * i, a, mu1, mu2)

            for j in range(1, y_count - 1):
                alpha[j] = -aa / (bb + cc * alpha[j - 1])
                xi_jk = fract_layer[j][k] + a * t_ / 2 * (
                            fract_layer[j][k + 1] - 2 * fract_layer[j][k] + fract_layer[j][k - 1]) / (dx ** 2)
                beta[j] = (xi_jk - cc * beta[j - 1]) / (bb + cc * alpha[j - 1])

            layer[y_count - 1][k] = u_ypi(k * dx, t_ * i, a, mu1, mu2)

            for j in range(y_count - 2, -1, -1):
                layer[j][k] = layer[j + 1][k] * alpha[j] + beta[j]

        layer[:, 0] = np.array([u_x0(dy * j, t_ * i, a, mu1, mu2) for j in range(y_count)])
        layer[:, x_count - 1] = np.array([u_xpi(dy * j, t_ * i, a, mu1, mu2) for j in range(y_count)])
        layers[i] = layer

    return layers


# Метод дробных шагов
def double_steps(a, dx, dy, t_, lc, mu1, mu2):
    x_count = int(pi / (2 * dx)) + 1
    y_count = int(log(2) / dy) + 1

    layers = np.zeros((lc, y_count, x_count))
    layers[0] = np.array([[u_xy(i * dx, j * dy, a, mu1, mu2) for i in range(x_count)] for j in range(y_count)])

    for i in range(1, lc):
        aa = -a * t_ / (dx ** 2)
        bb = 1 + 2 * a * t_ / (dx ** 2)
        cc = aa

        fract_layer = np.zeros((y_count, x_count))

        for k in range(1, y_count - 1):
            alpha = np.zeros(x_count - 1)
            beta = np.zeros(x_count - 1)
            beta[0] = u_x0(k * dy, t_ * (i - 0.5), a, mu1, mu2)

            for j in range(1, x_count - 1):
                alpha[j] = -aa / (bb + cc * alpha[j - 1])
                beta[j] = (layers[i - 1][k][j] - cc * beta[j - 1]) / (bb + cc * alpha[j - 1])

            fract_layer[k][x_count - 1] = u_xpi(k * dy, t_ * (i - 0.5), a, mu1, mu2)

            for j in range(x_count - 2, -1, -1):
                fract_layer[k][j] = fract_layer[k][j + 1] * alpha[j] + beta[j]

        fract_layer[0] = np.array([u_y0(k * dx, t_ * (i - 0.5), a, mu1, mu2) for k in range(x_count)])
        fract_layer[y_count - 1] = np.array([u_ypi(k * dx, t_ * (i - 0.5), a, mu1, mu2) for k in range(x_count)])

        aa = -a * t_ / (dy ** 2)
        bb = 1 + 2 * a * t_ / (dy ** 2)
        cc = aa

        layer = np.zeros((y_count, x_count))

        for k in range(1, x_count - 1):
            alpha = np.zeros(y_count - 1)
            beta = np.zeros(y_count - 1)
            beta[0] = u_y0(k * dx, t_ * i, a, mu1, mu2)

            for j in range(1, y_count - 1):
                alpha[j] = -aa / (bb + cc * alpha[j - 1])
                beta[j] = (fract_layer[j][k] - cc * beta[j - 1]) / (bb + cc * alpha[j - 1])

            layer[y_count - 1][k] = u_ypi(k * dx, t_ * i, a, mu1, mu2)

            for j in range(y_count - 2, -1, -1):
                layer[j][k] = layer[j + 1][k] * alpha[j] + beta[j]

        layer[:, 0] = np.array([u_x0(dy * j, t_ * i, a, mu1, mu2) for j in range(y_count)])
        layer[:, x_count - 1] = np.array([u_xpi(dy * j, t_ * i, a, mu1, mu2) for j in range(y_count)])
        layers[i] = layer
    return layers


print("введите параметры а, t_, dx, dy, lc ")
a = float(input())
t_ = float(input())
dx = float(input())
dy = float(input())
lc = int(input())

mu1 = 1
mu2 = 1

real_res = analitic_on_net(a, dx, dy, t_, lc, mu1, mu2)
var_res = double_steps(a, dx, dy, t_, lc, mu1, mu2)

layers_count = lc

x_count = int(pi / (2 * dx)) + 1
y_count = int(log(2) / dy) + 1

X = np.array([i * dx for i in range(x_count)])
Y = np.array([i * dy for i in range(y_count)])

fig = plt.figure(figsize=(50, 35))
fst = fig.add_subplot(2, 3, 1)
scd = fig.add_subplot(2, 3, 2)
thd = fig.add_subplot(2, 3, 3)

y_const = int(y_count - 1)
t_const = layers_count - 1

fst.set_title('Аналитическое решение.'.format(y_const * dy), fontsize=10)
fst.set_xlabel('x', fontsize=10)
fst.set_ylabel('U(x, y, t)'.format(y_const * dy), fontsize=10)

for layer in real_res[-3:]:
    fst.plot(X, layer[y_const])
fst.legend(['t = {}'.format(i * t_) for i in range(layers_count - 3, layers_count)], fontsize=5, loc='upper right')

scd.set_title('Схема переменных направлений.'.format(y_const * dy), fontsize=7)
scd.set_xlabel('x', fontsize=10)
scd.set_ylabel('U(x, y, t)'.format(y_const * dy), fontsize=10)

for layer in real_res[-3:]:
    scd.plot(X, layer[y_const])
scd.legend(['t = {}'.format(i * t_) for i in range(layers_count - 3, layers_count)], fontsize=5, loc='upper right')

T = np.array([i * t_ for i in range(layers_count)])

MSE_error = np.array([error(i[y_const], j[y_const]) for i, j in zip(var_res, real_res)])
thd.set_title('1e-5            MSE: dx = {}, dy = {}, tau = {}'.format(dx, dy, t_), fontsize=10, pad=5, loc='left')
thd.set_xlabel('t', fontsize=10)
thd.set_ylabel('mse_error', fontsize=10)
thd.plot(T, MSE_error)
thd.scatter(T, MSE_error, s=50)
thd.legend(['Значение ошибки при'.format(y_const * dy)], fontsize=5, loc='upper right')

# plt.show()

real_res = analitic_on_net(a, dx, dy, t_, lc, mu1, mu2)
var_res = variable_direction_method(a, dx, dy, t_, lc, mu1, mu2)

layers_count = lc

x_count = int(pi / (2 * dx)) + 1
y_count = int(log(2) / dy) + 1

X = np.array([i * dx for i in range(x_count)])
Y = np.array([i * dy for i in range(y_count)])

fig1 = plt.figure(figsize=(50, 35))
fst = fig1.add_subplot(2, 3, 1)
scd = fig1.add_subplot(2, 3, 2)
thd = fig1.add_subplot(2, 3, 3)

y_const = int(y_count - 1)
t_const = layers_count - 1

fst.set_title('Аналитическое решение.'.format(y_const * dy), fontsize=10)
fst.set_xlabel('x', fontsize=10)
fst.set_ylabel('U(x, y, t)'.format(y_const * dy), fontsize=10)

for layer in real_res[-3:]:
    fst.plot(X, layer[y_const])
fst.legend(['t = {}'.format(i * t_) for i in range(layers_count - 3, layers_count)], fontsize=5, loc='upper right')

scd.set_title('Схема дробных шагов.'.format(y_const * dy), fontsize=10)
scd.set_xlabel('x', fontsize=10)
scd.set_ylabel('U(x, y, t)'.format(y_const * dy), fontsize=10)

for layer in real_res[-3:]:
    scd.plot(X, layer[y_const])
scd.legend(['t = {}'.format(i * t_) for i in range(layers_count - 3, layers_count)], fontsize=5, loc='upper right')

T = np.array([i * t_ for i in range(layers_count)])

MSE_error = np.array([error(i[y_const], j[y_const]) for i, j in zip(var_res, real_res)])
thd.set_title('1e-5            MSE: dx = {}, dy = {}, tau = {}'.format(dx, dy, t_), fontsize=10, pad=5, loc='left')
thd.set_xlabel('t', fontsize=10)
thd.set_ylabel('mse_error', fontsize=10)
thd.plot(T, MSE_error)
thd.scatter(T, MSE_error, marker='o', s=50)
thd.legend(['Значение ошибки при'.format(y_const * dy)], fontsize=5, loc='upper right')

plt.show()
