from math import *

import matplotlib.pyplot as plt
import numpy as np


# аналитические решения
def analitic(x, t, ap):
    return exp(-4 * (pi ** 2) * ap * t) * sin(2 * pi * x)


# начальные условия
def u_x_0(x):
    return sin(2 * pi * x)


# граничные условия, в одном уравнении тк их возможно обобщить
def u_gr():
    return 0


# Аналитическое решение на сетке
def analitic_on_net(t_, lc, h, ap):
    x_c = int(pi / h) + 1
    group = np.zeros((lc, x_c))

    for i in range(lc):
        group[i] = np.array([analitic(j * h, i * t_, ap) for j in range(x_c)])
    return group


# Явная схема
def explicit_solution(t_, lc, h, a):
    x_count = int(pi / h) + 1
    group = np.zeros((lc, x_count))
    alpha = np.zeros(x_count - 1)
    beta = np.zeros(x_count - 1)

    aa = -a * t_ / (h ** 2)
    bb = 1 + 2 * a * t_ / (h ** 2)
    cc = aa  # -a * t_/ (h**2)

    for i in range(lc):
        if i == 0:
            group_ = [u_x_0(j * h) for j in range(x_count)]
        else:
            group_ = np.zeros(x_count)
            beta[0] = u_gr()
            for j in range(1, x_count - 1):
                alpha[j] = -aa / (bb + cc * alpha[j - 1])
                beta[j] = ((group[i - 1][j] - aa * (
                        group[i - 1][j + 1] - 2 * group[i - 1][j] + group[i - 1][j - 1]))) / (
                                  bb + cc * alpha[j - 1])
            group_[x_count - 1] = u_gr()
            for j in range(x_count - 2, -1, -1):
                group_[j] = group_[j + 1] * alpha[j] + beta[j]
        group[i] = np.array(group_)
    return group


# He явная схема
def implicit_solve(t_, lc, h, a):
    x_count = int(pi / h) + 1
    layers = np.zeros((lc, x_count))
    alpha = np.zeros(x_count - 1)
    beta = np.zeros(x_count - 1)

    aa = -a * t_ / (h ** 2)
    bb = 1 + 2 * a * t_ / (h ** 2)
    cc = -a * t_ / (h ** 2)

    for i in range(lc):
        if i == 0:
            layer = [u_x_0(j * h) for j in range(x_count)]
        else:
            layer = np.zeros(x_count)
            beta[0] = u_gr()
            for j in range(1, x_count - 1):
                alpha[j] = -aa / (bb + cc * alpha[j - 1])
                beta[j] = (layers[i - 1][j] - cc * beta[j - 1]) / (bb + cc * alpha[j - 1])
            layer[x_count - 1] = u_gr()
            for j in range(x_count - 2, -1, -1):
                layer[j] = layer[j + 1] * alpha[j] + beta[j]
        layers[i] = np.array(layer)
    return layers


# схема Кранка-Николсона
def KN_solve(t_, lc, h, a):
    x_count = int(pi / h) + 1
    group = np.zeros((lc, x_count))
    alpha = np.zeros(x_count - 1)
    beta = np.zeros(x_count - 1)

    aa = -a * t_ / (2 * h ** 2)
    bb = 1 + a * t_ / (h ** 2)
    cc = aa

    for i in range(lc):
        if i == 0:
            group_ = [u_x_0(j * h) for j in range(x_count)]
        else:
            group_ = np.zeros(x_count)
            beta[0] = u_gr()
            for j in range(1, x_count - 1):
                alpha[j] = -aa / (bb + cc * alpha[j - 1])
                beta[j] = ((group[i - 1][j] - aa * (
                        group[i - 1][j + 1] - 2 * group[i - 1][j] + group[i - 1][j - 1])) - cc * beta[j - 1]) / (
                                  bb + cc * alpha[j - 1])
            group_[x_count - 1] = u_gr()
            for j in range(x_count - 2, -1, -1):
                group_[j] = group_[j + 1] * alpha[j] + beta[j]
        group[i] = np.array(group_)
    return group


# погрешность
def error(explicit, analitic):
    er = sum((explicit - analitic) ** 2) / len(explicit)
    return er


def update_plot(analitic, t_, lc, h):
    cur_display_mode = display_modes[display_index]
    if cur_display_mode == 1:
        name = 'Явная схема.'
        res = explicit_solution(t_, lc, h, a)
        plot_function(name, analitic, res, t_, lc, h)
    elif cur_display_mode == 2:
        name = 'HeЯвная схема.'
        res = implicit_solve(t_, lc, h, a)
        plot_function(name, analitic, res, t_, lc, h)
    elif cur_display_mode == 3:
        name = 'Схема Кранка-Николсона.'
        res = KN_solve(t_, lc, h, a)
        plot_function(name, analitic, res, t_, lc, h)
    fig.canvas.draw()


# функция визуализации
def plot_function(name, analitic, res, time_p, lc, h):
    x_c = int(pi / h) + 1
    X = np.array([i * h for i in range(x_c)])

    for layer in res[-4:]:
        fst.plot(X, layer)
    fst.legend(['t = {}'.format(time_p - i * time_p) for i in range(4)], fontsize=5, loc='upper right')

    scd.set_title(name + 'Результаты сетки', fontsize=10)
    scd.set_xlabel('x', fontsize=5)
    scd.set_ylabel('U(x, t)', fontsize=5)

    for layer in res[-4:]:
        scd.plot(X, layer)
    scd.legend(['t = {}'.format(time_p - i * time_p) for i in range(4)], fontsize=5, loc='upper right')
    T = np.array([i * t_ for i in range(lc)])
    MSE_error = np.array([error(i, j) for i, j in zip(res, analitic)])
    thd.set_title('MSE: h = {}, tau = {}'.format(h, t_), fontsize=10)
    thd.set_xlabel('t', fontsize=5)
    thd.set_ylabel('mse_error', fontsize=5)
    thd.plot(T, MSE_error)
    thd.scatter(T, MSE_error, marker='o', c='r', s=50)
    thd.legend(['Значение ошибки в разные моменты времени'], fontsize=5, loc='upper right')


cur_display_mode = int(input("Выберите схему (1 - Явная схема, 2 - He явная схема, 3 - Схема Кранка-Николсона): "))
a = 0.5
t_ = 0.01
h = 0.1
lc = 100

print('\tdu/dt = a* d^2u/dx^2, a>0',
      'u(0,t) = 0,',
      'u(1,t)=0,',
      'u(x,0)=sin(2pix)', sep='\n\t')
print("")
o = False
analitic = analitic_on_net(t_, lc, h, a)
fig = plt.figure(figsize=(50, 35))
fst = fig.add_subplot(2, 3, 1)
scd = fig.add_subplot(2, 3, 2)
thd = fig.add_subplot(2, 3, 3)

display_modes = [1, 2, 3]
display_index = 0
# Первоначальное отображение графиков
update_plot(analitic, t_, lc, h)

# Отображение схемы
plt.show()
