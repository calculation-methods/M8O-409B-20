import matplotlib.pyplot as plt
import numpy as np


def analyt_func(x, y):
    return np.exp(x) * np.cos(y)


def func_borderX0(y):
    return np.cos(y)


def func_borderXl(y):
    return np.e * np.cos(y)


def func_borderY0(x):  # dU/dy
    return 0


def func_borderYl(x):  # dU/dy
    return -np.exp(x)


def norm(cur_u, prev_u):
    max = 0
    for i in range(cur_u.shape[0]):
        for j in range(cur_u.shape[1]):
            if abs(cur_u[i, j] - prev_u[i, j]) > max:
                max = abs(cur_u[i, j] - prev_u[i, j])

    return max


def liebman(x, y, h, eps):
    N = len(x)
    M = len(y)
    count = 0
    prev_u = np.zeros((N, M))
    cur_u = np.zeros((N, M))
    cur_u[0] = func_borderX0(y)
    cur_u[-1] = func_borderXl(y)

    for j in range(M):
        for i in range(1, N - 1):
            cur_u[i][j] = cur_u[i][0] + (cur_u[i][-1] - cur_u[i][0]) / (x[-1] - x[0]) * (x[i] - x[0])

    while norm(cur_u, prev_u) > eps:
        count += 1
        prev_u = np.copy(cur_u)
        for i in range(1, N - 1):
            for j in range(1, M - 1):
                cur_u[i][j] = (h[0] ** 2 * (prev_u[i - 1][j] + prev_u[i + 1][j]) +
                               h[1] ** 2 * (prev_u[i][j - 1] + prev_u[i][j + 1])) / (2 * (h[0] ** 2 + h[1] ** 2))
        cur_u[:, 0] = cur_u[:, 1] - h[1] * func_borderY0(x)
        cur_u[:, -1] = cur_u[:, -2] + h[1] * func_borderYl(x)

    U = np.copy(cur_u)
    return U, count


def relaxation(x, y, h, eps, w=1.8):
    N = len(x)
    M = len(y)
    count = 0
    prev_u = np.zeros((N, M))
    cur_u = np.zeros((N, M))

    cur_u[0] = func_borderX0(y)
    cur_u[-1] = func_borderXl(y)

    for j in range(M):
        for i in range(1, N - 1):
            cur_u[i][j] = cur_u[i][0] + (cur_u[i][-1] - cur_u[i][0]) / (x[-1] - x[0]) * (x[i] - x[0])

    while norm(cur_u, prev_u) > eps:
        count += 1
        prev_u = np.copy(cur_u)
        for i in range(1, N - 1):
            for j in range(1, M - 1):
                cur_u[i][j] = (h[0]**2 * (cur_u[i-1][j] + prev_u[i+1][j]) +
                               h[1]**2 * (cur_u[i][j-1] + prev_u[i][j+1])) / (2 * (h[0]**2 + h[1]**2))
                cur_u[i][j] *= w
                cur_u[i][j] += (1 - w) * prev_u[i][j]
        cur_u[:, 0] = cur_u[:, 1] - h[1] * func_borderY0(x)
        cur_u[:, -1] = cur_u[:, -2] + h[1] * func_borderYl(x)

    U = np.copy(cur_u)
    return U, count


def Zeidel(x, y, h, eps):
    return relaxation(x, y, h, eps, w=1)


def main(n, eps):
    hx = 0.1
    hy = 0.1
    h = [hx, hy]
    x = np.arange(0, 1 + h[0] / 2 - 1e-4, h[0])
    y = np.arange(0, np.pi / 2 + h[1] / 2 - 1e-4, h[1])

    while (1):
        print("Выберите метод:\n"
              "1 - метод Либмана\n"
              "2 - метод Зейделя\n"
              "3 - метод простых итераций с верхней релаксацией\n"
              "0 - выход из программы")
        method = int(input())
        if method == 0:
            break
        if method == 1:
            U, count = liebman(x, y, h, eps)
        if method == 2:
            U, count = Zeidel(x, y, h, eps)
        if method == 3:
            U, count = relaxation(x, y, h, eps)

        X, Y = np.meshgrid(x, y)
        U_analytic = np.zeros((len(x), len(y)))
        for i in range(len(x)):
            for j in range(len(y)):
                U_analytic[i][j] = analyt_func(x[i], y[j])

        error_x = []
        error_y = []
        for i in range(len(x)):
            error_x.append(max(abs(U_analytic[:, i] - U[:, i])))
            error_y.append(max(abs(U_analytic[i, :] - U[i, :])))
        plt.title("График ошибок")
        plt.plot(x, error_y, label="При фиксированном x", color="red")
        plt.plot(x, error_x, label="При фиксированном y", color="blue")
        plt.text(0, 2, "Кол-во итераций: " + str(count))
        plt.xlabel("x, y")
        plt.ylabel("error")
        plt.grid()
        plt.legend()

        n = 2
        m = 2
        x_step = x.size // (n * m)
        y_step = y.size // (n * m)
        p_x = [k for k in range(0, x.size - 1, x_step)]
        p_y = [k for k in range(0, y.size - 1, y_step)]
        fig, ax = plt.subplots(n, m)
        fig.suptitle('Сравнение решений по y')
        fig.set_figheight(8)
        fig.set_figwidth(16)
        k = 0
        for i in range(n):
            for j in range(m):
                ax[i][j].set_title(f'Решение при x = {y[p_y[k]]}')
                ax[i][j].plot(x, U_analytic[:, p_y[k]], label='Aналитическое решение')
                ax[i][j].plot(x, U[:, p_y[k]], label='Численный метод')
                ax[i][j].grid(True)
                ax[i][j].set_xlabel('y')
                ax[i][j].set_ylabel('u')
                k += 1
        plt.legend(bbox_to_anchor=(1.05, 2), loc='upper left', borderaxespad=0.)
        fig, ax = plt.subplots(n, m)
        fig.suptitle('Сравнение решений по x')
        fig.set_figheight(8)
        fig.set_figwidth(16)
        k = 0
        for i in range(n):
            for j in range(m):
                ax[i][j].set_title(f'Решение при y = {x[p_x[k]]}')
                ax[i][j].plot(y, U_analytic[p_x[k]], label='Aналитическое решение')
                ax[i][j].plot(y, U[p_x[k]], label='Численный метод')
                ax[i][j].grid(True)
                ax[i][j].set_xlabel('x')
                ax[i][j].set_ylabel('u')
                k += 1
        plt.legend(bbox_to_anchor=(1.05, 2), loc='upper left', borderaxespad=0.)
        plt.show()

    return 0


eps = 0.001
n = 10
main(n, eps)
