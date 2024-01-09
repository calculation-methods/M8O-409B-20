import matplotlib.pyplot as plt
import numpy as np


def analyt_func(x, t):
    return np.cos(x) * np.sin(2 * t)


def func_border1(t):
    return np.sin(2 * t)


def func_border2(t):
    return -np.sin(2 * t)


def run_through(a, b, c, d, s):
    P = np.zeros(s + 1)
    Q = np.zeros(s + 1)

    P[0] = -c[0] / b[0]
    Q[0] = d[0] / b[0]

    k = s - 1
    for i in range(1, s):
        P[i] = -c[i] / (b[i] + a[i] * P[i - 1])
        Q[i] = (d[i] - a[i] * Q[i - 1]) / (b[i] + a[i] * P[i - 1])
    P[k] = 0
    Q[k] = (d[k] - a[k] * Q[k - 1]) / (b[k] + a[k] * P[k - 1])

    x = np.zeros(s)
    x[k] = Q[k]

    for i in range(s - 2, -1, -1):
        x[i] = P[i] * x[i + 1] + Q[i]

    return x


def explicit(K, t, tau, h, x):
    N = len(x)
    U = np.zeros((K, N))
    t += tau
    for j in range(N):
        U[0, j] = 0

        U[1][j] = 2 * np.cos(x[j]) * tau

    for k in range(1, K - 1):
        t += tau
        for j in range(1, N - 1):
            U[k + 1, j] = (U[k, j + 1] * (tau ** 2 / h ** 2) + U[k, j] * (-2 * tau ** 2 / h ** 2 + 2 - 3 * tau ** 2)
                           + U[k, j - 1] * tau ** 2 / h ** 2 - U[k - 1, j])

        U[k + 1, 0] = func_border1(t)
        U[k + 1, N - 1] = func_border2(t)

    return U


def implicit(K, t, tau, h, x):
    N = len(x)
    U = np.zeros((K, N))
    t += tau
    for j in range(N):
        U[0, j] = 0

        U[1, j] = 2 * np.cos(x[j]) * tau

    for k in range(1, K - 1):
        a = np.zeros(N)
        b = np.zeros(N)
        c = np.zeros(N)
        d = np.zeros(N)
        t += tau

        for j in range(1, N - 1):
            a[j] = 1 / h ** 2
            b[j] = -2 / h ** 2 - 1 / tau ** 2
            c[j] = 1 / h ** 2
            d[j] = U[k, j] * (3 - 2 / tau ** 2) + U[k - 1, j] / tau ** 2

        b[0] = 1
        c[0] = 0
        d[0] = func_border1(t)

        a[N - 1] = 0
        b[N - 1] = 1
        d[N - 1] = func_border2(t)

        u_new = run_through(a, b, c, d, N)
        for i in range(N):
            U[k + 1, i] = u_new[i]

    return U


N = 25
K = 200
time = 3

h = np.pi / N
tau = time / K
print(h, tau)
x = np.arange(0, np.pi + h / 2 - 1e-4, h)
T = np.arange(0, time, tau)
t = 0

while True:
    print("Выберите метод:\n"
          "(1) Явная конечно-разностная схема\n"
          "(2) Неявная конечно-разностная схема\n"
          "(0) Выход из программы")
    method = int(input("=> "))
    if method == 0:
        break
    else:
        if method == 1:
            if tau / h ** 2 <= 1:
                print("Условие Куррента выполнено:", tau / h ** 2, "<= 1\n")
                U = explicit(K, t, tau, h, x)
            else:
                print("Условие Куррента не выполнено:", tau / h ** 2, "> 1")
                break
        elif method == 2:
            U = implicit(K, t, tau, h, x)

    dt = int(input("Введите момент времени: "))
    print()

    U_analytic = analyt_func(x, T[dt])
    error = abs(U_analytic - U[dt, :])
    plt.title("График точного и численного решения задачи")
    plt.plot(x, U_analytic, label="Точное решение", color="red")
    plt.scatter(x, U[dt, :], label="Численное решение")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.text(0.2, -0.8, "Максимальная ошибка метода: " + str(max(error)))
    # plt.axis([-0.2, 3.3, -0.005, 0.005])
    plt.grid()
    plt.legend()
    plt.show()

    plt.title("График зависимости ошибок по времени и пространству")
    error_time = np.zeros(len(T))
    for i in range(len(T)):
        error_time[i] = max(abs(analyt_func(x, T[i]) - U[i, :]))
    plt.plot(T, error_time, label="По времени")
    plt.plot(x, error, label="По пространству в выбранный момент времени")
    plt.legend()
    plt.grid()
    plt.show()
