import numpy as np
import matplotlib.pyplot as plt
from math import sqrt


def psi(x):
    return 0

def dpsidt(x):
    return np.exp(-x)*np.sin(x)

def dpsidx(x):
    return 0

def dpsidxx(x):
    return 0


def phi0(t):
    return 0


def phi1(t):
    return 0


def true_fval(x, t):
    return 0.5*np.exp(-x)*np.sin(x)*np.sin(2*t)



def norma(a):
    norm = 0
    for i in range(len(a)):
        norm += a[i] ** 2
    return sqrt(norm)


# метод прогонки
def tridiagonal(a, b, c, d):
    n = len(d)
    x = np.zeros(n)
    p = [-c[0] / b[0]]
    q = [d[0] / b[0]]
    for i in range(1, n):
        p.append(-c[i] / (b[i] + a[i] * p[i - 1]))
        q.append((d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1]))
    x[-1] = q[-1]
    for i in reversed(range(n - 1)):
        x[i] = p[i] * x[i + 1] + q[i]
    return x


def ExScheme(a,b, lb, ub, h, tau, T,apr):
    # разбиение осей
    x = np.arange(lb, ub + h, h)
    t = np.arange(0, T + tau, tau)
    # строим конечно-разностную сетку
    U = np.zeros((len(t), len(x)))
    # заполним первый уровень
    for j in range(len(x)):
        U[0, j] = psi(x[j])
    #заполним второй уровень
    if apr == 1:
        for j in range (len(x)):
            U[1,j] = U[0, j]+ tau*dpsidt(x[j])
    if apr == 2:
        for j in range (len(x)):
            U[1,j] = U[0, j]+ tau*dpsidt(x[j])+(tau**2)*dpsidxx(x[j])+2*(tau**2)*dpsidx(x[j])
    # прямая схема
    for i in range(2, len(t)):
        for j in range(1, len(x)-1):
            U[i, j] = a*(tau**2)/(h**2)*U[i-1,j+1]+(2-2*a*(tau**2)/(h**2))*U[i-1,j]+a*(tau**2)/(h**2)*U[i-1,j-1]+b*(tau**2)/(2*h)*U[i-1,j+1]-b*(tau**2)/(2*h)*U[i-1,j-1]-U[i-2,j]
            U[i, 0] = phi0(t[i])
            U[i, -1] = phi1(t[i])

    return U



def ImScheme(a,b, lb, ub, h, tau, T,apr):
    x = np.arange(lb, ub + h, h)
    t = np.arange(0, T + tau, tau)
    U = np.zeros((len(t), len(x)))
     # заполним первый уровень
    for j in range(len(x)):
        U[0, j] = psi(x[j])
    #заполним второй уровень
    if apr == 1:
        for j in range (len(x)):
            U[1,j] = U[0, j]+ tau*dpsidt(x[j])
    if apr == 2:
        for j in range (len(x)):
            U[1,j] = U[0, j]+ tau*dpsidt(x[j])+(tau**2)*dpsidxx(x[j])+2*(tau**2)*dpsidx(x[j])

    for i in range(2, len(t)):
        aa = np.zeros(len(x)-2)
        bb = np.zeros(len(x)-2)
        cc = np.zeros(len(x)-2)
        dd = np.zeros(len(x)-2)
        dd[0] = -2*U[i-1,1]+U[i-2,1] - (a * (tau**2) / (h ** 2) - b*(tau**2)/(2*h)) * phi0(t[i])
        dd[-1] = -2*U[i-1,len(x)-1]+U[i-2,len(x)-1] - (a * (tau**2) / (h ** 2) + b*(tau**2)/(2*h)) * phi1(t[i])
        bb[0] = -(1 + 2 * a * (tau**2) / (h ** 2))
        bb[-1] = -(1 + 2 * a * (tau**2) / (h ** 2))
        cc[0] = a * (tau**2) / (h ** 2) + b*(tau**2)/(2*h)
        aa[-1] = a * (tau**2) / (h ** 2) - b*(tau**2)/(2*h)
        for j in range(1, len(x) - 2):
            aa[j] = a * (tau**2) / (h ** 2) - b*(tau**2)/(2*h)
            bb[j] = -(1 + 2 * a * (tau**2) / (h ** 2))
            cc[j] = a * (tau**2) / (h ** 2) + b*(tau**2)/(2*h)
            dd[j] = -2*U[i - 1, j+1] + U[i-2,j+1]
        xx = tridiagonal(aa, bb, cc, dd)
        for j in range(1,len(x)-1):
            U[i, j] = xx[j-1]

    return U



def plot_ex(a,b, lb, ub, h, tau, T,k,apr):
    x = np.arange(lb, ub + h, h)
    t = np.arange(0, T + tau, tau)
    sigma = a*(tau**2)/(h**2)
    plt.figure(1)
    plt.title('Явная схема, t = ' + str(t[k])+ ' sigma = ' + str(sigma))
    plt.grid()
    plt.plot(x, true_fval(x, t[k]), color='red', label='аналитическое решение')
    U = ExScheme(a,b, lb, ub, h, tau, T,apr)
    plt.plot(x, U[k, :], color='blue', label='численное решение')
    plt.legend()
    plt.xlim((0, ub))
 
    plt.figure(3)
    plt.title('с.к.погрешность от t')
    plt.grid()
    eps = []
    for i in range(len(t)):
        a = true_fval(x,t[i]) - U[i, :]
        eps = np.append(eps, norma(a))
    plt.plot(t, eps, color='blue')

    plt.show()
    return

def plot_im(a,b, lb, ub, h, tau, T,k,apr):
    x = np.arange(lb, ub + h, h)
    t = np.arange(0, T + tau, tau)
    sigma = a*(tau**2)/(h**2)
    plt.figure(1)
    plt.title('Невная схема, t = ' + str(t[k])+ ' sigma = ' + str(sigma))
    plt.grid()
    plt.plot(x, true_fval(x, t[k]), color='red', label='аналитическое решение')
    U = ImScheme(a,b, lb, ub, h, tau, T,apr)
    plt.plot(x, U[k, :], color='blue', label='численное решение')
    plt.legend()
    plt.xlim((0, ub))

    plt.figure(3)
    plt.title('с.к.погрешность от t')
    plt.grid()
    eps = []
    for i in range(len(t)):
        a = true_fval(x, t[i]) - U[i, :]
        eps = np.append(eps, norma(a))
    plt.plot(t, eps, color='blue')

    plt.show()
    return



plot_ex(2,4,0,np.pi,0.05,0.00005,1,10000,2)
plot_im(2,4,0,np.pi,0.005,0.00005,1,10000,2)



