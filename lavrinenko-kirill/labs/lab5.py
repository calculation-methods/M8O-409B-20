import math
import matplotlib.pyplot as plt
import numpy as np
import ipywidgets as widgets
from IPython.display import display

N = 10
K = 500
T = 5
l = (np.pi) / 2
h = l/N
tau = T/K
sigma = tau/(h ** 2)
x = np.linspace(0,l,N)
t = np.linspace(0,T,K)
Xp, Tp = np.meshgrid(x, t)

if sigma <= 0.5:
    print("OK")
else:
    print(f'Kurant > 1/2: {sigma}')

def psi(x):
    return 0
def f(x, t):
    return np.cos(x)*(np.cos(t)+np.sin(t))
def phi0(t):
    return np.sin(t)
def phi1(t):
    return -np.sin(t)
def solution(x,t):
    return np.sin(t)*np.cos(x)

def analitic():
    u = [0]*K
    for i in range(K):
        u[i] = [0]*N
    for i in range(K):
        for j in range(N):
            u[i][j] = solution(x[j], t[i])
    return u

def update_plot(index, u_array, t_array):
    plt.figure(figsize=(10, 6))
    plt.plot(x, u_array[index])
    plt.xlabel('x')
    plt.ylabel('u')
    plt.title(f'u от x при t = {t_array[index]:.2f}')
    plt.grid()
    plt.xlim(0, (np.pi) / 2)
    plt.show()

u = analitic()
index_slider = widgets.IntSlider(value=95, min=0, max=len(u) - 1, description='Index u:')
widgets.interactive(update_plot, index=index_slider, u_array=widgets.fixed(u), t_array=widgets.fixed(t))

def explicit_solve(l, N, K, T, app):
    u = [0]*K
    for i in range(K):
        u[i] = [0]*N
    for j in range(N):
        u[0][j] = psi(j * h)
    for k in range(K):
        u[k][0] = phi0(tau * k)

    for k in range(1, K):
        for j in range(1, N-1):
            u[k][j] = u[k-1][j]+tau*(u[k-1][j-1]-2*u[k-1][j]+u[k-1][j+1])/h**2+tau*f(j*h, tau*k)
    if app == 1:
        for k in range(K):
            u[k][-1] = phi1(tau * k)*h + u[k][-2]
    if app == 2:
        for k in range(K):
            u[k][-1] = (phi1(k * tau) * 2 * h - u[k][-3] + 4 * u[k][-2]) / 3
    if app == 3:
        for k in range(K):
            u[k][-1] = (phi1(k * tau) + u[k][-2] / h + 2 * tau * u[k - 1][-1] / h) / (1 / h + 2 * tau / h)

    return u

exp1 = explicit_solve(l, N, K, T, 1)
index_slider = widgets.IntSlider(value=95, min=0, max=len(exp1) - 1, description='Index u:')
widgets.interactive(update_plot, index=index_slider, u_array=widgets.fixed(exp1), t_array=widgets.fixed(t))

exp2 = explicit_solve(l, N, K, T, 2)
index_slider = widgets.IntSlider(value=95, min=0, max=len(exp2) - 1, description='Index u:')
widgets.interactive(update_plot, index=index_slider, u_array=widgets.fixed(exp2), t_array=widgets.fixed(t))

exp3 = explicit_solve(l, N, K, T, 3)
index_slider = widgets.IntSlider(value=95, min=0, max=len(exp3) - 1, description='Index u:')
widgets.interactive(update_plot, index=index_slider, u_array=widgets.fixed(exp3), t_array=widgets.fixed(t))

def tma(a, b, c, d):
    n = len(a)
    p, q = [], []
    p.append(-c[0] / b[0])
    q.append(d[0] / b[0])
    for i in range(1, n):
        p.append(-c[i] / (b[i] + a[i] * p[i - 1]))
        q.append((d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1]))
    x = [0 for _ in range(n)]
    x[n - 1] = q[n - 1]
    for i in range(n-2, -1, -1):
        x[i] = p[i] * x[i+1] + q[i]
    return x

def implicit_solve(l, N, K, T, app):
    a = np.zeros(N)
    b = np.zeros(N)
    c = np.zeros(N)
    d = np.zeros(N)
    u = [0]*K
    for i in range(K):
        u[i] = [0]*N
    for j in range(N):
        u[0][j] = psi(j * h)
    for k in range(K):
        u[k][0] = phi0(tau*k)
    for k in range(1, K):
        a[0] = 0
        b[0] = -(1 + 2 * sigma)
        c[0] = sigma
        d[0] = -u[k-1][0]-sigma*phi0((k)*tau)
        for j in range(1, N):
            a[j] = sigma
            b[j] = -(1 + 2 * sigma)
            c[j] = sigma
            d[j] = -u[k-1][j] - tau * f(j * h, (k-1) * tau)
        if app == 1:
            a[-1] = sigma
            b[-1] = -(1 + 2 * sigma)
            c[-1] = 0
            d[-1] = -(u[k][-1] - 2 * tau * phi1(k * tau) / h)
        elif app == 2:
            a[-1] = sigma
            b[-1] = -(1 + 2 * sigma)
            c[-1] = 0
            d[-1] = -h*phi1(tau*(k))*h-u[k][-1]-tau*f(N*h,tau*(k+1))
        elif app == 3:
            a[-1] = sigma
            b[-1] = -(1 + 2 * sigma)
            c[-1] = 0
            d[-1] = -((1 - sigma) * u[k - 1][-1] + sigma / 2 * u[k - 1][-2]) - sigma * phi1(k * tau)
        u[k] = tma(a, b, c, d)
    return u

imp = implicit_solve(l, N, K, T, 1)
index_slider = widgets.IntSlider(value=95, min=0, max=len(imp) - 1, description='Index u:')
widgets.interactive(update_plot, index=index_slider, u_array=widgets.fixed(imp), t_array=widgets.fixed(t))

imp = implicit_solve(l, N, K, T, 2)
index_slider = widgets.IntSlider(value=95, min=0, max=len(imp) - 1, description='Index u:')
widgets.interactive(update_plot, index=index_slider, u_array=widgets.fixed(imp), t_array=widgets.fixed(t))

imp = implicit_solve(l, N, K, T, 3)
index_slider = widgets.IntSlider(value=95, min=0, max=len(imp) - 1, description='Index u:')
widgets.interactive(update_plot, index=index_slider, u_array=widgets.fixed(imp), t_array=widgets.fixed(t))

def crank_nicholson_solve(l, N, K, T, app):
    omega = 0.5
    a = np.zeros(N)
    b = np.zeros(N)
    c = np.zeros(N)
    d = np.zeros(N)
    u = [0]*K
    for i in range(K):
        u[i] = [0]*N
    for j in range(N):
        u[0][j] = psi(j * h)
    for k in range(K):
        u[k][0] = phi0(tau*k)
    for k in range(1, K):
        a[0] = 0
        b[0] = -(1 + 2 * sigma)
        c[0] = sigma
        d[0] = -u[k-1][0]-sigma*phi0((k)*tau)
        for j in range(1, N):
            a[j] = sigma
            b[j] = -(1 + 2 * sigma)
            c[j] = sigma
            d[j] = -u[k-1][j] - tau * f(j * h, (k-1) * tau)
        if app == 1:
            a[-1] = sigma
            b[-1] = -(1 + 2 * sigma)
            c[-1] = 0
            d[-1] = -(u[k][-1] - 2 * tau * phi1(k * tau) / h)
        elif app == 2:
            a[-1] = sigma
            b[-1] = -(1 + 2 * sigma)
            c[-1] = 0
            d[-1] = -h*phi1(tau*(k))*h-u[k][-1]-tau*f(N*h,tau*(k+1))
        elif app == 3:
            a[-1] = sigma
            b[-1] = -(1 + 2 * sigma)
            c[-1] = 0
            d[-1] = -((1 - sigma) * u[k - 1][-1] + sigma / 2 * u[k - 1][-2]) - sigma * phi1(k * tau)
        u[k] = tma(a, b, c, d)
    return u

kn = crank_nicholson_solve(l, N, K, T, 3)

index_slider = widgets.IntSlider(value=95, min=0, max=len(kn) - 1, description='Index u:')
widgets.interactive(update_plot, index=index_slider, u_array=widgets.fixed(kn), t_array=widgets.fixed(t))

def calculate_error(analytic, num):
    errors = []
    for i in range(len(analytic)):
        max_error = np.max(np.abs(np.array(analytic[i]) - np.array(num[i])))
        errors.append(max_error)
    return errors

analytic_result = analitic()
explicit_result = explicit_solve(l, N, K, T, 3)
implicit_result = implicit_solve(l, N, K, T, 3)
kn_res = kn

# Вычисление ошибки между результатами
errors1 = calculate_error(analytic_result, explicit_result)
errors2 = calculate_error(analytic_result, implicit_result)
errors3 = calculate_error(analytic_result, kn_res)

# Построение графика ошибки от времени
time_steps = np.linspace(0, T, K)  # Массив временных шагов
plt.figure(figsize=(10, 6))
plt.plot(time_steps, errors1, label='Явный метод')
plt.plot(time_steps, errors2, label='Неявный метод')
plt.plot(time_steps, errors3, label='Метод Кранка-Николсона')
plt.xlabel('Время')
plt.ylabel('Максимальная ошибка')
plt.title('Зависимость ошибки от времени')
plt.grid()
plt.legend()
plt.show()

