from solver import Tridiagonal_matrix_algorithm as tma
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
plt.style.use("Solarize_Light2")


x0 = 0
x1 = 1
t = 2
alpha = 0.05
h = 0.1
tau = 0.1
n = int((x1 - x0) / h) + 1
m = int((t) / tau) + 1


def initial_conditions(x):
    return np.sin(2 * np.pi * x)

def boundary_conditions(t=t, x=0, x0=x0, x1=x1, u0=0, u1=0):
    return (x == x0) * u0 + (x == x1) * u1


def ex_step(u, k):
    l = []
    l.append(boundary_conditions(t=k*tau, x=x0))
    for j in range(1, n-1):
        u_c = (alpha * tau / h**2) * (u[j - 1] - 2 * u[j] + u[j + 1]) + u[j]
        l.append(u_c)
    l.append(boundary_conditions(t=k*tau, x=x1))
    return np.array(l)

def explcit_scheme():
    x = np.linspace(x0, x1, n)
    u = []
    u.append(initial_conditions(x))
    for k in range(1, m):
        u_c = ex_step(u[-1], k = k)
        u.append(u_c)
    return np.array(u)


def im_step(u, k):

    n = len(u)
    A = np.zeros((n, n))
    b = np.zeros(n)
    b[0] = boundary_conditions(t=k*tau, x=x0)
    A[0][0] = 1
    b[-1] = boundary_conditions(t=k*tau, x=x1)
    A[n-1][n-1] = 1

    for j in range(1,n-1):
        A[j][j-1] = alpha / h**2
        A[j][j] = -1 / tau - 2 * alpha / h**2
        A[j][j+1] = alpha / h**2
        b[j] = -u[j] / tau

    return np.array(tma(A, b))

def implicit_scheme(h, tau, x0, x1, t):
    x = np.linspace(x0, x1, int((x1 - x0) / h) + 1)
    u = []
    u.append(initial_conditions(x))

    for k in range(1, m):
        u.append(np.array(im_step(u[-1], k=k)))
    return np.array(u)


def ex_im_scheme(theta=0.5):
    x = np.linspace(x0, x1, int((x1 - x0) / h) + 1)
    u = []
    u.append(initial_conditions(x))

    for k in range(1, m):
        u.append(theta * np.array(im_step(u[-1], k=k)) + (1 - theta) * np.array(ex_step(u[-1], k=k)))
    return np.array(u)


fig, ax = plt.subplots(1, 2, subplot_kw={"projection": "3d"})

fig.canvas.set_window_title('lab 5')

# Make data.
X = np.linspace(x0, x1, n)
T = np.linspace(0, t, m)
X, T = np.meshgrid(X, T)
U_true = np.exp(-4*np.pi**2 * alpha * T) * np.sin(2 * np.pi * X)
U = ex_im_scheme(theta=0.5)


ax[0].set_title('calculated\n solution')
ax[1].set_title('analitycal\n solution')
# Plot the surface.
surf = ax[0].plot_surface(X, T, U, cmap=cm.coolwarm, linewidth=0, antialiased=False)
surf1 = ax[1].plot_surface(X, T, U_true, cmap=cm.coolwarm, linewidth=0, antialiased=False,)

plt.show()