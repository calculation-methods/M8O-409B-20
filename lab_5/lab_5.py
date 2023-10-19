from solver import Tridiagonal_matrix_algorithm as tma
import numpy as np


def initial_conditions(x):
    return np.sin(2 * np.pi * x)

def boundary_conditions(t, x, b0=0, b1=1, u0=0, u1=0):
    return (x == b0) * u0 + (x == b1) * u1

'''
def explcit_scheme(h, tau, b0, b1, t):
    alpha = 0.01
    n = int((b1-b0) / h)
    m = int(t / tau)
    x = np.arange(b0, b1, h)
    u = []
    u.append(initial_conditions(x))

    for k in range(1, m):

        l = []
        l.append(boundary_conditions(t=k*tau, x=b0))
        for j in range(1, n-1):
            u_c = (alpha * tau / h**2) * (u[k-1][j - 1] - 2 * u[k-1][j] + u[k-1][j + 1]) + u[k-1][j]
            l.append(u_c)
        l.append(boundary_conditions(t=k*tau, x=b1))
        u.append(l)
    return np.array(u)


import matplotlib.pyplot as plt

from matplotlib import cm
from matplotlib.ticker import LinearLocator

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

b0=0
b1=1
h = 0.05
tau = 0.05
t = 2

# Make data.
X = np.arange(b0, b1, h)
T = np.arange(0, t, tau)
X, T = np.meshgrid(X, T)
U = explcit_scheme(h=h, tau=tau, b0=b0, b1=b1, t=t)

# Plot the surface.
surf = ax.plot_surface(X, T, U, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

plt.show()
plt.imsave()
'''

def implicit_scheme(h, tau, b0, b1, t):
    alpha = 0.01
    n = int((b1-b0) / h) + 1
    m = int(t / tau) + 1
    x = np.linspace(b0, b1, int((b1 - b0) / h) + 1)
    u = []
    u.append(initial_conditions(x))
    for k in range(1, m):
        A = np.zeros((n, n))
        b = np.zeros(n)
        b[0] = boundary_conditions(t=k*tau, x=b0)
        A[0][0] = 1
        b[-1] = boundary_conditions(t=k*tau, x=b1)
        A[n-1][n-1] = 1
        print(len(u[-1]))
        for j in range(1,n-1):
            A[j][j-1] = alpha / h**2
            A[j][j] = -1 / tau - 2 * alpha / h**2
            A[j][j+1] = alpha / h**2
            b[j] = -u[-1][j] / tau
        u.append(np.array(tma(A, b)))
        
    return np.array(u)


import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

b0=0
b1=1
h = 0.05
tau = 0.05
t = 2


# Make data.
X = np.linspace(b0, b1, int((b1 - b0) / h) + 1)
T = np.linspace(0, t, int((t) / tau) + 1)
X, T = np.meshgrid(X, T)
U = implicit_scheme(h=h, tau=tau, b0=b0, b1=b1, t=t)


# Plot the surface.
surf = ax.plot_surface(X, T, U, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

plt.show()