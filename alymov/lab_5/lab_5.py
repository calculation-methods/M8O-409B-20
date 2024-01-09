from solver import Tridiagonal_matrix_algorithm as tma
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import argparse

plt.style.use("Solarize_Light2")

parser = argparse.ArgumentParser(description='Solve PDE and draw it')

parser.add_argument('--theta', type=float,
                    help=' 0 < theta < 1 - implicit/explicit; 0.5-default', default=0.5, required=False )
parser.add_argument("-d", "--draw", action="store_true",
                    help="show and save solution")
parser.add_argument("-s", "--save", action="store_true",
                    help="save solution to file")
args = parser.parse_args()

#paramters
x0 = 0
x1 = 1
t = 5
alpha = 0.01
h = 0.05
tau = 0.05
n = int((x1 - x0) / h) + 1
m = int((t) / tau) + 1
c_0 = 1
theta = args.theta

#initial conditions
def initial_conditions(x):
    return np.sin(2 * np.pi * x)

#boundary conditions
def boundary_conditions(t=t, x=0, x0=x0, x1=x1, u0=0, u1=0):
    return (x == x0) * u0 + (x == x1) * u1

#computation
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
    A[0][0] = c_0

    b[-1] = boundary_conditions(t=k*tau, x=x1)
    A[n-1][n-1] = c_0

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


#draw
fig, ax = plt.subplots(1, 2, subplot_kw={"projection": "3d"})

fig.canvas.manager.set_window_title('lab 5')
ax[0].set_title('calculated\n solution')
ax[1].set_title('analitycal\n solution')

# Make data.
X = np.linspace(x0, x1, n)
T = np.linspace(0, t, m)
X, T = np.meshgrid(X, T)
U_true = np.exp(-4*np.pi**2 * alpha * T) * np.sin(2 * np.pi * X)
U = ex_im_scheme(theta=theta)


# Plot the surface.
surf = ax[0].plot_surface(X, T, U, cmap=cm.coolwarm, linewidth=0, antialiased=False)
surf1 = ax[1].plot_surface(X, T, U_true, cmap=cm.coolwarm, linewidth=0, antialiased=False,)

if(args.draw):
    plt.show()
    plt.imsave('solutions.png', U_true)
if(args.draw):
    np.savetxt('solution.csv', U, delimiter=',')