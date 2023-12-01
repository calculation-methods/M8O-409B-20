import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import csv
import os
import cv2

MU1 = 1
MU2 = 1
ALPHA_SQ = 0.1
L = 3.141593
T = 1.0
N = 6
M = 6
K = 12
H_X  = L / (N - 1)
H_Y = L / (M - 1)
TAU = T / K

os.system("make prog")
os.system("./prog")

U_calc = []
with open('sol.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        U_calc.append(np.array([float(x) for x in row]))

U_calc = np.array(U_calc)
U_calc = U_calc.reshape(K, M, N)

x = np.linspace(0, L, N)
y = np.linspace(0, L, M)

X, Y = np.meshgrid(x, y)

U_true = lambda k: np.cos(MU1 * X) * np.cos(MU2 * Y) * np.exp(-(MU1**2 + MU2**2) * ALPHA_SQ * TAU * k)

print(U_true(4))
fig, ax = plt.subplots(1, 2, subplot_kw={"projection": "3d"})
fig.canvas.manager.set_window_title('lab 5')
ax[0].set_title('calculated\n solution')
ax[1].set_title('analitycal\n solution')

sol = np.concatenate((U_true(0), U_calc[0]), axis=1)
# Plot the surface.
for k in range(1, K):
    #surf = ax[0].plot_surface(X, Y, U_calc[k], cmap=cm.coolwarm, linewidth=0, antialiased=False)
    #surf1 = ax[1].plot_surface(X, Y, U_true(k), cmap=cm.coolwarm, linewidth=0, antialiased=False,)
    
    a = np.concatenate((U_true(K), U_calc[k]), axis=1)
    sol = np.concatenate((sol, a), axis=0)

plt.imsave('sol.png', sol, cmap='inferno')
im = cv2.imread('sol.png')
im = cv2.resize(im, (10*M*2, 10*K*M))
cv2.imwrite('sol.png', im)
