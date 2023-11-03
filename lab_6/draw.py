import numpy as np
import matplotlib.pyplot as plt
import csv

sol = []
with open('sol.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        sol.append(np.array([float(x) for x in row]))

sol = np.array(sol)

def solution(x, t):
    return np.sin(x - 0.5*t)

fig, ax = plt.subplots(1, 2, subplot_kw={"projection": "3d"})
fig.canvas.manager.set_window_title('lab 6')
ax[0].set_title('calculated\n solution')
ax[1].set_title('analitycal\n solution')

x = np.linspace(0, np.pi, 101)
t = np.linspace(0, 2, 101)
X, T = np.meshgrid(x, t)
ax[0].plot_surface(X, T, sol)
zs = np.array(solution(np.ravel(X), np.ravel(T)))
Z = zs.reshape(X.shape)
ax[1].plot_surface(X, T, Z)
plt.show()