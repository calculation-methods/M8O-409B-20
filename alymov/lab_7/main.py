#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import csv
import argparse
import subprocess

plt.style.use("Solarize_Light2")

parser = argparse.ArgumentParser(description='Solve PDE and draw it')

parser.add_argument("-s", "--sol", action="store_true",
                    help="resolve eq")
args = parser.parse_args()

if(args.sol):
    subprocess.call(["./build/prog"])

sol = []
with open('sol.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        sol.append(np.array([float(x) for x in row]))

sol = np.array(sol)

def solution(x, y):
    return x + y

fig, ax = plt.subplots(1, 2, subplot_kw={"projection": "3d"})
fig.canvas.manager.set_window_title('lab 6')
ax[0].set_title('calculated\n solution')
ax[1].set_title('analitycal\n solution')

x = np.linspace(0, 1, 11)
t = np.linspace(0, 1, 11)
X, T = np.meshgrid(x, t)
ax[0].plot_surface(X, T, sol)
zs = np.array(solution(np.ravel(X), np.ravel(T)))
Z = zs.reshape(X.shape)
ax[1].plot_surface(X, T, Z)

plt.savefig("lab_7")
plt.show()