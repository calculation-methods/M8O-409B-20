#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import csv
import os
import cv2
import argparse

MU1 = 1
MU2 = 1
ALPHA_SQ = 0.1
L = 3.141593
T = 1.0
N = 6
M = 6
K = 5
H_X  = L / (N - 1)
H_Y = L / (M - 1)
TAU = T / K


parser = argparse.ArgumentParser(description='Solve PDE and draw it')

parser.add_argument("-r", "--recalc", action="store_true",
                    help="calculate solution again")

parser.add_argument("-s", "--scale", action="store",
                    help="scale of solution", default=30)
args = parser.parse_args()
IMAGE_SCALE = int(args.scale)
if(args.recalc):
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


sol = np.concatenate((U_true(0), U_calc[0]), axis=1)

#getting sollution in k-th time step
for k in range(1, K):
    a = np.concatenate((U_true(K), U_calc[k]), axis=1)
    sol = np.concatenate((sol, a), axis=0)

#save to file and resize
plt.imsave('sol.png', sol, cmap='inferno')
im = cv2.imread('sol.png')
im = cv2.resize(im, (IMAGE_SCALE*M*2, IMAGE_SCALE*K*M))
cv2.imwrite('sol.png', im)
