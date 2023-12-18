import numpy as np
from numpy.linalg import norm
from scipy.sparse import csc_matrix
from time import time


def read_matrix(filename: str):
    with open(filename) as f:
        shape = int(f.readline())
        matrix = [[float(num) for num in line.split()]
                  for _, line in zip(range(shape), f)]
        matrix = csc_matrix(matrix)
        b = np.array([float(num) for num in f.readline().split()])
        return matrix, b


class Solver:
    def __init__(self, matrix, b, x0=None, eps=1e-5):
        self.matrix = matrix
        self.b = b
        self.eps = eps
        self.shape = matrix.shape[0]
        self.x0 = np.array([0] * self.shape) if x0 is None else x0
        self.k = 0

    def solve(self, max_iter=100000):
        x0 = self.x0
        r0 = self.b - self.matrix.dot(x0)
        p0 = np.copy(r0)
        for _ in range(max_iter):
            temp = self.matrix @ p0
            norm_0 = np.dot(r0, r0)
            alpha_i = norm_0 / (temp @ p0)
            x_new = x0 + p0 * alpha_i
            r_new = r0 - temp * alpha_i
            norm_new = r_new @ r_new
            beta_i = norm_new/norm_0
            p_new = r_new + p0 * beta_i

            r0 = r_new
            p0 = p_new
            x0 = x_new

            self.k += 1
            if norm(r_new) < self.eps:
                break

        return x0

    def solve_and_print(self):
        start = time()
        x = self.solve()
        end = time()
        start2 = time()
        x2 = np.linalg.solve(self.matrix.toarray(), self.b)
        end2 = time()
        print('Метод сопряженных градиентов:\n')
        print(f'Ответ:\nX = {x.round(6)}\n')
        print(f'eps = {self.eps}; Размерность задачи = {self.shape}\n'
              f'Количество итераций = {self.k}; Время решения = {round(end - start, 8)} сек.\n')
        print('Встроенный метод NumPy:\n')
        print(f'Ответ:\nX = {x2.round(6)}\n')
        print(f'Время решения = {round(end2 - start2, 8)} сек.\n')


epsilon = float(input("Введите погрешность => "))
print()
print("Размерность задачи 10х10:")
mat10, B10 = read_matrix("matrix10.txt")
solver = Solver(mat10, B10, eps=epsilon)
solver.solve_and_print()

print("Размерность задачи 50х50:")
mat50, B50 = read_matrix("matrix50.txt")
solver = Solver(mat50, B50, eps=epsilon)
solver.solve_and_print()

print("Размерность задачи 100х100:")
mat100, B100 = read_matrix("matrix100.txt")
solver = Solver(mat100, B100, eps=epsilon)
solver.solve_and_print()
