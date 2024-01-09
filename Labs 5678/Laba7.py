import math

import matplotlib.pyplot as plt
import numpy as np

a = 0
b = 1.005

c = 0
d = 1.005

steps_x = 0.1
steps_y = 0.1

eps = 0.0000001


def ux_0(y):
    return y


def ux_1(y):
    return 1 + y


def uy_0(x):
    return x


def uy_1(x):
    return 1 + x


def analitic(x, y):
    return x + y


print("______________________________________________________________")
print('\td^2u/dx^2 + d^2u/dy^2 = 0, ',
      'u(0, y) = y,',
      'u(1, y) = 1 + y,',
      'u(x ,0) = x,',
      'u(x, 1) = 1 + x,',
      'Аналитическое решение: U(x, y) = x + y', sep='\n\t')
print("______________________________________________________________")


def build_analitic_solution(x_range, y_range, steps_x, steps_y):
    x_values = np.arange(x_range[0], x_range[1], steps_x)
    y_values = np.arange(y_range[0], y_range[1], steps_y)
    res = np.zeros((len(x_values), len(y_values)))

    for i in range(len(x_values)):
        for j in range(len(y_values)):
            res[i, j] = analitic(x_values[i], y_values[j])
    return res


def error(arr1, arr2):
    return np.max(np.abs(arr1 - arr2))


def normal(vector):
    norm = 0
    for i in vector:
        norm += i ** 2
    return math.sqrt(norm)


# Конечно-разностная схема
def Shema(x_range, y_range, steps_x, steps_y, method, phi_0=ux_0, phi_1=ux_1, uy_0=uy_0, uy_1=uy_1):
    x = np.arange(x_range[0], x_range[1], steps_x)
    y = np.arange(y_range[0], y_range[1], steps_y)

    # Инициализируем сетку с граничными условиями.
    res = np.zeros((len(x), len(y)))

    for cur_x_id in range(len(x)):
        res[cur_x_id][0] = uy_0(x[cur_x_id])
        res[cur_x_id][-1] = uy_1(x[cur_x_id])

    for cur_y_id in range(len(y)):
        res[0][cur_y_id] = phi_0(y[cur_y_id])
        res[-1][cur_y_id] = phi_1(y[cur_y_id])

    # Система уровнений
    mapping = np.zeros((len(x), len(y)), dtype='int')
    cur_eq_id = 0
    for cur_x_id in range(1, len(x) - 1):
        for cur_y_id in range(1, len(y) - 1):
            mapping[cur_x_id][cur_y_id] = cur_eq_id
            cur_eq_id += 1

    nums_of_equations = (len(x) - 2) * (len(y) - 2)
    A = np.zeros((nums_of_equations, nums_of_equations))
    b = np.zeros((nums_of_equations))
    for cur_x_id in range(1, len(x) - 1):
        for cur_y_id in range(1, len(y) - 1):
            cur_eq_id = mapping[cur_x_id][cur_y_id]

            A[cur_eq_id][mapping[cur_x_id][cur_y_id]] = 1  # u_{i, j}

            if cur_y_id - 1 == 0:
                b[cur_eq_id] += uy_0(x[cur_x_id]) * steps_x ** 2 / (2 * (steps_x ** 2 + steps_y ** 2))
            else:
                A[cur_eq_id][mapping[cur_x_id][cur_y_id - 1]] = -steps_x ** 2 / (2 * (steps_x ** 2 + steps_y ** 2))

            if cur_y_id + 1 == len(y) - 1:
                b[cur_eq_id] += uy_1(x[cur_x_id]) * steps_x ** 2 / (2 * (steps_x ** 2 + steps_y ** 2))
            else:
                A[cur_eq_id][mapping[cur_x_id][cur_y_id + 1]] = -steps_x ** 2 / (
                            2 * (steps_x ** 2 + steps_y ** 2))  # u_{i, j+1}

            if cur_x_id - 1 == 0:
                b[cur_eq_id] += phi_0(y[cur_y_id]) * steps_y ** 2 / (2 * (steps_x ** 2 + steps_y ** 2))
            else:
                A[cur_eq_id][mapping[cur_x_id - 1][cur_y_id]] = -steps_y ** 2 / (
                            2 * (steps_x ** 2 + steps_y ** 2))  # u_{i-1, j}

            if cur_x_id + 1 == len(x) - 1:
                b[cur_eq_id] += phi_1(y[cur_y_id]) * steps_y ** 2 / (2 * (steps_x ** 2 + steps_y ** 2))
            else:
                A[cur_eq_id][mapping[cur_x_id + 1][cur_y_id]] = -steps_y ** 2 / (
                            2 * (steps_x ** 2 + steps_y ** 2))  # u_{i+1, j}

    # СЛАУ
    ans, iters = method(A, b, eps)
    for cur_x_id in range(1, len(x) - 1):
        for cur_y_id in range(1, len(y) - 1):
            res[cur_x_id][cur_y_id] = ans[mapping[cur_x_id][cur_y_id]]

    return res, iters


analytical_solution = build_analitic_solution(x_range=(a, b), y_range=(c, d), steps_x=steps_x, steps_y=steps_y, )
solutions = dict()
solutions["analytical solution"] = analytical_solution


# Метод простых итераций для решения СЛАУ
def iterative(A, b, eps):
    n = len(b)
    # Ax=b -> x = alpha * x + beta
    alpha = np.zeros_like(A)
    beta = np.zeros_like(b)

    for i in range(n):
        alpha[i, :] = -A[i, :] / A[i, i]
        alpha[i, i] = 0
        beta[i] = b[i] / A[i, i]

    # счет итераций
    k = 0
    cur_x = np.copy(beta)
    converge = False
    while not converge:
        prev_x = np.copy(cur_x)
        cur_x = np.zeros_like(prev_x)
        for i in range(alpha.shape[0]):
            cur_x[i] = beta[i]
            for j in range(alpha.shape[1]):
                cur_x[i] += alpha[i][j] * prev_x[j]
        k += 1
        converge = normal(prev_x - cur_x) <= eps
    return cur_x, k


iterative_solution, iterative_iters = Shema(x_range=(a, b), y_range=(c, d), steps_x=steps_x, steps_y=steps_y,
                                            method=iterative)

solutions["iterative solution"] = iterative_solution

print(f'max abs error iterative = {error(iterative_solution, analytical_solution)}')
print(f'iter iterative = {iterative_iters}')


def multiplication(alpha, x, beta):
    res = np.copy(x)
    for i in range(alpha.shape[0]):
        res[i] = beta[i]
        for j in range(alpha.shape[1]):
            res[i] += alpha[i][j] * res[j]
    return res


# Перемножений матриц
def seidel_multiplication(alpha, x, beta):
    res = np.copy(x)
    for i in range(alpha.shape[0]):
        res[i] = beta[i]
        for j in range(alpha.shape[1]):
            res[i] += alpha[i][j] * res[j]
    return res


def Zeidel(A, b, eps):
    n = len(b)
    # Ax=b -> x = alpha * x + beta
    alpha = np.zeros_like(A)
    beta = np.zeros_like(b)

    for i in range(n):
        alpha[i, :] = -A[i, :] / A[i, i]
        alpha[i, i] = 0
        beta[i] = b[i] / A[i, i]

    # Счет итераций
    k = 0
    cur_x = np.copy(beta)
    converge = False
    while not converge:
        prev_x = np.copy(cur_x)
        cur_x = multiplication(alpha, prev_x, beta)
        k += 1
        converge = normal(prev_x - cur_x) <= eps
    return cur_x, k


seidel_solution, seidel_iters = Shema(x_range=(a, b), y_range=(c, d), steps_x=steps_x, steps_y=steps_y, method=Zeidel)

solutions["Zeidel solution"] = seidel_solution
print(f'max abs error Zeidel= {error(seidel_solution, analytical_solution)}')
print(f'Zeidel iter = {seidel_iters}')


# релаксация для решения СЛАУ
def relaxation(A, b, eps, w=1.5):
    n = len(b)
    # Ax=b -> x = alpha * x + beta
    alpha = np.zeros_like(A)
    beta = np.zeros_like(b)

    for i in range(n):
        alpha[i, :] = -A[i, :] / A[i, i]
        alpha[i, i] = 0
        beta[i] = b[i] / A[i, i]

    k = 0
    cur_x = np.copy(beta)
    converge = False
    while not converge:
        prev_x = np.copy(cur_x)
        cur_x = multiplication(alpha, prev_x, beta)
        cur_x = w * cur_x + (1 - w) * prev_x
        k += 1
        converge = normal(prev_x - cur_x) <= eps
    return cur_x, k


relaxation_solution, relaxation_iters = Shema(x_range=(a, b), y_range=(c, d), steps_x=steps_x, steps_y=steps_y,
                                              method=relaxation)
solutions["relaxation solution"] = relaxation_solution
print(f'max abs error relaxation = {error(relaxation_solution, analytical_solution)}')
print(f'iter relaxation = {relaxation_iters}')


def plot_results_and_errors(solutions, cur_y, x_range, y_range, steps_x, steps_y, analytical_solution_name):
    x_values = np.arange(*x_range, steps_x)
    y_values = np.arange(*y_range, steps_y)
    cur_y_id = np.argmin(np.abs(y_values - cur_y))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    colors = iter(['blue', 'red', 'green'])
    ax1.plot(x_values, solutions[analytical_solution_name][:, cur_y_id], label=analytical_solution_name, color='blue')
    for method_name, solution in solutions.items():
        if method_name != analytical_solution_name:
            ax1.plot(x_values, solution[:, cur_y_id], label=method_name, color=next(colors))

    ax1.set_title('Solution Profiles')
    ax1.set_xlabel('x')
    ax1.set_ylabel('Solution Value')
    ax1.legend()
    ax1.grid()

    colors = iter(['red', 'green', 'purple'])
    for method_name, solution in solutions.items():
        if method_name != analytical_solution_name:
            max_abs_errors = np.array([
                error(solution[:, i], solutions[analytical_solution_name][:, i])
                for i in range(len(y_values))
            ])
            ax2.plot(y_values, max_abs_errors, label=method_name, color=next(colors))

    ax2.set_title('Max Absolute Errors')
    ax2.set_xlabel('y')
    ax2.set_ylabel('Max Abs Error')
    ax2.legend()
    ax2.grid()

    plt.tight_layout()
    plt.show()


plot_results_and_errors(solutions=solutions, cur_y=0.5, x_range=(a, b), y_range=(c, d), steps_x=steps_x,
                        steps_y=steps_y, analytical_solution_name="analytical solution")
