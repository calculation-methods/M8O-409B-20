import numpy as np

def Tridiagonal_matrix_algorithm(L, y):
    a = [0]
    b =[L[0][0]]
    c = [L[0][1]]
    n = len(y)

    x = []
    alpha = []
    beta = []

    for i in range(1, n - 1):
        a.append(L[i][i-1])
        b.append(L[i][i])
        c.append(L[i][i+1])

    a.append(L[n-1][n-2])
    b.append(L[n-1][n-1])
    c.append(0)

    a = np.array(a)
    b = np.array(b)
    c = np.array(c)

    alpha.append(-c[0]/b[0])
    beta.append(y[0]/b[0])

    for i in range(1, n-1):
        d = a[i] * alpha[i-1] + b[i]
        alpha.append(-c[i] / d)
        beta.append((y[i] - a[i] * beta[i-1]) / d)
    alpha.append(0)
    beta.append((y[n - 1] - a[n-1] * beta[n-2]) / (a[n-1] * alpha[n-2] + b[n-1]))

    x.append(beta[n-1])

    for i in range(1, n):
        x.append(x[i - 1] * alpha[n - 1 - i] + beta[n - 1 - i])

    return x[::-1]