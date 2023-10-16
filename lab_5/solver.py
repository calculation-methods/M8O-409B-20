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

    print(a, b, c)

    alpha.append(-b[0]/c[0])
    beta.append(y[0]/c[0])

    for i in range(0, n-1):
        alpha.append(-c[i] / (a[i] * alpha[i] + b[i]))
        beta.append((y[i] - a[i] * beta[i]) / (a[i] * alpha[i] + b[i]))

    x.append((y[n-1] - a[n-1] * beta[n-1]) / (a[n-1] * alpha[n-1] + b[n-1]))

    for i in range(0, n - 1):
        x.append(x[i] * alpha[n - i - 1] + beta[n - i - 1])

    return x[::-1]
