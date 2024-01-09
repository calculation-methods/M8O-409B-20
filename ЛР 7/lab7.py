import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
import matplotlib.cm

def stop(L,U,x,y):
    maxx = 0
    for i in range(len(x)):
        for j in range(len(y)):
            if abs(U[i,j] - L[i,j]) > maxx:
                maxx = abs(U[i,j] - L[i,j])
    return maxx

def true_fval(x, y):
    return x*np.cos(y)


lbx = 0
ubx = 1
nx = 40
hx = (ubx - lbx)/nx


lby = 0
uby = np.pi/2
ny = 40
hy = (uby - lby)/ny




x = np.arange(lbx, ubx + hx, hx)
y = np.arange(lby, uby + hy, hy)


eps = 0.000001


def MPI(hx, hy, eps, lbx, lby, ubx, uby):
    x = np.arange(lbx, ubx + hx, hx)
    y = np.arange(lby, uby + hy, hy)
    U = np.zeros((len(x),len(y)))
    for i in range(len(x)):
        U[i,0] = x[i]
        U[i,-1] = 0    
    k = 0
    for j in range(1,len(y)-1):
        for i in range(len(x)):     
            U[i,j] = U[i,0] * (y[-1]-y[j])/y[-1]        
    while True:
        k = k+1
        L = np.copy(U)
        U = np.zeros((len(x),len(y)))
        for i in range(len(x)):
            U[i,0] = x[i]
            U[i,-1] = 0  
        for j in range(1, len(y)-1):
            U[0,j] = L[1,j]-hx*np.cos(y[j])
            U[-1,j] = L[-2,j]/(1-hx)

        for j in range(1, len(y)- 1):
            for i in range(1, len(x)- 1):
                U[i,j] = (L[i+1,j] + L[i-1,j] + L[i,j+1] + L[i,j-1])/(4-hx*hy)
        if stop(L,U,x,y) <= eps:
            print('MPI eps = ',stop(L,U,x,y))
            print('MPI k = ', k)
            break

    return U

def ZEI(hx, hy, eps, lbx, lby, ubx, uby):
    x = np.arange(lbx, ubx + hx, hx)
    y = np.arange(lby, uby + hy, hy)
    U = np.zeros((len(x),len(y)))
    for i in range(len(x)):
        U[i,0] = x[i]
        U[i,-1] = 0    
    k = 0
    for j in range(1,len(y)-1):
        for i in range(len(x)):     
            U[i,j] = U[i,0] * (y[-1]-y[j])/y[-1]        
    while True:
        k = k+1
        L = np.copy(U)
        U = np.zeros((len(x),len(y)))
        for i in range(len(x)):
            U[i,0] = x[i]
            U[i,-1] = 0  
        for j in range(1, len(y)-1):
            U[0,j] = L[1,j]-hx*np.cos(y[j])
            U[-1,j] = L[-2,j]/(1-hx)

        for j in range(1, len(y)- 1):
            for i in range(1, len(x)- 1):
                U[i,j] = (L[i+1,j] + U[i-1,j] + L[i,j+1] + U[i,j-1])/(4-hx*hy)
        if stop(L,U,x,y) <= eps:
            print('ZEI eps = ',stop(L,U,x,y))
            print('ZEI k = ', k)
            break

    return U

def MPI_relax(hx, hy, eps, lbx, lby, ubx, uby, C):
    x = np.arange(lbx, ubx + hx, hx)
    y = np.arange(lby, uby + hy, hy)
    U = np.zeros((len(x),len(y)))
    for i in range(len(x)):
        U[i,0] = x[i]
        U[i,-1] = 0    
    k = 0
    for j in range(1,len(y)-1):
        for i in range(len(x)):     
            U[i,j] = U[i,0] * (y[-1]-y[j])/y[-1]        
    while True:
        k = k+1
        L = np.copy(U)
        U = np.zeros((len(x),len(y)))
        for i in range(len(x)):
            U[i,0] = x[i]
            U[i,-1] = 0  
        for j in range(1, len(y)-1):
            U[0,j] = L[1,j]-hx*np.cos(y[j])
            U[-1,j] = L[-2,j]/(1-hx)

        for j in range(1, len(y)- 1):
            for i in range(1, len(x)- 1):
                U[i,j] = ((L[i+1,j] + L[i-1,j] + L[i,j+1] + L[i,j-1])/(4-hx*hy))*C + L[i,j]*(1-C)
        if stop(L,U,x,y) <= eps:
            print('MPI_relax eps = ',stop(L,U,x,y))
            print('MPI_relax k = ', k)
            break

    return U


def ZEI_relax(hx, hy, eps, lbx, lby, ubx, uby,C):
    x = np.arange(lbx, ubx + hx, hx)
    y = np.arange(lby, uby + hy, hy)
    U = np.zeros((len(x),len(y)))
    for i in range(len(x)):
        U[i,0] = x[i]
        U[i,-1] = 0    
    k = 0
    for j in range(1,len(y)-1):
        for i in range(len(x)):     
            U[i,j] = U[i,0] * (y[-1]-y[j])/y[-1]        
    while True:
        k = k+1
        L = np.copy(U)
        U = np.zeros((len(x),len(y)))
        for i in range(len(x)):
            U[i,0] = x[i]
            U[i,-1] = 0  
        for j in range(1, len(y)-1):
            U[0,j] = L[1,j]-hx*np.cos(y[j])
            U[-1,j] = L[-2,j]/(1-hx)

        for j in range(1, len(y)- 1):
            for i in range(1, len(x)- 1):
                U[i,j] = ((L[i+1,j] + U[i-1,j] + L[i,j+1] + U[i,j-1])/(4-hx*hy))*C + L[i,j]*(1-C)
        if stop(L,U,x,y) <= eps:
            print('ZEI_relax eps = ',stop(L,U,x,y))
            print('ZEI_relax k = ', k)
            break

    return U


def plot_U(z, u):
    x = np.arange(lbx, ubx + hx, hx)
    y = np.arange(lby, uby + hy, hy)
    xx, yy = np.meshgrid(x, y)
    fig = plt.figure(figsize = (20, 7))
    d = abs(z - u)
        
    
    ax1 = fig.add_subplot(1, 2, 1)
    
    plt.xlabel("x")
    plt.ylabel("error")
    
    ax2 = fig.add_subplot(1, 2, 2, projection = "3d")

    plt.xlabel("x")
    plt.ylabel("y")
    
    surf = ax2.plot_surface(xx, yy, d,
                             edgecolors = ["black"], linewidth = 1,
                             cmap = matplotlib.cm.Spectral, shade = True, antialiased = True)
    fig.colorbar(surf)

    ppp = []
    
    for i in range(len(y)):
        mmm = d[i][0]
        for j in range(len(x)):
            mmm = max(d[i][j], mmm)
        ppp.append(mmm)

    ax1.plot(y, ppp)  

    fig1, ax = plt.subplots(1, 1, figsize = (20, 20), subplot_kw = {"projection": "3d"})
    ax.view_init(30, 30)
    surf = ax.plot_surface(xx, yy, z,
                         edgecolors = ["black"], linewidth = 1,
                         cmap = matplotlib.cm.Spectral, shade = True, antialiased = True)
    surf1 = ax.plot_surface(xx, yy, u,
                         edgecolors = ["red"], linewidth = 1,
                         cmap = matplotlib.cm.Spectral, shade = False, antialiased = True)
    plt.xlabel("x")
    plt.ylabel("y")

    #fig1.colorbar(surf)
    #fig1.colorbar(surf1)
    
    plt.show()
    return


z = np.zeros((len(x),len(y)))
for i in range(len(x)):
        for j in range(len(y)):
            z[i,j]= true_fval(x[i],y[j])

#u1 = MPI(hx, hy, eps, lbx, lby, ubx, uby)
#u2 = ZEI(hx, hy, eps, lbx, lby, ubx, uby)
#u3 = MPI_relax(hx, hy, eps, lbx, lby, ubx, uby,0.5)
u4 = ZEI_relax(hx, hy, eps, lbx, lby, ubx, uby,1.91)

#plot_U(z,u1)
#plot_U(z,u2)
#plot_U(z,u3)
plot_U(z,u4)

