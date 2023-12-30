import KP_function1 as kp1
import KP_function2 as kp2
import numpy as np
import random
import KP_function3 as kp3
from inspect import signature

def Max_Min(razb, tx, bx, ty, by, tz, bz, function):
    dlt = (tx(2, 0) - bx(2, 0)) / razb
    x = np.arange(bx(2, 0), tx(2, 0) + dlt, dlt)
    if kp3.order()[0]=='z' or (kp3.order()[0]=='y' and kp3.order()[1]=='z'):
        y_max = np.max(ty(0, x))
        y_min = np.min(by(0, x))
    else:
        y_max = np.max(ty(x, 0))
        y_min = np.min(by(x, 0))
    x_max = tx(0, 0)
    x_min = bx(0, 0)
    X = np.linspace(x_min, x_max, razb)
    Y = np.linspace(y_min, y_max, razb)
    if kp3.order()[0]=='y' and kp3.order()[1]=='x':
        z_max1, z_min1 = Max_Min3(Y, X, function)
        z_max = max(z_max1, Max_Min2(Y, X, tz)[0])
        z_min = min(z_min1, Max_Min2(Y, X, bz)[1])
    else:
        z_max1, z_min1 = Max_Min3(X, Y, function)
        z_max=max(z_max1, Max_Min2(X, Y, tz)[0])
        z_min=min(z_min1, Max_Min2(X, Y, bz)[1])
    return x_max, x_min, y_max, y_min, z_max, z_min

def Border(X, Y, func):
    sig = signature(func)
    param = sig.parameters
    if len(param)==2:
        Z = np.zeros((len(X), len(Y)))
        for i in range(len(X)):
            for j in range(len(Y)):
                Z[i][j]=func(X[i][j], Y[i][j])
    else:
        Z = np.zeros((len(X), len(Y)))
        for i in range(len(X)):
            for j in range(len(Y)):
                Z[i][j] = func(X[i][j], Y[i][j], 0)
    return Z

def Max_Min2(X, Y, function):
    mmax = function(X[0], Y[0])
    mmin = function(X[0], Y[0])
    for i in range(len(X)):
        for j in range(len(Y)):
            if function(X[i], Y[j])>mmax:
                mmax=function(X[i], Y[j])
            if function(X[i], Y[j])<mmin:
                mmin=function(X[i], Y[j])
    return mmax, mmin

def Max_Min3(X, Y, function):
    mmax = function(X[0], Y[0], 0)
    mmin = function(X[0], Y[0], 0)
    for i in range(len(X)):
        for j in range(len(Y)):
            if function(X[i], Y[j], 0)>mmax:
                mmax=function(X[i], Y[j], 0)
            if function(X[i], Y[j], 0)<mmin:
                mmin=function(X[i], Y[j], 0)
    return mmax, mmin

def Monte_Karlo_1(N, razb):
    a=kp1.bx()
    b=kp1.tx()
    dlt=(b-a)/razb
    x=np.arange(a, b+dlt, dlt)
    mx = np.max(kp1.function(x))
    mn = np.min(kp1.function(x))
    print(mx, mn)
    x_max = b
    x_min = a
    y_max = mx
    y_min = 0
    if mn<0:
        y_min=mn
    cntP=0
    cntO=0
    Pnts=np.ones((3,N))
    for i in range(N):
        xi=x_min+(x_max-x_min)*random.random()
        yi=y_min+(y_max-y_min)*random.random()
        Pnts[0][i]=xi
        Pnts[1][i]=yi
        #print(xi)
        if kp1.function(xi)*yi>0:
            if kp1.function(xi)>0 and kp1.function(xi)>yi:
                Pnts[2][i]=0
                cntP+=1
            elif kp1.function(xi)<0 and kp1.function(xi)<yi:
                cntO+=1
                Pnts[2][i]=0
    dolP=cntP/N
    dolO=cntO/N
    S=(x_max-x_min)*(y_max-y_min)
    print(cntP, cntO)
    res=S*dolP-S*dolO
    return res, Pnts

def Monte_Karlo2(N, razb):
    if kp2.ty(2)==kp2.ty(1):
        dlt = (kp2.ty(2) - kp2.by(2)) / razb
        y = np.arange(kp2.by(2), kp2.ty(2) + dlt, dlt)
        x_max= np.max(kp2.tx(y))
        x_min= np.min(kp2.bx(y))
        y_max=kp2.ty(0)
        y_min=kp2.by(0)
    else:
        dlt = (kp2.tx(2) - kp2.bx(2)) / razb
        x = np.arange(kp2.bx(2), kp2.tx(2) + dlt, dlt)
        y_max = np.max(kp2.ty(x))
        y_min = np.min(kp2.by(x))
        x_max = kp2.tx(0)
        x_min = kp2.bx(0)

    cnt = 0
    Pnts=np.ones((3,N))
    print(x_max, x_min, y_max, y_min)
    for i in range(N):
        xi=x_min+(x_max-x_min)*random.random()
        yi=y_min+(y_max-y_min)*random.random()
        Pnts[0][i]=xi
        Pnts[1][i]=yi
        if kp2.ty(xi)>yi and kp2.by(xi)<yi:
            #print(kp2.ty(xi), yi, kp2.by(xi), yi)
            if kp2.tx(yi)>xi and kp2.bx(yi)<xi:
                cnt+=1
                Pnts[2][i]=0
    dol=cnt/N
    S=(x_max-x_min)*(y_max-y_min)
    print(cnt)
    res=S*dol
    return res, Pnts

def Monte_Karlo2_3D(N, razb):
    if kp2.ty(2)==kp2.ty(1):
        dlt = (kp2.ty(2) - kp2.by(2)) / razb
        y = np.arange(kp2.by(2), kp2.ty(2) + dlt, dlt)
        x_max= np.max(kp2.tx(y))
        x_min= np.min(kp2.bx(y))
        y_max=kp2.ty(0)
        y_min=kp2.by(0)
    else:
        dlt = (kp2.tx(2) - kp2.bx(2)) / razb
        x = np.arange(kp2.bx(2), kp2.tx(2) + dlt, dlt)
        y_max = np.max(kp2.ty(x))
        y_min = np.min(kp2.by(x))
        x_max = kp2.tx(0)
        x_min = kp2.bx(0)

    X = np.linspace(x_min, x_max, razb)
    Y = np.linspace(y_min, y_max, razb)
    z_max, z_min = Max_Min2(X, Y, kp2.function)
    cntP = 0
    cntO = 0
    Pnts=np.ones((4,N))
    print(x_max, x_min, y_max, y_min, z_max, z_min)
    for i in range(N):
        xi=x_min+(x_max-x_min)*random.random()
        yi=y_min+(y_max-y_min)*random.random()
        zi=z_min+(z_max-z_min)*random.random()
        Pnts[0][i]=xi
        Pnts[1][i]=yi
        Pnts[2][i]=zi
        if kp2.ty(xi)>yi and kp2.by(xi)<yi:
            if kp2.tx(yi)>xi and kp2.bx(yi)<xi:
                if kp2.function(xi, yi)>0 and kp2.function(xi, yi)>zi and zi>0:
                    cntP+=1
                    Pnts[3][i] = 0
                elif kp2.function(xi, yi)<0 and kp2.function(xi, yi)<zi and zi<0:
                    cntO += 1
                    Pnts[3][i] = 0
    dol1=cntP/N
    dol2=cntO/N
    S=(x_max-x_min)*(y_max-y_min)*(z_max-z_min)
    print(S)
    res=S*dol1-S*dol2
    return res, Pnts

def Monte_Karlo2_3D_3(N, razb):
    if kp3.order()[0]=='x':
        if kp3.order()[1]=='y':
            x_max, x_min, y_max, y_min, z_max, z_min = Max_Min(razb, kp3.tx, kp3.bx, kp3.ty, kp3.by, kp3.tz, kp3.bz, kp3.function)
        else:
            x_max, x_min, z_max, z_min, y_max, y_min = Max_Min(razb, kp3.tx, kp3.bx, kp3.tz, kp3.bz, kp3.ty, kp3.by, kp3.function)
    if kp3.order()[0]=='y':
        if kp3.order()[1]=='x':
            y_max, y_min, x_max, x_min, z_max, z_min = Max_Min(razb, kp3.ty, kp3.by, kp3.tx, kp3.bx, kp3.tz, kp3.bz, kp3.function)
        else:
            y_max, y_min, z_max, z_min, x_max, x_min = Max_Min(razb, kp3.ty, kp3.by, kp3.tz, kp3.bz, kp3.tx, kp3.bx, kp3.function)
    if kp3.order()[0]=='z':
        if kp3.order()[1]=='x':
            z_max, z_min, x_max, x_min, y_max, y_min = Max_Min(razb, kp3.tz, kp3.bz, kp3.tx, kp3.bx, kp3.ty, kp3.by, kp3.function)
        else:
            z_max, z_min, y_max, y_min, x_max, x_min = Max_Min(razb, kp3.tz, kp3.bz, kp3.ty, kp3.by, kp3.tx, kp3.bx, kp3.function)
    cntP = 0
    cntO = 0
    Pnts = np.ones((4, N))
    print(x_max, x_min, y_max, y_min, z_max, z_min)
    for i in range(N):
        xi=x_min+(x_max-x_min)*random.random()
        yi=y_min+(y_max-y_min)*random.random()
        zi=z_min+(z_max-z_min)*random.random()
        Pnts[0][i]=xi
        Pnts[1][i]=yi
        Pnts[2][i]=zi
        if kp3.ty(xi, zi)>yi and kp3.by(xi, zi)<yi:
            if kp3.tx(yi, zi)>xi and kp3.bx(yi, zi)<xi:
                if kp3.function(xi, yi, 0)==kp3.function(xi+1, yi+1, 0) and kp3.function(xi, yi, zi)==1:
                    if zi>0 and kp3.tz(xi, yi)>zi:
                        cntP+=1
                        Pnts[3][i]=0
                    elif zi<0 and kp3.bz(xi, yi)<zi:
                        cntO += 1
                        Pnts[3][i] = 0
    dol1=cntP/N
    dol2=cntO/N
    S=(x_max-x_min)*(y_max-y_min)*(z_max-z_min)
    print(S)
    res=S*dol1-S*dol2
    return res, Pnts

def Master_element_1(x0, x1, psi):
    return ((x1-x0)*(psi+1))/2+x0

def Gausse_1(n):
    if n==3:
        psi=[-np.sqrt(3/5), 0, np.sqrt(3/5)]
        w = [5/9, 8/9, 5/9]
    elif n==4:
        psi = [-np.sqrt(3/7 - 2/7*np.sqrt(6/5)), np.sqrt(3/7 - 2/7*np.sqrt(6/5)), -np.sqrt(3/7 + 2/7*np.sqrt(6/5)), np.sqrt(3/7 + 2/7*np.sqrt(6/5))]
        w = [(18+np.sqrt(30))/36, (18+np.sqrt(30))/36, (18-np.sqrt(30))/36, (18-np.sqrt(30))/36]
    res=0
    for i in range(len(psi)):
        res+=w[i]*(kp1.tx()-kp1.bx())/2*kp1.function(Master_element_1(kp1.bx(), kp1.tx(), psi[i]))
    return res



def Gausse_2(n):
    if n==4:
        psi = [-np.sqrt(1/3), -np.sqrt(1/3), np.sqrt(1/3), np.sqrt(1/3)]
        eta = [-np.sqrt(1/3), np.sqrt(1/3), -np.sqrt(1/3), np.sqrt(1/3)]
        w = [1, 1, 1, 1]
    elif n==12:
        psi = [-np.sqrt(6/3), np.sqrt(6/3), 0, 0, -np.sqrt((114-3*np.sqrt(583))/287), np.sqrt((114-3*np.sqrt(583))/287), -np.sqrt((114-3*np.sqrt(583))/287), np.sqrt((114-3*np.sqrt(583))/287), -np.sqrt((114+3*np.sqrt(583))/287), np.sqrt((114+3*np.sqrt(583))/287), -np.sqrt((114+3*np.sqrt(583))/287), np.sqrt((114+3*np.sqrt(583))/287)]
        eta = [0, 0, -np.sqrt(6/3), np.sqrt(6/3), -np.sqrt((114-3*np.sqrt(583))/287), -np.sqrt((114-3*np.sqrt(583))/287), np.sqrt((114-3*np.sqrt(583))/287), np.sqrt((114-3*np.sqrt(583))/287), -np.sqrt((114+3*np.sqrt(583))/287), -np.sqrt((114+3*np.sqrt(583))/287), np.sqrt((114+3*np.sqrt(583))/287), np.sqrt((114+3*np.sqrt(583))/287)]
        w = [98/405, 98/405, 98/405, 98/405, 307/810+923/(270*np.sqrt(583)), 307/810+923/(270*np.sqrt(583)), 307/810+923/(270*np.sqrt(583)), 307/810+923/(270*np.sqrt(583)), 307/810-923/(270*np.sqrt(583)), 307/810-923/(270*np.sqrt(583)), 307/810-923/(270*np.sqrt(583)), 307/810-923/(270*np.sqrt(583))]
    res=0
    for i in range(len(psi)):
        res += w[i] * (kp2.tx(0)-kp2.bx(0))*(kp2.ty(0)-kp2.by(0))/4 * kp2.function(Master_element_1(kp2.bx(0), kp2.tx(0), psi[i]), Master_element_1(kp2.by(0), kp2.ty(0), eta[i]))
        #print(Master_element_1(kp2.bx(0), kp2.tx(0), psi[i]))
        #print(Master_element_1(kp2.by(0), kp2.ty(0), psi[i]))
    return res

def Gausse_3(n):
    if n==8:
        psi = [-np.sqrt(1 / 3), -np.sqrt(1 / 3), -np.sqrt(1 / 3), -np.sqrt(1 / 3), np.sqrt(1 / 3), np.sqrt(1 / 3), np.sqrt(1 / 3), np.sqrt(1 / 3)]
        eta = [-np.sqrt(1 / 3), -np.sqrt(1 / 3), np.sqrt(1 / 3), np.sqrt(1 / 3), -np.sqrt(1 / 3), -np.sqrt(1 / 3), np.sqrt(1 / 3), np.sqrt(1 / 3)]
        ksi = [-np.sqrt(1 / 3), np.sqrt(1 / 3), -np.sqrt(1 / 3), np.sqrt(1 / 3), -np.sqrt(1 / 3), np.sqrt(1 / 3), -np.sqrt(1 / 3), np.sqrt(1 / 3)]
        w = [1, 1, 1, 1, 1, 1, 1, 1]
    elif n==14:
        psi = [-np.sqrt(19/30), np.sqrt(19/30), 0, 0, 0, 0, -np.sqrt(19/33), -np.sqrt(19/33), -np.sqrt(19/33), -np.sqrt(19/33), np.sqrt(19/33), np.sqrt(19/33), np.sqrt(19/33), np.sqrt(19/33)]
        eta = [0, 0, -np.sqrt(19/30), np.sqrt(19/30), 0, 0, -np.sqrt(19/33), -np.sqrt(19/33), np.sqrt(19/33), np.sqrt(19/33), -np.sqrt(19/33), -np.sqrt(19/33), np.sqrt(19/33), np.sqrt(19/33)]
        ksi = [0, 0, 0, 0, -np.sqrt(19/30), np.sqrt(19/30), -np.sqrt(19/33), np.sqrt(19/33), -np.sqrt(19/33), np.sqrt(19/33), -np.sqrt(19/33), np.sqrt(19/33), -np.sqrt(19/33), np.sqrt(19/33)]
        w = [320/361, 320/361, 320/361, 320/361, 320/361, 320/361, 121/361, 121/361, 121/361, 121/361, 121/361, 121/361, 121/361, 121/361]
    res = 0
    for i in range(len(psi)):
        res += w[i] * (kp3.tx(0, 0) - kp3.bx(0, 0)) * (kp3.ty(0, 0) - kp3.by(0, 0))*(kp3.tz(0,0)-kp3.bz(0,0))/ 8 * kp3.function(Master_element_1(kp3.bx(0, 0), kp3.tx(0, 0), psi[i]), Master_element_1(kp3.by(0, 0), kp3.ty(0, 0), eta[i]), Master_element_1(kp3.bz(0,0), kp3.tz(0,0), ksi[i]))
    return res

#print(Gausse_2(1000))
#print(Monte_Karlo2(5000, 1000))
#print(Monte_Karlo_1(50000000, 100))