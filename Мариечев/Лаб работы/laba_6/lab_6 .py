import sympy as sp
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline 

N=1000 #  x
K=100 #  t
T=0.001

L=0
R=math.pi
# �����
h=(R-L)/N 
tau=(T-0)/K

# f1=1  ; f1=2
# f2=1  ; f2=2

print('\ntau <= h')
print(tau,'<=', h)

def phi1(t):
    return math.exp(-t)

def phi2(t):
    return -math.exp(-t)

def psi1(x):
    return math.cos(x)

def psi2(x):
    return -math.cos(x)

def f(x,t):
    return math.exp(-t)*math.sin(x)

def dpsi1_dx(x):
    return -math.sin(x)
def d2psi1_dx2(x):
    return -math.cos(x)


def draw(u,name):

    TT=[j*tau for j in range(K+1)]
    H=[L+j*h for j in range(N+1)]

    plt.figure(name)
    plt.subplot(1,2,1)
    x = np.linspace(0,R, 1000)
    plt.plot(x, np.exp(0)*(np.cos(x)),label='U(x,t)=exp(-t)*cos(x)  t=0')
    plt.plot(x, np.exp(round(K/2)*tau)*(np.cos(x)),label='U(x,t)=exp(-t)*cos(x)  t=T/2')
    plt.plot(x, np.exp(-T)*(np.cos(x)),label='U(x,t)=exp(-t)*cos(x)  t=T')
    plt.title('u(x), t=const')
    #plt.plot(H,u[0],'ro')
    xnew = np.linspace(0, R, 1000) 
    gfg = make_interp_spline(H, u[0], k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='t=0')   
    gfg = make_interp_spline(H, u[round(K/2)], k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='t=T/2')  
    gfg = make_interp_spline(H, u[K], k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='t=T')  
    plt.legend(fontsize=7, loc=0)

    plt.subplot(1,2,2)
    t = np.linspace(0,T, 300)
    plt.plot(t, np.exp(-t)*(np.cos(0)),label='exp(-a*t)*cos(x+b*t)  t=0')
    plt.plot(t, np.exp(-t)*(np.cos(round(N/2)*h)),label='exp(-a*t)*cos(x+b*t)  t=R/2')
    plt.plot(t, np.exp(-t)*(np.cos(R)),label='exp(-a*t)*cos(x+b*t)  t=R')
    plt.title('u(t), x=const')
    #plt.plot(H,u[0],'ro')
    xnew = np.linspace(0, T, 300) 

    ut0=[u[j][0] for j in range(K+1)]
    ut1=[u[j][round(N/2)] for j in range(K+1)]
    ut2=[u[j][N] for j in range(K+1)]

    gfg = make_interp_spline(TT, ut0, k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='x=0')   
    gfg = make_interp_spline(TT, ut1, k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='x=R/2')  
    gfg = make_interp_spline(TT, ut2, k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='x=R')  
    plt.legend(fontsize=7, loc=0)
    



def running_method(A,B,n):
    P=np.zeros(n)
    Q=np.zeros(n)

    P[0]=-A[0][1]/A[0][0]
    Q[0]=B[0]/A[0][0]

    for i in range(1,n-1):
        P[i]=-A[i][i+1]/(A[i][i]+A[i][i-1]*P[i-1])
        Q[i]=(B[i]-A[i][i-1]*Q[i-1])/(A[i][i]+A[i][i-1]*P[i-1])

    P[n-1]=0
    Q[n-1]=(B[n-1]-A[n-1][n-2]*Q[n-2])/(A[n-1][n-1]+A[n-1][n-2]*P[n-2])
    x=np.zeros(n)
    x[n-1]=Q[n-1]

    for i in range(n-2,-1,-1):
        x[i]=P[i]*x[i+1]+Q[i]
    return(x)

def solve(f1,f2):
    u=np.zeros((K+1,N+1))


    k=0
    for i in range(N+1):
        x=L+i*h
        u[k][i]=psi1(x)


    k=1
    for i in range(N+1):
        if f1==1:
            x=L+i*h
            u[k][i]=psi1(x)+tau*psi2(x)
        if f1==2:
            x=L+i*h
            t=k*tau
            #u[k][i]=psi1(x)+tau*psi2(x)+tau**2*(-math.cos(x)-math.sin(x)-psi1(x)+f(x,0)-psi2(x))/2
            if i==0:
                u[k][i]=phi1(t)
            if i==N:
               u[k][i]=phi2(t)
            else :
               u[k][i]=psi1(x)+tau*psi2(x)+tau**2*(dpsi1_dx(x)+d2psi1_dx2(x)-psi1(x)+f(x,0)-3*psi2(x))/2


    #  k+1

    for k in range(1,K):
        if f2==1:
            for i in range(1,N):    
                u[k+1][i]=((u[k][i+1]-2*u[k][i]+u[k][i-1])/h**2+(u[k][i+1]-u[k][i-1])/(2*h)-u[k][i]+f(L+i*h,k*tau)+(2*u[k][i]-u[k-1][i])/tau**2+3*u[k-1][i]/(2*tau))*((2*tau**2)/(2+3*tau))
    
            t=tau*(k+1)
            u[k+1][0]=phi1(t) 
            u[k+1][N]=phi2(t) 
        if f2==2:
            A=np.zeros((N+1,N+1))
            B=np.zeros((N+1,1))
            t=tau*(k+1)

            A[0][0]=1 
            B[0][0]=phi1(t)

            A[N][N]=1 
            B[N][0]=phi2(t)

            for i in range(1,N):
                x=L+i*h
                A[i][i-1]=-1/h**2+1/(2*h)
                A[i][i]=1/tau**2+3/(2*tau)+2/h**2+1
                A[i][i+1]=-1/h**2-1/(2*h)
                B[i][0]=(2*u[k][i]-u[k-1][i])/tau**2+3*u[k-1][i]/(2*tau)+f(x,t)
            X=running_method(A,B,N+1)

            for i in range(0,N+1):
                u[k+1][i]=X[i]
    return(u)

 

U=np.zeros((K+1,N+1))   
for k in range(0,K+1):
    for i in range(0, N+1):
        x=L+i*h 
        t=k*tau 
        U[k][i]=math.exp(-t)*math.cos(x)

def error(u):
    error=0 
    for k in range(0,K+1):
        for i in range(0, N+1):
            error+=(U[k][i]-u[k][i])**2
    print(math.sqrt(error/N/K))


u11=solve(1,1)
draw(u11,'explicit method 1')
print('\neroor explicit method 1')
error(u11)

u21=solve(2,1)
draw(u21,'explicit method 2')
print('\neroor explicit method 2')
error(u21)

u12=solve(1,2)
draw(u12,'implicit_method 1')
print('\neroor implicit_method 1')
error(u12)

u22=solve(2,2)
draw(u22,'implicit_method 2')
print('\neroor implicit_method 2')
error(u22)




plt.show() 
