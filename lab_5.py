import sympy as sp
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline 

alpha=1
betta=-1
gamma=1
delta=-1

n=10 # разбиения x 30
k=100# разбиения t 1500
T=0.1
a=1
b=1
l=0
R=math.pi
# сетка
h=(R-l)/n 
tau=(T-0)/k

print(tau,'<=', h**2/2*a)

#f1=-math.exp(-a*t)*(math.cos(b*t)+math.sin(b*t))
#f2=math.exp(-a*t)*(math.cos(b*t)+math.sin(b*t))
#f3=math.cos(x)

def draw(u,name):

    TT=[j*tau for j in range(k+1)]
    H=[l+j*h for j in range(n+1)]

    plt.figure(name)
    plt.subplot(1,2,1)
    x = np.linspace(0,R, 1000)
    plt.plot(x, np.exp(-a*0)*(np.cos(x+b*0)),label='exp(-a*t)*cos(x+b*t)  t=0')
    plt.plot(x, np.exp(-a*T/2)*(np.cos(x+b*T/2)),label='exp(-a*t)*cos(x+b*t)  t=T/2')
    plt.plot(x, np.exp(-a*T)*(np.cos(x+b*T)),label='exp(-a*t)*cos(x+b*t)  t=T')
    plt.title('u(x), t=const')
    #plt.plot(H,u[0],'ro')
    xnew = np.linspace(0, R, 1000) 
    gfg = make_interp_spline(H, u[0], k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='t=0')   
    gfg = make_interp_spline(H, u[round(k/2)], k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='t=T/2')  
    gfg = make_interp_spline(H, u[k], k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='t=T')  
    plt.legend(fontsize=7, loc=0)

    plt.subplot(1,2,2)
    t = np.linspace(0,T, 300)
    plt.plot(t, np.exp(-a*t)*(np.cos(0+b*t)),label='exp(-a*t)*cos(x+b*t)  t=0')
    plt.plot(t, np.exp(-a*t)*(np.cos(R/2+b*t)),label='exp(-a*t)*cos(x+b*t)  t=R/2')
    plt.plot(t, np.exp(-a*t)*(np.cos(R+b*t)),label='exp(-a*t)*cos(x+b*t)  t=R')
    plt.title('u(t), x=const')
    #plt.plot(H,u[0],'ro')
    xnew = np.linspace(0, T, 300) 

    ut0=[u[j][0] for j in range(k+1)]
    ut1=[u[j][round(n/2)] for j in range(k+1)]
    ut2=[u[j][n] for j in range(k+1)]

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

# явный метод 
# 1 способ апроксимации


# t=0
K=0
u=np.zeros((k+1,n+1))

for i in range(n+1):
    x=l+i*h
    u[0][i]=math.cos(x)

for K in range(0,k):
    for i in range(1,n):    
        u[K+1][i]=tau*(a*(u[K][i+1]-2*u[K][i]+u[K][i-1])/h**2+b*(u[K][i+1]-u[K][i-1])/(2*h))+u[K][i]
    t=tau*(K+1)
    u[K+1][0]=(-alpha/h)*u[K+1][1]/(betta-alpha/h)-math.exp(-a*t)*(math.cos(b*t)+math.sin(b*t))/(betta-alpha/h)
    u[K+1][n]=(gamma/h)*u[K+1][n-1]/(delta+gamma/h)+math.exp(-a*t)*(math.cos(b*t)+math.sin(b*t))/(delta+gamma/h)

print('eroor explicit method 1')
xx=0
time=0
error1=math.exp(-a*time*tau)*math.cos(xx*h+b*time*tau)-u[time][xx]
xx=round(n/2)
time=round(k/2)
error2=math.exp(-a*time*tau)*math.cos(xx*h+b*time*tau)-u[time][xx]
xx=n
time=k
error3=math.exp(-a*time*tau)*math.cos(xx*h+b*time*tau)-u[time][xx]

print(error1,error2,error3)


draw(u,'explicit method 1')



# 2 способ апроксимации

K=0
u=np.zeros((k+1,n+1))

for i in range(n+1):
    x=l+i*h
    u[0][i]=math.cos(x)

for K in range(0,k):
    for i in range(1,n):    
        u[K+1][i]=tau*(a*(u[K][i+1]-2*u[K][i]+u[K][i-1])/h**2+b*(u[K][i+1]-u[K][i-1])/(2*h))+u[K][i]
    t=tau*(K+1)
    u[K+1][0]=(-alpha/(2*h))*(4*u[K+1][1]-u[K+1][2])/(betta-3*alpha/(2*h))-math.exp(-a*t)*(math.cos(b*t)+math.sin(b*t))/(betta-3*alpha/(2*h))
    u[K+1][n]=(-gamma/(2*h))*(u[K+1][n-2]-4*u[K+1][n-1])/(delta+3*gamma/(2*h))+math.exp(-a*t)*(math.cos(b*t)+math.sin(b*t))/(delta+3*gamma/(2*h))

print('eroor explicit method 2')
xx=0
time=0
error1=math.exp(-a*time*tau)*math.cos(xx*h+b*time*tau)-u[time][xx]
xx=round(n/2)
time=round(k/2)
error2=math.exp(-a*time*tau)*math.cos(xx*h+b*time*tau)-u[time][xx]
xx=n
time=k
error3=math.exp(-a*time*tau)*math.cos(xx*h+b*time*tau)-u[time][xx]

print(error1,error2,error3)

draw(u,'explicit method 2')





# 3 способ апроксимации

K=0
u=np.zeros((k+1,n+1))

for i in range(n+1):
    x=l+i*h
    u[0][i]=math.cos(x)

D=2*a*tau*h*betta-tau*b*betta*h**2-2*a*alpha*tau-b*h**2
C=2*a*gamma*tau+gamma*h**2+2*a*h*tau*delta+b*tau*delta*h**2

for K in range(0,k):
    for i in range(1,n):    
        u[K+1][i]=tau*(a*(u[K][i+1]-2*u[K][i]+u[K][i-1])/h**2+b*(u[K][i+1]-u[K][i-1])/(2*h))+u[K][i]
    t=tau*(K+1)

    u[K+1][0]=(-2*alpha*a*tau*u[K+1][1])/D-(alpha*u[K][0]*h**2)/D+(tau*h*(2*a-b*h)*(-math.exp(-a*t)*(math.cos(b*t)+math.sin(b*t))))/D
    u[K+1][n]=(2*gamma*a*tau*u[K+1][n-1])/C+(gamma*u[K][n]*h**2)/C+(tau*h*(2*a+b*h)*(math.exp(-a*t)*(math.cos(b*t)+math.sin(b*t))))/C

print('eroor explicit method 3')
xx=0
time=0
error1=math.exp(-a*time*tau)*math.cos(xx*h+b*time*tau)-u[time][xx]
xx=round(n/2)
time=round(k/2)
error2=math.exp(-a*time*tau)*math.cos(xx*h+b*time*tau)-u[time][xx]
xx=n
time=k
error3=math.exp(-a*time*tau)*math.cos(xx*h+b*time*tau)-u[time][xx]

print(error1,error2,error3)

draw(u,'explicit method 3')

# метод прогонки из прошлого семестра

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

#гибридная схема
def hybrid_method(tetta, name1,name2,name3):
    # 1 способ апроксимации
    
    K=0
    u=np.zeros((k+1,n+1))

    for i in range(n+1):
        x=l+i*h
        u[0][i]=math.cos(x)

    A=np.zeros((n+1,n+1))

    for K in range(0,k):
        A=np.zeros((n+1,n+1))
        B=np.zeros((n+1,1))
        t=tau*(K+1)
        A[0][0]=betta-alpha/h 
        A[0][1]=alpha/h 
        B[0][0]=-math.exp(-a*t)*(math.cos(b*t)+math.sin(b*t))

        A[n][n-1]=-gamma/h 
        A[n][n]=delta+gamma/h 
        B[n][0]=math.exp(-a*t)*(math.cos(b*t)+math.sin(b*t))

        for i in range(1,n):
            A[i][i-1]=-a*tetta*tau/h**2-b*tetta*tau/(2*h)
            A[i][i]=1+2*a*tetta*tau/h**2
            A[i][i+1]=-a*tetta*tau/h**2+b*tetta*tau/(2*h)
            B[i][0]=tau*(1-tetta)*(a*(u[K][i+1]-2*u[K][i]+u[K][i-1])/h**2+b*(u[K][i+1]-u[K][i-1])/(2*h))+u[K][i]
        x=running_method(A,B,n+1)
    
        for i in range(0,n+1):
            u[K+1][i]=x[i]
    
    print('eroor',name1)
    xx=0
    time=0
    error1=math.exp(-a*time*tau)*math.cos(xx*h+b*time*tau)-u[time][xx]
    xx=round(n/2)
    time=round(k/2)
    error2=math.exp(-a*time*tau)*math.cos(xx*h+b*time*tau)-u[time][xx]
    xx=n
    time=k
    error3=math.exp(-a*time*tau)*math.cos(xx*h+b*time*tau)-u[time][xx]

    print(error1,error2,error3)

    draw(u,name1)



    # 2 способ апроксимации

    K=0
    u=np.zeros((k+1,n+1))

    for i in range(n+1):
        x=l+i*h
        u[0][i]=math.cos(x)

    A=np.zeros((n+1,n+1))

    for K in range(0,k):
        A=np.zeros((n+1,n+1))
        B=np.zeros((n+1,1))

        for i in range(1,n):
            A[i][i-1]=-a*tetta*tau/h**2-b*tetta*tau/(2*h)
            A[i][i]=1+2*a*tetta*tau/h**2
            A[i][i+1]=-a*tetta*tau/h**2+b*tetta*tau/(2*h)
            B[i][0]=tau*(1-tetta)*(a*(u[K][i+1]-2*u[K][i]+u[K][i-1])/h**2+b*(u[K][i+1]-u[K][i-1])/(2*h))+u[K][i]
    
        t=tau*(K+1)
    
        A[0][0]= (-3*alpha+2*h*betta)-A[1][0]*(-alpha)/A[1][2]
        A[0][1]=4*alpha-A[1][1]*(-alpha)/A[1][2] 
        B[0][0]=-2*h*math.exp(-a*t)*(math.cos(b*t)+math.sin(b*t))-B[1][0]*(-alpha)/A[1][2]

        A[n][n-1]=-4*gamma-A[n-1][n-1]*gamma/A[n-1][n-2]
        A[n][n]=3*gamma+2*h*delta-A[n-1][n]*gamma/A[n-1][n-2]
        B[n][0]=2*h*math.exp(-a*t)*(math.cos(b*t)+math.sin(b*t))-B[n-1][0]*gamma/A[n-1][n-2]
        x=running_method(A,B,n+1)
    
        for i in range(0,n+1):
            u[K+1][i]=x[i]
    
    print('eroor',name2)
    xx=0
    time=0
    error1=math.exp(-a*time*tau)*math.cos(xx*h+b*time*tau)-u[time][xx]
    xx=round(n/2)
    time=round(k/2)
    error2=math.exp(-a*time*tau)*math.cos(xx*h+b*time*tau)-u[time][xx]
    xx=n
    time=k
    error3=math.exp(-a*time*tau)*math.cos(xx*h+b*time*tau)-u[time][xx]

    print(error1,error2,error3)
    draw(u,name2)

    # 3 способ апроксимации

    K=0
    u=np.zeros((k+1,n+1))

    for i in range(n+1):
        x=l+i*h
        u[0][i]=math.cos(x)

    A=np.zeros((n+1,n+1))

    for K in range(0,k):
        A=np.zeros((n+1,n+1))
        B=np.zeros((n+1,1))
        t=tau*(K+1)
        A[0][0]=-2*a*alpha/(h*(2*a-b*h))-alpha*h/(tau*(2*a-b*h))+betta
        A[0][1]=2*a*alpha/(h*(2*a-b*h)) 
        B[0][0]=-alpha*h/(tau*(2*a-b*h))*u[K][0]-math.exp(-a*t)*(math.cos(b*t)+math.sin(b*t))

        A[n][n-1]=-2*a*gamma/(h*(2*a+b*h)) 
        A[n][n]=2*a*gamma/(h*(2*a+b*h))+gamma*h/(tau*(2*a+b*h))+delta
        B[n][0]=gamma*h/(tau*(2*a+b*h))*u[K][n]+math.exp(-a*t)*(math.cos(b*t)+math.sin(b*t))

        for i in range(1,n):
            A[i][i-1]=-a*tetta*tau/h**2-b*tetta*tau/(2*h)
            A[i][i]=1+2*a*tetta*tau/h**2
            A[i][i+1]=-a*tetta*tau/h**2+b*tetta*tau/(2*h)
            B[i][0]=tau*(1-tetta)*(a*(u[K][i+1]-2*u[K][i]+u[K][i-1])/h**2+b*(u[K][i+1]-u[K][i-1])/(2*h))+u[K][i]
        x=running_method(A,B,n+1)
    
        for i in range(0,n+1):
            u[K+1][i]=x[i]
    
    print('eroor',name3)
    xx=0
    time=0
    error1=math.exp(-a*time*tau)*math.cos(xx*h+b*time*tau)-u[time][xx]
    xx=round(n/2)
    time=round(k/2)
    error2=math.exp(-a*time*tau)*math.cos(xx*h+b*time*tau)-u[time][xx]
    xx=n
    time=k
    error3=math.exp(-a*time*tau)*math.cos(xx*h+b*time*tau)-u[time][xx]

    print(error1,error2,error3)
    draw(u,name3)

    

hybrid_method(1,'implicit_method 1', 'implicit_method 2', 'implicit_method 3' )

hybrid_method(0.5,'Krank Nikolson  1', 'Krank Nikolson 2', 'Krank Nikolson 3' )




plt.show() 


