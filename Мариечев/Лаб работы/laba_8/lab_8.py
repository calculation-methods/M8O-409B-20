# noinspection LossyEncoding
import sympy as sp
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline 

a=2
b=1
nu=1


N=30 # ��������� x 
M=30 # ��������� y 
K=30 # ��������� t 
T=0.01

L1=0
R1=math.pi/2

L2=0
R2=math.pi

# �����
hx=(R1-L1)/N 
hy=(R2-L2)/M 
tau=(T-0)/K

# f=1 - ������������� ������� ���������� ������� � ������ �������� ; f=2 - ������������� ������� ���������� ������� �� ������ �������� 


def phi1(y,t):
    return 0

def phi2(y,t):
    return math.sin(y)*math.sin(nu*t)

def phi3(x,t):
    return 0

def phi4(x,t):
    return -math.sin(x)*math.sin(nu*t)


def psi1(x,y):
    return 0


def f(x,y,t):
    return math.sin(x)*math.sin(y)*(nu*math.cos(nu*t)+(a+b)*math.sin(nu*t))

U=np.zeros((N+1,M+1,K+1)) 
for i in range(0,N+1):
    for j in range(0, M+1):
        for k in range(0, K+1):
            x=i*hx
            y=j*hy
            t=k*tau
            U[i][j][k]=math.sin(x)*math.sin(y)*math.sin(nu*t)

def error1(u,U):
    error=0 
    max1=0
    for i in range(0,N+1):
        for j in range(0, M+1):
            for k in range(0,K+1):
                if abs(u[i][j][k]-U[i][j][k])>max1:
                    max1=abs(u[i][j][k]-U[i][j][k])
            #error+=(u_new[k][i]-u_old[k][i])**2
    error=max1
    return(error)

def error2(u,U):
    error=0 
    for i in range(0,N+1):
        for j in range(0, M+1):
            for k in range(0,K+1):
                error+=(u[i][j][k]-U[i][j][k])**2
            #error+=(u_new[k][i]-u_old[k][i])**2
    
    return (error/(M+1)/(N+1)/(K+1))


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

def draw(u,name):

    TT=[j*tau for j in range(K+1)]
    Hx=[j*hx for j in range(N+1)]
    Hy=[j*hy for j in range(M+1)]

    plt.figure(name)
    plt.subplot(1,3,1)
    x = np.linspace(0,R1, 1000)
    plt.plot(x, np.sin(x)*np.sin(round(M/2)*hy)*np.sin(nu*0),label='U(x,y,t)  t=0, y=R2/2')
    plt.plot(x, np.sin(x)*np.sin(round(M/2)*hy)*np.sin(nu*round(K/2)*tau),label='U(x,y,y)  t=T/2, y=R2/2')
    plt.plot(x, np.sin(x)*np.sin(round(M/2)*hy)*np.sin(nu*T),label='U(x,y,t)  t=T, y=R2/2')
    '''
    plt.plot(x, np.sin(x)*np.sin(R2)*np.sin(nu*0),label='U(x,y,t)  t=0, y=R2')
    plt.plot(x, np.sin(x)*np.sin(R2)*np.sin(nu*round(K/2)*tau),label='U(x,y)  t=T/2, y=R2')
    plt.plot(x, np.sin(x)*np.sin(R2)*np.sin(nu*T),label='U(x,y)  t=T, y=R2')
    '''

    plt.title('u(x), y=const, t=const')
    #plt.plot(H,u[0],'ro')
    xnew = np.linspace(0, R1, 1000) 

    ut0=[u[j][round(M/2)][0] for j in range(N+1)]
    ut1=[u[j][round(M/2)][round(K/2)] for j in range(N+1)]
    ut2=[u[j][round(M/2)][K] for j in range(N+1)]

    
    
    gfg = make_interp_spline(Hx, ut0, k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='t=0, y=R2/2')   
    gfg = make_interp_spline(Hx, ut1, k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='t=T/2, y=R2/2')  
    gfg = make_interp_spline(Hx, ut2, k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='t=T, y=R2/2')  
    '''
    gfg = make_interp_spline(Hx, ut3, k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='t=0, y=R2')   
    gfg = make_interp_spline(Hx, ut4, k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='t=T/2, y=R2')  
    gfg = make_interp_spline(Hx, ut5, k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='t=T, y=R2')  
    '''
    plt.legend(fontsize=7, loc=0)
    
    plt.subplot(1,3,2)
    y = np.linspace(0,R2, 300)
    plt.plot(y, np.sin(y)*np.sin(round(N/2)*hx)*np.sin(nu*0),label='U(x,y,t)  t=0, x=R1/2')
    plt.plot(y, np.sin(y)*np.sin(round(N/2)*hx)*np.sin(nu*round(K/2)*tau),label='U(x,y,t)  t=T/2, x=R1/2')
    plt.plot(y, np.sin(y)*np.sin(round(N/2)*hx)*np.sin(nu*T),label='U(x,y,t)  t=T, x=R1/2')
    plt.title('u(y), x=const, t=const')
    #plt.plot(H,u[0],'ro')
    xnew = np.linspace(0, R2, 300) 

    ut3=[u[round(N/2)][j][0] for j in range(M+1)]
    ut4=[u[round(N/2)][j][round(K/2)] for j in range(M+1)]
    ut5=[u[round(N/2)][j][K] for j in range(M+1)]
   
    gfg = make_interp_spline(Hy, ut3, k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='t=0, x=R1/2')   
    gfg = make_interp_spline(Hy, ut4, k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='t=T/2, x=R1/2')  
    gfg = make_interp_spline(Hy, ut5, k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='t=T, x=R1/2')  
    plt.legend(fontsize=7, loc=0)

    plt.subplot(1,3,3)
    t = np.linspace(0,T, 300)
    plt.plot(t, np.sin(round(M/2)*hy)*np.sin(round(N/2)*hx)*np.sin(nu*t),label='U(x,y,t)  y=R2/2, x=R1/2')
    plt.plot(t, np.sin(R1)*np.sin(round(M/2)*hy)*np.sin(nu*t),label='U(x,y,t)  y=R2/2, x=R1')
    plt.plot(t, np.sin(R1)*np.sin(0)*np.sin(nu*t),label='U(x,y,t)  y=0, x=R1')
    plt.title('u(y), x=const, t=const')
    #plt.plot(H,u[0],'ro')
    xnew = np.linspace(0, T, 300) 

    ut6=[u[round(N/2)][round(M/2)][j] for j in range(K+1)]
    ut7=[u[N][round(M/2)][j] for j in range(K+1)]
    ut8=[u[N][0][j] for j in range(K+1)]
   
    gfg = make_interp_spline(TT, ut6, k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='y=R2/2, x=R1/2')   
    gfg = make_interp_spline(TT, ut7, k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='y=R2/2, x=R1')  
    gfg = make_interp_spline(TT, ut8, k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='y=0, x=R1')  
    plt.legend(fontsize=7, loc=0)





   

def solve(f1,name):
    u=np.zeros((N+1,M+1,K+1))

    for i in range(0,N+1):
        for j in range(0, M+1):
            x=i*hx
            y=j*hy
            u[i][j][0]=psi1(x,y)


    k=0
    for k in range(0,K):
        u2=np.zeros((N+1,M+1))
        for j in range(1, M):
    
            y=j*hy
            t=(k+1/2)*tau
            A=np.zeros((N+1,N+1))
            B=np.zeros((N+1,1))
            A[0][0]=1
            A[0][1]=0 
            B[0][0]=phi1(y,(k+1/2)*tau)

            A[N][N-1]=0
            A[N][N]=1
            B[N][0]=phi2(y,(k+1/2)*tau)
            if f1==1:
                for i in range(1,N):
                    x=i*hx
                    A[i][i-1]=-a/hx**2
                    A[i][i]=2/tau+2*a/hx**2
                    A[i][i+1]=-a/hx**2
                    B[i][0]=2*u[i][j][k]/tau+b*(u[i][j+1][k]-2*u[i][j][k]+u[i][j-1][k])/hy**2+f(x,y,(k+1/2)*tau)
            elif f1==2:
                for i in range(1,N):
                    x=i*hx
                    A[i][i-1]=-a/hx**2
                    A[i][i]=1/tau+2*a/hx**2
                    A[i][i+1]=-a/hx**2
                    B[i][0]=u[i][j][k]/tau+f(x,y,k*tau)/2
            x=running_method(A,B,N+1)
        
            for i in range(0,N+1):
                u2[i][j]=x[i]

        for i in range(0,N+1):
            u2[i][0]=phi3(i*hx,(k+1/2)*tau)
            u2[i][M]=phi4(i*hx,(k+1/2)*tau)*hy+u2[i][N-1]
        

        for i in range(1, N):
    
            x=i*hx
            t=(k+1)*tau
            A=np.zeros((M+1,M+1))
            B=np.zeros((M+1,1))
            A[0][0]=1
            A[0][1]=0 
            B[0][0]=phi3(x,t)

            A[M][M-1]=-1/hy
            A[M][M]=1/hy
            B[M][0]=phi4(x,t)
            if f1==1:
                for j in range(1,M):
                    y=j*hy
                    A[j][j-1]=-b/hy**2
                    A[j][j]=2/tau+2*b/hy**2
                    A[j][j+1]=-b/hy**2
                    B[j][0]=2*u2[i][j]/tau+a*(u2[i+1][j]-2*u2[i][j]+u2[i-1][j])/hx**2+f(x,y,(k+1/2)*tau)
            elif f1==2:
                for j in range(1,M):
                    y=j*hy
                    A[j][j-1]=-b/hy**2
                    A[j][j]=1/tau+2*b/hy**2
                    A[j][j+1]=-b/hy**2
                    B[j][0]=u2[i][j]/tau++f(x,y,(k+1)*tau)/2
        
            y=running_method(A,B,M+1)
        
            for j in range(0,M+1):
                u[i][j][k+1]=y[j]

        for j in range(0,M+1):
            u[0][j][k+1]=phi1(j*hy,t)
            u[N][j][k+1]=phi2(j*hy,t)



    print('max error  ',error1(u,U))
    print('avarege error  ', error2(u,U))
    draw(u,name)

print('Method 1')
solve(1,'Method 1')

print('Method 2')
solve(2,'Method 2')

plt.show()