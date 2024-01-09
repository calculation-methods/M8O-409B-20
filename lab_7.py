import sympy as sp
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline 

N=10 # разбиения x 
M=10 # разбиения y 

L_x=0
R_x=math.pi/2

L_y=0
R_y=math.pi/2

# сетка
h_x=(R_x-L_x)/N 
h_y=(R_y-L_y)/M 

epsilon=0.00001

f=4
c=1.1
# f=1 - метод Либмана ; f=2 - метод Зейделя; f=3 - метод верхних релаксаций


def phi1(y):
    return math.exp(-y)*math.cos(y)

def phi2(y):
    return 0

def psi1(x):
    return math.cos(x)

def psi2(x):
    return 0

U=np.zeros((N+1,M+1)) 
for i in range(0,N+1):
    for j in range(0, M+1):
        x=L_x+i*h_x
        y=L_y+j*h_y
        U[i][j]=math.exp(-y)*math.cos(x)*math.cos(y)


def error(u_new,u_old):
    error=0 
    max1=0
    for k in range(0,N+1):
        for i in range(0, M+1):
            if abs(u_new[k][i]-u_old[k][i])>max1:
                max1=abs(u_new[k][i]-u_old[k][i])
            #error+=(u_new[k][i]-u_old[k][i])**2
    error=max1
    return(error)

def error3(u_new,u_old):
    error=0 
    max1=0
    for k in range(0,N+1):
        for i in range(0, M+1):
            
                error+=(u_new[k][i]-u_old[k][i])**2
            #error+=(u_new[k][i]-u_old[k][i])**2
    
    return (error/(M+1)/(N+1))

u=np.zeros((N+1,M+1))

def draw(u,name):

    TT=[j*h_y for j in range(M+1)]
    H=[L_x+j*h_x for j in range(N+1)]

    plt.figure(name)
    plt.subplot(1,2,1)
    x = np.linspace(0,R_x, 1000)
    plt.plot(x, np.exp(0)*np.cos(x)*np.cos(0),label='U(x,y)  y=0')
    plt.plot(x, np.exp(-round(M/2)*h_y)*np.cos(x)*np.cos(round(M/2)*h_y),label='U(x,y)  y=R_y/2')
    plt.plot(x, np.exp(-R_y)*(np.cos(x))*np.cos(R_y),label='U(x,y)  y=R_y')
    plt.title('u(x), y=const')
    #plt.plot(H,u[0],'ro')
    xnew = np.linspace(0, R_y, 1000) 

    ut0=[u[j][0] for j in range(N+1)]
    ut1=[u[j][round(M/2)] for j in range(N+1)]
    ut2=[u[j][M] for j in range(N+1)]

    gfg = make_interp_spline(H, ut0, k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='y=0')   
    gfg = make_interp_spline(H, ut1, k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='y=R_y/2')  
    gfg = make_interp_spline(H, ut2, k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='y=R_y')  
    plt.legend(fontsize=7, loc=0)

    plt.subplot(1,2,2)
    t = np.linspace(0,R_y, 300)
    plt.plot(t, np.exp(-t)*np.cos(0)*np.cos(t),label='U(x,y)  x=0')
    plt.plot(t, np.exp(-t)*np.cos(round(N/2)*h_x)*np.cos(t),label='U(x,y)  x=R_x/2')
    plt.plot(t, np.exp(-t)*np.cos(R_x)*np.cos(t),label='U(x,y)  x=R_x')
    plt.title('u(y), x=const')
    #plt.plot(H,u[0],'ro')
    xnew = np.linspace(0, R_y, 300) 

   
    gfg = make_interp_spline(TT, u[0], k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='x=0')   
    gfg = make_interp_spline(TT, u[round(N/2)], k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='x=R_x/2')  
    gfg = make_interp_spline(TT, u[N], k=3)   
    y_new = gfg(xnew)   
    plt.plot(xnew, y_new,label='x=R_x')  
    plt.legend(fontsize=7, loc=0)


k=0

for i in range(0,N+1):
    x=L_x+i*h_x
    u[i][0]=psi1(x) 
    u[i][M]=psi2(x)


for j in range(0,M+1):
    y=L_y+j*h_y
    u[0][j]=phi1(y)
    u[N][j]=phi2(y)



for i in range(1,N):
    for j in range(1,M):
        x=L_x+i*h_x
        y=L_y+j*h_y
        Cx=phi1(x)*(R_x-x)/(R_x-L_x)+phi2(x)*(x-L_x)/(R_x-L_x) 
        Cy=psi1(x)*(R_y-y)/(R_y-L_y)+psi2(x)*(y-L_y)/(R_y-L_y) 
        u[i][j]=(Cy+Cx)/2

u_new=np.array(u) 

error1=1






def solve(u,f):
    u_new=np.array(u) 
    k=0
    error1=1
    while epsilon<error1:

        k+=1
        u_old=np.array(u_new)
    
    
        for i in range(1,N):
            if f==1:
                for j in range(1,M):
                    u_new[i][j]=(h_y**2*(u_old[i+1][j]+u_old[i-1][j])+h_x**2*(u_old[i][j-1]+u_old[i][j+1])+h_x**2*h_y*(u_old[i][j+1]-u_old[i][j-1]))/(2*h_y**2+2*h_y**2-3*h_x**2*h_y**2)
            elif f==2:
                for j in range(1,M):
                    u_new[i][j]=(h_y**2*(u_old[i+1][j]+u_new[i-1][j])+h_x**2*(u_new[i][j-1]+u_old[i][j+1])+h_x**2*h_y*(u_old[i][j+1]-u_new[i][j-1]))/(2*h_y**2+2*h_y**2-3*h_x**2*h_y**2)
            elif f==3:
                for j in range(1,M):
                    uk=(h_y**2*(u_old[i+1][j]+u_new[i-1][j])+h_x**2*(u_new[i][j-1]+u_old[i][j+1])+h_x**2*h_y*(u_old[i][j+1]-u_new[i][j-1]))/(2*h_y**2+2*h_y**2-3*h_x**2*h_y**2)
                    u_new[i][j]=c*uk+(1-c)*u_old[i][j]
            

        error1=error(u_new,u_old)
   
    error2=error(u_new,U)
    print("k=",k)
    print("max error",error2)
    error4=error3(u_new,U)
    print("error",error4)
    return u_new

print("Simple iteration")
draw(solve(u,1),'Simple iteration')

print("\nZeydel ")
draw(solve(u,2),'Zeydel')

print("\nRelaxation")
draw(solve(u,3),'Relaxation')

plt.show()