import tkinter as tk
import matplotlib.pyplot as plt
import numpy as np
import KP_function1 as kp1
import KP_function2 as kp2
import KP_function3 as kp3
import KP as kp
from tkinter import ttk, Tk
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
NavigationToolbar2Tk)


razb=0
N=1000
flag=0
legend=0
flagP=0
choice_int=""
choice_meth=''
sc=0
Pnts=[]
Methods = []


def Scale(ax, maxx, minx, maxy, miny, sc):
    if sc==0:
        dx=(maxx-minx)*0.05
        dy=(maxy-miny)*0.05
        ax.set_xlim([minx-dx, maxx+dx])
        ax.set_ylim([miny-dy, maxy+dy])
    else:
        if maxx-minx>maxy-miny:
            d=(maxx-minx)*0.05
            half=miny+(maxy-miny)/2
            ax.set_xlim([minx-d, maxx+d])
            ax.set_ylim([half-(maxx-minx)/2-d, half+(maxx-minx)/2+d])
        else:
            d = (maxy - miny) * 0.05
            half=minx+(maxx-minx)/2
            print(half)
            ax.set_ylim([miny - d, maxy + d])
            ax.set_xlim([half - (maxy - miny) / 2 - d, half + (maxy - miny) / 2 + d])

def Scale3d(ax, maxx, minx, maxy, miny, maxz, minz, sc):
    if sc==0:
        dx=(maxx-minx)*0.05
        dy=(maxy-miny)*0.05
        dz=(maxz-minz)*0.05
        ax.set_xlim([minx-dx, maxx+dx])
        ax.set_ylim([miny-dy, maxy+dy])
        ax.set_zlim([minz-dz, maxz+dz])
    else:
        if maxx-minx>maxy-miny and maxx-minx>maxz-minz:
            d=(maxx-minx)*0.05
            half=miny+(maxy-miny)/2
            half2=minz+(maxz-minz)/2
            ax.set_xlim([minx-d, maxx+d])
            ax.set_ylim([half-(maxx-minx)/2-d, half+(maxx-minx)/2+d])
            ax.set_zlim([half2-(maxx-minx)/2-d, half2+(maxx-minx)/2+d])
        elif maxy-miny>maxx-minx and maxy-miny>maxz-minz:
            d = (maxy - miny) * 0.05
            half=minx+(maxx-minx)/2
            half2 = minz + (maxz - minz) / 2
            ax.set_ylim([miny - d, maxy + d])
            ax.set_xlim([half - (maxy - miny) / 2 - d, half + (maxy - miny) / 2 + d])
            ax.set_zlim([half2-(maxy-miny)/2-d, half2+(maxy-miny)/2+d])
        else:
            d = (maxz - minz) * 0.05
            half = minx + (maxx - minx) / 2
            half2 = miny + (maxy - miny) / 2
            ax.set_zlim([minz - d, maxz + d])
            ax.set_xlim([half - (maxz - minz) / 2 - d, half + (maxz - minz) / 2 + d])
            ax.set_ylim([half2 - (maxz - minz) / 2 - d, half2 + (maxz - minz) / 2 + d])

def plot(razb, ax, Pnts, flag, flagP, choice_int):
    ax.clear()
    fig.clear()
    print(choice_int)
    if choice_int=="Одномерный интеграл":
        ax = fig.add_subplot(111)
        x = np.linspace(kp1.bx(), kp1.tx(), razb)
        ax.axhline(y=0, color='k', linewidth=0.5)
        ax.axvline(x=0, color='k', linewidth=0.5)
        y = kp1.function(x)
        ax.plot(x, y, label="Подинтегральная функция")
        ax.plot([kp1.bx(), kp1.bx()], [np.min(y), np.max(y)], color="olive", label="Ограничения по x")
        ax.plot([kp1.tx(), kp1.tx()], [np.min(y), np.max(y)], color="olive")
        Scale(ax, np.max(x), np.min(x), np.max(y), np.min(y), sc)
    elif choice_int=="Двумерный интеграл":
        print("yes")
        if kp2.function(2.4, 2.4)==1:
            ax = fig.add_subplot(111)
            if kp2.ty(1)==kp2.ty(2):
                ax.axhline(y=0, color='k', linewidth=0.5)
                ax.axvline(x=0, color='k', linewidth=0.5)
                ax.axhline(y=kp2.ty(0), color="orange", label="Ограничения по y")
                ax.axhline(y=kp2.by(0), color="orange")
                y = np.linspace(kp2.by(0), kp2.ty(0), razb)
                if kp2.tx(2)==kp2.tx(3):
                    x = np.full(shape=razb, fill_value=kp2.tx(2))
                else:
                    x = kp2.tx(y)
                maxx1=np.max(x); minx1=np.min(x); maxy1=np.max(y); miny1=np.min(y)
                ax.plot(x, y, color="olive", label="Ограничения по x")
                if kp2.bx(2)==kp2.bx(3):
                    x = np.full(shape=razb, fill_value=kp2.bx(2))
                else:
                    x = kp2.bx(y)
                ax.plot(x, y, color="olive")
                maxx2 = np.max(x); minx2 = np.min(x); maxy2 = np.max(y); miny2 = np.min(y)
                Scale(ax, max(maxx1, maxx2), min(minx1, minx2), max(maxy1, maxy2), min(miny1, miny2), sc)
            elif kp2.tx(1)==kp2.tx(2):
                ax.axhline(y=0, color='k', linewidth=0.5)
                ax.axvline(x=0, color='k', linewidth=0.5)
                ax.axvline(x=kp2.tx(0), color="olive", label="Ограничения по x")
                ax.axvline(x=kp2.bx(0), color="olive")
                x = np.linspace(kp2.bx(0), kp2.tx(0), razb)
                if kp2.ty(2)==kp2.ty(3):
                    y = np.full(shape=razb, fill_value=kp2.ty(2))
                else:
                    y = kp2.ty(x)
                maxx1 = np.max(x); minx1 = np.min(x); maxy1 = np.max(y); miny1 = np.min(y)
                ax.plot(x, y, color="orange", label="Ограничения по y")
                if kp2.by(2)==kp2.by(3):
                    y=np.full(shape=razb, fill_value=kp2.by(2))
                else:
                    y = kp2.by(x)
                ax.plot(x, y, color="orange")
                maxx2 = np.max(x); minx2 = np.min(x); maxy2 = np.max(y); miny2 = np.min(y)
                Scale(ax, max(maxx1, maxx2), min(minx1, minx2), max(maxy1, maxy2), min(miny1, miny2), sc)
        else:
            if kp2.ty(0)==kp2.ty(2):
                y = np.linspace(kp2.by(0), kp2.ty(0), int(razb))
                x = np.linspace(np.min(kp2.bx(y)), np.max(kp2.tx(y)), int(razb))
                X, Y = np.meshgrid(x, y)
                ax = fig.add_subplot(111, projection='3d')
                ax.plot_wireframe(X, Y, np.zeros((len(x), len(y))), color="black")
                ax.plot_wireframe(X, Y, kp2.function(X, Y), color="blue", label="Подинтегральная функция")
                z = np.linspace(np.min(kp2.function(X,Y)), np.max(kp2.function(X,Y)), int(razb))
                x1 = kp2.tx(y)
                Y, Z = np.meshgrid(y, z)
                maxx1 = np.max(x1); minx1 = np.min(x1); maxy1 = np.max(y); miny1 = np.min(y)
                ax.plot_wireframe(kp2.tx(Y), Y, Z, color="olive", label="Ограничение по x")
                ax.plot_wireframe(kp2.bx(Y), Y, Z, color="olive")
                x2 = kp2.bx(y)
                X, Z = np.meshgrid(x, z)
                ax.plot_wireframe(X, kp2.ty(X), Z, color="orange", label="Ограничение по y")
                ax.plot_wireframe(X, kp2.by(X), Z, color="orange")
                maxx2 = np.max(x2); minx2 = np.min(x2); maxy2 = np.max(y); miny2 = np.min(y)
                maxz=kp.Max_Min2(x, y, kp2.function)[0]
                minz=kp.Max_Min2(x, y, kp2.function)[1]
                Scale3d(ax, max(maxx1, maxx2), min(minx1, minx2), max(maxy1, maxy2), min(miny1, miny2), maxz, minz, sc)
            elif kp2.tx(2)==kp2.tx(3):
                x = np.linspace(kp2.bx(0), kp2.tx(0), int(razb))
                y = np.linspace(np.min(kp2.by(x)), np.max(kp2.ty(x)), int(razb))
                X, Y = np.meshgrid(x, y)
                ax = fig.add_subplot(111, projection='3d')
                ax.plot_wireframe(X, Y, np.zeros((len(x), len(y))), color="black")
                ax.plot_wireframe(X, Y, kp2.function(X, Y), color="blue", label="Подинтегральная функция")
                z = np.linspace(np.min(kp2.function(X, Y)), np.max(kp2.function(X, Y)), int(razb))
                y1 = kp2.ty(x)
                X, Z = np.meshgrid(y, z)
                maxx1 = np.max(x); minx1 = np.min(x); maxy1 = np.max(y1); miny1 = np.min(y1)
                ax.plot_wireframe(X, kp2.ty(X), Z, color="orange", label="Ограничение по y")
                y2 = kp2.by(x)
                ax.plot_wireframe(X, kp2.by(X), Z, color="orange")
                maxx2 = np.max(x); minx2 = np.min(x); maxy2 = np.max(y2); miny2 = np.min(y2)
                maxz = kp.Max_Min2(x, y, kp2.function)[0]
                minz = kp.Max_Min2(x, y, kp2.function)[1]
                Scale3d(ax, max(maxx1, maxx2), min(minx1, minx2), max(maxy1, maxy2), min(miny1, miny2), maxz, minz, sc)
    elif choice_int=="Трехмерный интеграл":
        #ax = fig.add_subplot(111)
        ax = fig.add_subplot(111, projection='3d')
        if kp3.order()[0]=='x':
            if kp3.order()[1]=='y':
                x = np.linspace(kp3.bx(0, 0), kp3.tx(0, 0), int(razb))
                y = np.linspace(np.min(kp3.by(x, 0)), np.max(kp3.ty(x, 0)), int(razb))
                X, Y = np.meshgrid(x, y)
                ax.plot_wireframe(X, Y, np.zeros((len(x), len(y))), color="black")
                z = np.linspace(np.min(kp3.bz(X, Y)), np.max(kp3.tz(X, Y)), int(razb))
                y1 = kp3.ty(x, 0)
                Y, Z = np.meshgrid(y, z)
                ax.plot_wireframe(kp3.tx(Y, 0), Y, Z, color="olive", label="Ограничение по x")
                ax.plot_wireframe(kp3.bx(Y, 0), Y , Z, color="olive")
                X, Z = np.meshgrid(x, z)
                maxx1 = np.max(x)
                minx1 = np.min(x)
                maxy1 = np.max(y1)
                miny1 = np.min(y1)
                ax.plot_wireframe(X, kp3.ty(X, 0), Z, color="orange", label="Ограничение по y")
                y2 = kp3.by(x, 0)
                ax.plot_wireframe(X, kp3.by(X, 0), Z, color="orange")
                X, Y = np.meshgrid(x, y)
                ax.plot_wireframe(X, Y, kp.Border(X, Y, kp3.tz), color="green", label="Ограничение по z")
                ax.plot_wireframe(X, Y, kp.Border(X, Y, kp3.bz), color="green")
                maxx2 = np.max(x)
                minx2 = np.min(x)
                maxy2 = np.max(y2)
                miny2 = np.min(y2)
                maxz = np.max(kp3.tz(X, Y))
                minz = np.min(kp3.bz(X, Y))
                Scale3d(ax, max(maxx1, maxx2), min(minx1, minx2), max(maxy1, maxy2), min(miny1, miny2), maxz, minz, sc)
            else:
                x = np.linspace(kp3.bx(0, 0), kp3.tx(0, 0), int(razb))
                z = np.linspace(np.min(kp3.bz(x, 0)), np.max(kp3.tz(x, 0)), int(razb))
                X, Z = np.meshgrid(x, z)
                #ax.plot_wireframe(X, Y, np.zeros((len(x), len(y))), color="black")
                y = np.linspace(np.min(kp3.by(X, Z)), np.max(kp3.ty(X, Z)), int(razb))
                z1 = kp3.tz(x, 0)
                Y, Z = np.meshgrid(y, z)
                ax.plot_wireframe(kp3.tx(Y, 0), Y, Z, color="olive", label="Ограничение по x")
                ax.plot_wireframe(kp3.bx(Y, 0), Y, Z, color="olive")
                X, Z = np.meshgrid(x, z)
                maxx1 = np.max(x)
                minx1 = np.min(x)
                maxz1 = np.max(z1)
                minz1 = np.min(z1)
                ax.plot_wireframe(X, kp3.ty(X, Z), Z, color="orange", label="Ограничение по y")
                ax.plot_wireframe(X, kp3.by(X, Z), Z, color="orange")
                z2 = kp3.bz(x, 0)
                X, Y = np.meshgrid(x, y)
                ax.plot_wireframe(X, Y, kp.Border(X, Y, kp3.tz), color="green", label="Ограничение по z")
                ax.plot_wireframe(X, Y, kp.Border(X, Y, kp3.bz), color="green")
                maxx2 = np.max(x)
                minx2 = np.min(x)
                maxz2 = np.max(z2)
                minz2 = np.min(z2)
                maxy = np.max(kp3.ty(X, Z))
                miny = np.min(kp3.by(X, Z))
                Scale3d(ax, max(maxx1, maxx2), min(minx1, minx2), maxy, miny, max(maxz1, maxz2), min(minz1, minz2), sc)
        if kp3.order()[0]=='y':
            if kp3.order()[1]=='x':
                y = np.linspace(kp3.by(0, 0), kp3.ty(0, 0), int(razb))
                x = np.linspace(np.min(kp3.bx(y, 0)), np.max(kp3.tx(y, 0)), int(razb))
                X, Y = np.meshgrid(x, y)
                ax.plot_wireframe(X, Y, np.zeros((len(x), len(y))), color="black")
                z = np.linspace(np.min(kp3.bz(X, Y)), np.max(kp3.tz(X, Y)), int(razb))
                x1 = kp3.tx(y, 0)
                Y, Z = np.meshgrid(y, z)
                ax.plot_wireframe(X, kp3.ty(X, 0), Z, color="orange", label="Ограничение по y")
                ax.plot_wireframe(X, kp3.by(X, 0), Z, color="orange")
                X, Z = np.meshgrid(y, z)
                maxy1 = np.max(y)
                miny1 = np.min(y)
                maxx1 = np.max(x1)
                minx1 = np.min(x1)
                ax.plot_wireframe(kp3.tx(Y, 0), Y, Z, color="olive", label="Ограничение по x")
                ax.plot_wireframe(kp3.bx(Y, 0), Y, Z, color="olive")
                x2 = kp3.bx(y, 0)
                X, Y = np.meshgrid(x, y)
                ax.plot_wireframe(X, Y, kp.Border(X, Y, kp3.tz), color="green", label="Ограничение по z")
                ax.plot_wireframe(X, Y, kp.Border(X, Y, kp3.bz), color="green")
                maxy2 = np.max(y)
                miny2 = np.min(y)
                maxx2 = np.max(x2)
                minx2 = np.min(x2)
                maxz = np.max(kp3.tz(X, Y))
                minz = np.min(kp3.bz(X, Y))
                Scale3d(ax, max(maxx1, maxx2), min(minx1, minx2), max(maxy1, maxy2), min(miny1, miny2), maxz, minz, sc)
            else:
                y = np.linspace(kp3.by(0, 0), kp3.ty(0, 0), int(razb))
                z = np.linspace(np.min(kp3.bz(0, y)), np.max(kp3.tz(0, y)), int(razb))
                Y, Z = np.meshgrid(y, z)
                #ax.plot_wireframe(X, Y, np.zeros((len(x), len(y))), color="black")
                x = np.linspace(np.min(kp3.bx(Y, Z)), np.max(kp3.tx(Y, Z)), int(razb))
                z1 = kp3.tz(0, y)
                X, Z = np.meshgrid(x, z)
                ax.plot_wireframe(X, kp3.ty(X, 0), Z, color="orange", label="Ограничение по y")
                ax.plot_wireframe(X, kp3.by(X, 0), Z, color="orange")
                Y, Z = np.meshgrid(y, z)
                maxy1 = np.max(y)
                miny1 = np.min(y)
                maxz1 = np.max(z1)
                minz1 = np.min(z1)
                ax.plot_wireframe(kp3.tx(Y, Z), Y, Z, color="olive", label="Ограничение по x")
                ax.plot_wireframe(kp3.bx(Y, Z), Y, Z, color="olive")
                z2 = kp3.bx(0, y)
                X, Y = np.meshgrid(x, y)
                ax.plot_wireframe(X, Y, kp.Border(X, Y, kp3.tz), color="green", label="Ограничение по z")
                ax.plot_wireframe(X, Y, kp.Border(X, Y, kp3.bz), color="green")
                maxy2 = np.max(y)
                miny2 = np.min(y)
                maxz2 = np.max(z2)
                minz2 = np.min(z2)
                maxx = np.max(kp3.tx(Y, Z))
                minx = np.min(kp3.bx(Y, Z))
                Scale3d(ax, maxx, minx, max(maxy1, maxy2), min(miny1, miny2), max(maxz1, maxz2), min(minz1, minz2), sc)
        if kp3.order()[0]=='z':
            if kp3.order()[1]=='y':
                z = np.linspace(kp3.bz(0, 0), kp3.tz(0, 0), int(razb))
                y = np.linspace(np.min(kp3.by(0, z)), np.max(kp3.ty(0, z)), int(razb))
                Y, Z = np.meshgrid(y, z)
                x = np.linspace(np.min(kp3.bx(Y, Z)), np.max(kp3.tx(Y, Z)), int(razb))
                y1 = kp3.ty(0, z)
                Y, Z = np.meshgrid(y, z)
                ax.plot_wireframe(kp3.tx(Y, Z), Y, Z, color="olive", label="Ограничение по x")
                ax.plot_wireframe(kp3.bx(Y, Z), Y, Z, color="olive")
                X, Z = np.meshgrid(x, z)
                maxz1 = np.max(z)
                minz1 = np.min(z)
                maxy1 = np.max(y1)
                miny1 = np.min(y1)
                ax.plot_wireframe(X, kp3.ty(0, Z), Z, color="orange", label="Ограничение по y")
                ax.plot_wireframe(X, kp3.by(0, Z), Z, color="orange")
                y2 = kp3.by(x, 0)
                X, Y = np.meshgrid(x, y)
                ax.plot_wireframe(X, Y, np.zeros((len(x), len(y))), color="black")
                ax.plot_wireframe(X, Y, kp.Border(X, Y, kp3.tz), color="green", label="Ограничение по z")
                ax.plot_wireframe(X, Y, kp.Border(X, Y, kp3.bz), color="green")
                maxz2 = np.max(z)
                minz2 = np.min(z)
                maxy2 = np.max(y2)
                miny2 = np.min(y2)
                maxx = np.max(kp3.tx(X, Y))
                minx = np.min(kp3.bx(X, Y))
                Scale3d(ax, maxx, minx, max(maxy1, maxy2), min(miny1, miny2), max(maxz1, maxz2), min(minz1, minz2), sc)
            else:
                z = np.linspace(kp3.bz(0, 0), kp3.tz(0, 0), int(razb))
                x = np.linspace(np.min(kp3.bx(0, z)), np.max(kp3.tx(0, z)), int(razb))
                print(x[0], x[-1])
                X, Z = np.meshgrid(x, z)
                y = np.linspace(np.min(kp3.by(X, Z)), np.max(kp3.ty(X, Z)), int(razb))
                x1 = kp3.tx(0, z)
                ax.plot_wireframe(X, kp3.ty(X, Z), Z, color="orange", label="Ограничение по y")
                ax.plot_wireframe(X, kp3.by(X, Z), Z, color="orange")
                X, Y = np.meshgrid(x, y)
                ax.plot_wireframe(X, Y, np.zeros((len(x), len(y))), color="black")
                Y, Z = np.meshgrid(y, z)
                ax.plot_wireframe(kp3.tx(0, Z), Y, Z, color="olive", label="Ограничение по x")
                ax.plot_wireframe(kp3.bx(0, Z), Y, Z, color="olive")
                X, Z = np.meshgrid(y, z)
                maxz1 = np.max(z)
                minz1 = np.min(z)
                maxx1 = np.max(x1)
                minx1 = np.min(x1)
                x2 = kp3.bx(0, z)
                X, Y = np.meshgrid(x, y)
                ax.plot_wireframe(X, Y, kp.Border(X, Y, kp3.tz), color="green", label="Ограничение по z")
                ax.plot_wireframe(X, Y, kp.Border(X, Y, kp3.bz), color="green")
                maxz2 = np.max(z)
                minz2 = np.min(z)
                maxx2 = np.max(x2)
                minx2 = np.min(x2)
                maxy = np.max(kp3.ty(X, Y))
                miny = np.min(kp3.by(X, Y))
                Scale3d(ax, max(maxx1, maxx2), min(minx1, minx2), maxy, miny, max(maxz1, maxz2), min(minz1, minz2), sc)


    if flag==1:
        if len(Pnts)==3:
            plt.ion()
            X=[]; Y=[]
            for i in range(N):
                if Pnts[2][i]==1:
                    X.append(Pnts[0][i])
                    Y.append(Pnts[1][i])
            ax.scatter(X, Y, s=1, color="gray")
            plt.ioff()
        else:
            plt.ion()
            X = []; Y = []; Z=[]
            for i in range(N):
                if Pnts[3][i] == 1:
                    X.append(Pnts[0][i])
                    Y.append(Pnts[1][i])
                    Z.append(Pnts[2][i])
            ax.scatter(X, Y, Z, s=2, color="gray")
            plt.ioff()
    if flagP==1:
        if len(Pnts) == 3:
            plt.ion()
            X = []; Y = []
            for i in range(N):
                if Pnts[2][i] == 0:
                    X.append(Pnts[0][i])
                    Y.append(Pnts[1][i])
            ax.scatter(X, Y, s=1, color="red")
            plt.ioff()
        else:
            plt.ion()
            X = []; Y = []; Z = []
            for i in range(N):
                if Pnts[3][i] == 0:
                    X.append(Pnts[0][i])
                    Y.append(Pnts[1][i])
                    Z.append(Pnts[2][i])
            ax.scatter(X, Y, Z, s=2, color="red")
            plt.ioff()
    if legend==1:
        ax.legend()
    #ax.set_title("График подинтегрального выражения")
    ax.grid()
    canvas.draw()


def checkbutton_changed():
    global flag
    if enabled.get() == 1:
        flag=1
        plot(razb, ax, Pnts, flag, flagP, choice_int)
    else:
        flag=0
        plot(razb, ax, Pnts, flag, flagP, choice_int)

def checkbutton_changed2():
    global flagP
    if enabled2.get() == 1:
        flagP=1
        plot(razb, ax, Pnts, flag, flagP, choice_int)
    else:
        flagP=0
        plot(razb, ax, Pnts, flag, flagP, choice_int)

def display_selected(choice):
    global choice_int
    global choice_1
    global Methods
    choice_int=choice
    if choice_int=="Одномерный интеграл":
        #Methods = ["Метод Монте-Карло", "Метод Гаусса-2"]
        menu=combobox2['menu']
        menu.delete(0, 'end')
        combobox2.set_menu("None", "Метод Монте-Карло", "Метод Гаусса-3", "Метод Гаусса-4")
    elif choice_int=="Двумерный интеграл":
        menu = combobox2['menu']
        menu.delete(0, 'end')
        combobox2.set_menu("None", "Метод Монте-Карло", "Метод Гаусса-4", "Метод Гаусса-12")
    elif choice_int=="Трехмерный интеграл":
        menu = combobox2['menu']
        menu.delete(0, 'end')
        combobox2.set_menu("None", "Метод Монте-Карло", "Метод Гаусса-8", "Метод Гаусса-14")
    print(choice_int)

def display_selected_meth(choice):
    global choice_meth
    global choice_1
    choice_meth=choice
    if choice_meth!="Метод Монте-Карло":
        PointsNumber.place_forget()
        label2.place_forget()
    else:
        PointsNumber.place(relx=0.7, rely=0.5)
        label2.place(relx=0.55, rely=0.5)
    print(choice_meth)

def click_button():
    global Pnts
    global N
    global razb
    labelWarn.place_forget()
    N = int(PointsNumber.get())
    razb = int(DiscretNumber.get())
    checkBox_Points.place(relx=0.05, rely=0.95)
    checkBox_Points2.place(relx=0.4, rely=0.95)
    checkBox_Scale.place(relx=0.7, rely=0.95)
    checkBox_Legend.place(relx=0.85, rely=0.95)
    print(razb)
    if choice_int=="Одномерный интеграл" and choice_meth=="Метод Монте-Карло":
        checkBox_Points.place(relx=0.05, rely=0.95)
        checkBox_Points2.place(relx=0.4, rely=0.95)
        res, Pnts = kp.Monte_Karlo_1(N, razb)
        labelRes["text"] = str(res)
        plot(razb, ax, Pnts, flag, flagP, choice_int)
    elif choice_int=="Двумерный интеграл" and choice_meth=="Метод Монте-Карло":
        checkBox_Points.place(relx=0.05, rely=0.95)
        checkBox_Points2.place(relx=0.4, rely=0.95)
        if kp2.function(2, 2)==kp2.function(3, 3) and kp2.function(2, 2)==1:
            res, Pnts = kp.Monte_Karlo2(N, razb)
            labelRes["text"] = str(res)
            plot(razb, ax, Pnts, flag, flagP, choice_int)
        else:
            res, Pnts = kp.Monte_Karlo2_3D(N, razb)
            labelRes["text"]=str(res)
            plot(razb, ax, Pnts, flag, flagP, choice_int)
    elif choice_int=="Трехмерный интеграл" and choice_meth=="Метод Монте-Карло":
        checkBox_Points.place(relx=0.05, rely=0.95)
        checkBox_Points2.place(relx=0.4, rely=0.95)
        if kp3.function(2, 1 ,3)!=kp3.function(4, 3, 2) and kp3.function(2, 1,3)!=1:
            labelWarn['text'] = "Этот интеграл нельзя вычислить данным методом"
            labelWarn.place(relx=0.55, rely=0.2)
            checkBox_Points.place_forget()
            checkBox_Points2.place_forget()
            checkBox_Legend.place_forget()
            checkBox_Scale.place_forget()
            return
        res, Pnts = kp.Monte_Karlo2_3D_3(N, razb)
        labelRes["text"]=str(res)
        plot(razb, ax, Pnts, flag, flagP, choice_int)


    if choice_int=="Одномерный интеграл" and choice_meth=="Метод Гаусса-3":
        res = kp.Gausse_1(3)
        labelRes["text"]=str(res)
        checkBox_Points.place_forget()
        checkBox_Points2.place_forget()
        plot(razb, ax, Pnts, flag, flagP, choice_int)
    if choice_int == "Одномерный интеграл" and choice_meth == "Метод Гаусса-4":
        res = kp.Gausse_1(4)
        labelRes["text"] = str(res)
        checkBox_Points.place_forget()
        checkBox_Points2.place_forget()
        plot(razb, ax, Pnts, flag, flagP, choice_int)

    if choice_int=="Двумерный интеграл" and choice_meth=="Метод Гаусса-4":
        if kp2.tx(1) != kp2.tx(2) or kp2.bx(1)!=kp2.bx(2) or kp2.by(1)!=kp2.by(2) or kp2.ty(1)!=kp2.ty(2):
            labelWarn['text']="Этот интеграл нельзя вычислить данным методом"
            labelWarn.place(relx=0.55, rely=0.2)
            checkBox_Points.place_forget()
            checkBox_Points2.place_forget()
            checkBox_Legend.place_forget()
            checkBox_Scale.place_forget()
            return
        res = kp.Gausse_2(4)
        labelRes["text"] = str(res)
        checkBox_Points.place_forget()
        checkBox_Points2.place_forget()
        plot(razb, ax, Pnts, flag, flagP, choice_int)
    if choice_int=="Двумерный интеграл" and choice_meth=="Метод Гаусса-12":
        if kp2.tx(1) != kp2.tx(2) or kp2.bx(1)!=kp2.bx(2) or kp2.by(1)!=kp2.by(2) or kp2.ty(1)!=kp2.ty(2):
            labelWarn['text']="Этот интеграл нельзя вычислить данным методо"
            labelWarn.place(relx=0.55, rely=0.2)
            checkBox_Points.place_forget()
            checkBox_Points2.place_forget()
            checkBox_Legend.place_forget()
            checkBox_Scale.place_forget()
            return
        res = kp.Gausse_2(12)
        labelRes["text"] = str(res)
        checkBox_Points.place_forget()
        checkBox_Points2.place_forget()
        plot(razb, ax, Pnts, flag, flagP, choice_int)
    if choice_int=="Трехмерный интеграл" and choice_meth=="Метод Гаусса-8":
        if kp3.tx(1, 1) != kp3.tx(2, 2) or kp3.bx(1, 1)!=kp3.bx(2, 2) or kp3.by(1, 1)!=kp3.by(2, 2) or kp3.ty(1, 1)!=kp3.ty(2, 2) or kp3.tz(1, 1)!=kp3.tz(2, 2) or kp3.bz(1, 1)!=kp3.bz(2, 2):
            labelWarn['text']="Этот интеграл нельзя вычислить данным методо"
            checkBox_Points.place_forget()
            checkBox_Points2.place_forget()
            checkBox_Legend.place_forget()
            checkBox_Scale.place_forget()
            return
        res = kp.Gausse_3(8)
        labelRes["text"] = str(res)
        checkBox_Points.place_forget()
        checkBox_Points2.place_forget()
        if kp3.function(2, 1 ,3)!=kp3.function(4, 3, 2) and kp3.function(2, 1,3)!=1:
            labelWarn.place(relx=0.55, rely=0.2)
            labelWarn["text"]="У данного интеграла нет геометрической интерпретации"
            checkBox_Points.place_forget()
            checkBox_Points2.place_forget()
            checkBox_Legend.place_forget()
            checkBox_Scale.place_forget()
            return
        plot(razb, ax, Pnts, flag, flagP, choice_int)
    if choice_int=="Трехмерный интеграл" and choice_meth=="Метод Гаусса-14":
        res = kp.Gausse_3(14)
        labelRes["text"] = str(res)
        checkBox_Points.place_forget()
        checkBox_Points2.place_forget()
        if kp3.function(2, 1 ,3)!=kp3.function(4, 3, 2) and kp3.function(2, 1,3)!=1:
            labelWarn.place(relx=0.55, rely=0.2)
            labelWarn["text"] = "У данного интеграла нет геометрической интерпретации"
            checkBox_Points.place_forget()
            checkBox_Points2.place_forget()
            checkBox_Legend.place_forget()
            checkBox_Scale.place_forget()
            return
        plot(razb, ax, Pnts, flag, flagP, choice_int)

def checkbutton_Scale_changed():
    global sc
    if enabled1.get() == 1:
        sc=1
        plot(razb, ax, Pnts, flag, flagP, choice_int)
    else:
        sc=0
        plot(razb, ax, Pnts, flag, flagP, choice_int)

def checkbutton_Legend_changed():
    global legend
    if enabled_L.get()==1:
        legend=1
        plot(razb, ax, Pnts, flag, flagP, choice_int)
    else:
        legend=0
        plot(razb, ax, Pnts, flag, flagP, choice_int)

root: Tk = tk.Tk()
root.geometry("1000x600")
root.title('Курсовой проект')
root.resizable(0, 0)

frame = tk.Frame(root, borderwidth=1, relief="solid")
frame.pack(anchor="nw", padx=20, pady=10)
fig, ax = plt.subplots(figsize=(5, 5), dpi=100)
canvas = FigureCanvasTkAgg(fig, master=frame)
canvas.get_tk_widget().pack()

toolbar = NavigationToolbar2Tk(canvas, root)
toolbar.update()
canvas.get_tk_widget().pack()


enabled = tk.IntVar()
checkBox_Points=ttk.Checkbutton(text="Точки вне области", variable=enabled, command=checkbutton_changed, master=frame)
#checkBox_Points.place(relx=0.05, rely=0.95)

enabled2 = tk.IntVar()
checkBox_Points2=ttk.Checkbutton(text="Точки в области", variable=enabled2, command=checkbutton_changed2, master=frame)
#checkBox_Points2.place(relx=0.4, rely=0.95)

enabled1 = tk.IntVar()
checkBox_Scale=ttk.Checkbutton(text="1:1", variable=enabled1, command=checkbutton_Scale_changed, master=frame)
#checkBox_Scale.place(relx=0.7, rely=0.95)

enabled_L = tk.IntVar()
checkBox_Legend=ttk.Checkbutton(text="Легенда", variable=enabled_L, command=checkbutton_Legend_changed, master=frame)
#checkBox_Legend.place(relx=0.85, rely=0.95)

choice_=tk.StringVar()
choice_.set("None")
Types=["None", "Одномерный интеграл", "Двумерный интеграл", "Трехмерный интеграл"]
combobox = ttk.OptionMenu(root, choice_, *Types, command=display_selected)
combobox.place(relx=0.7, rely=0.02)

choice_1=tk.StringVar()
choice_1.set("None")
combobox2 = ttk.OptionMenu(root, choice_1, *Methods, command=display_selected_meth)
combobox2.place(relx=0.7, rely=0.07)

cmbLabel=ttk.Label(text="Тип интеграла")
cmbLabel.place(relx=0.61, rely=0.021)

cmbLabel=ttk.Label(text="Метод")
cmbLabel.place(relx=0.61, rely=0.071)

btn = ttk.Button(text="Вычислить интеграл", command=click_button, master=root)
btn.place(relx=0.7, rely=0.8)

PointsNumber = ttk.Entry()
PointsNumber.place(relx=0.7, rely=0.5)
PointsNumber.insert(0, "1000")

DiscretNumber = ttk.Entry()
DiscretNumber.place(relx=0.7, rely=0.4)
DiscretNumber.insert(0, "10")

label1 = ttk.Label(text="Дискретизация области")
label1.place(relx=0.55, rely=0.4)

label2 = ttk.Label(text="Кол-во точек")
label2.place(relx=0.55, rely=0.5)

labelRes = ttk.Label(font=("Arial", 16))
labelRes.place(relx=0.7, rely=0.9)

labelWarn = ttk.Label(text="У данного интеграла нет геометрической интерпретации", foreground='red')


#plot(1000, ax, Px, Py, flag, choice_int)
root.mainloop()