import tkinter as tk
import matplotlib.pyplot as plt
import numpy as np
import KP_function1 as kp1
import KP_function2 as kp2
import KP_function3 as kp3
import KP as kp
from tkinter import ttk
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
NavigationToolbar2Tk)

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

def Plotting_3D1(ax, razb, tx, bx, ty, by, tz, bz, sc):
    x = np.linspace(bx(0, 0), tx(0, 0), int(razb))
    if kp3.order()[0]=='z':
        y = np.linspace(np.min(by(0, x)), np.max(ty(0, x)), int(razb))
    else:
        y = np.linspace(np.min(by(x, 0)), np.max(ty(x, 0)), int(razb))
    X1, Y1 = np.meshgrid(x, y)
    if kp3.order()[0]=='y':
        z = np.linspace(np.min(bz(Y1, X1)), np.max(tz(Y1, X1)), int(razb))
    else:
        z = np.linspace(np.min(bz(X1, Y1)), np.max(tz(X1, Y1)), int(razb))
    y1 = ty(x, 0)
    Y2, Z1 = np.meshgrid(x, z)
    #ax.plot_wireframe(Y, kp3.tx(Y2), Z1, color="olive", label="Ограничение по x")
    #ax.plot_wireframe(Y, kp3.bx(Y2), Z1, color="olive")
    X2, Z2 = np.meshgrid(y, z)
    maxx1 = np.max(x)
    minx1 = np.min(x)
    maxy1 = np.max(y1)
    miny1 = np.min(y1)
    #ax.plot_wireframe(kp3.ty(X), X, Z, color="orange", label="Ограничение по y")
    y2 = by(x, 0)
    #ax.plot_wireframe(kp3.by(X), X, Z, color="orange")
    #ax.plot_wireframe(X, Y, kp.BorderZ(X, Y, kp3.tz), color="green", label="Ограничение по z")
    #ax.plot_wireframe(X, Y, kp.BorderZ(X, Y, kp3.bz), color="green")
    maxx2 = np.max(x)
    minx2 = np.min(x)
    maxy2 = np.max(y2)
    miny2 = np.min(y2)

    #Scale(ax, max(maxx1, maxx2), min(minx1, minx2), max(maxy1, maxy2), min(miny1, miny2), sc)
    return X1, X2, Y1, Y2, Z1, Z2
'''''
def plot(fig, razb, ax, Pnts, flag, flagP, choice_int):
    ax.clear()
    fig.clear()
    print(choice_int)
    if choice_int=="Одномерный интеграл":
        
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
                if kp2.tx(2)==kp2.tx(2):
                    x = np.full(shape=razb, fill_value=kp2.tx(2))
                else:
                    x = kp2.tx(y)
                maxx1=np.max(x); minx1=np.min(x); maxy1=np.max(y); miny1=np.min(y)
                ax.plot(x, y, color="red", label="Ограничения по x")
                if kp2.bx(2)==kp2.bx(2):
                    x = np.full(shape=razb, fill_value=kp2.bx(2))
                else:
                    x = kp2.bx(y)
                ax.plot(x, y, color="red")
                maxx2 = np.max(x); minx2 = np.min(x); maxy2 = np.max(y); miny2 = np.min(y)
                #Scale(ax, max(maxx1, maxx2), min(minx1, minx2), max(maxy1, maxy2), min(miny1, miny2), sc)
            elif kp2.tx(1)==kp2.tx(2):
                ax.axhline(y=0, color='k', linewidth=0.5)
                ax.axvline(x=0, color='k', linewidth=0.5)
                ax.axvline(x=kp2.tx(0), color="red", label="Ограничения по x")
                ax.axvline(x=kp2.bx(0), color="red")
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
                #Scale(ax, max(maxx1, maxx2), min(minx1, minx2), max(maxy1, maxy2), min(miny1, miny2), sc)
        else:
            if kp2.ty(0)==kp2.ty(2):
                y = np.linspace(kp2.by(0), kp2.ty(0), int(razb/100))
                x = np.linspace(np.min(kp2.bx(y)), np.max(kp2.tx(y)), int(razb/100))
                X, Y = np.meshgrid(x, y)
                ax = fig.add_subplot(111, projection='3d')
                ax.plot_wireframe(X, Y, np.zeros((len(x), len(y))), color="black")
                ax.plot_wireframe(X, Y, kp2.function(X, Y), color="blue", label="Подинтегральная функция")
                z = np.linspace(np.min(kp2.function(X,Y)), np.max(kp2.function(X,Y)), int(razb / 100))
                x1 = kp2.tx(y)
                Y, Z = np.meshgrid(y, z)
                maxx1 = np.max(x1); minx1 = np.min(x1); maxy1 = np.max(y); miny1 = np.min(y)
                ax.plot_wireframe(kp2.tx(Y), Y, Z, color="red", label="Ограничение по x")
                x2 = kp2.bx(y)
                ax.plot_wireframe(kp2.bx(Y), Y, Z, color="red")
                maxx2 = np.max(x2); minx2 = np.min(x2); maxy2 = np.max(y); miny2 = np.min(y)
                #Scale(ax, max(maxx1, maxx2), min(minx1, minx2), max(maxy1, maxy2), min(miny1, miny2), sc)
            elif kp2.tx(2)==kp2.tx(3):
                x = np.linspace(kp2.bx(0), kp2.tx(0), int(razb / 100))
                y = np.linspace(np.min(kp2.by(x)), np.max(kp2.ty(x)), int(razb / 100))
                X, Y = np.meshgrid(x, y)
                ax = fig.add_subplot(111, projection='3d')
                ax.plot_wireframe(X, Y, np.zeros((len(x), len(y))), color="black")
                ax.plot_wireframe(X, Y, kp2.function(X, Y), color="blue", label="Подинтегральная функция")
                z = np.linspace(np.min(kp2.function(X, Y)), np.max(kp2.function(X, Y)), int(razb / 100))
                y1 = kp2.ty(x)
                X, Z = np.meshgrid(y, z)
                maxx1 = np.max(x); minx1 = np.min(x); maxy1 = np.max(y1); miny1 = np.min(y1)
                ax.plot_wireframe(X, kp2.ty(X), Z, color="orange", label="Ограничение по y")
                y2 = kp2.by(x)
                ax.plot_wireframe(X, kp2.by(X), Z, color="orange")
                maxx2 = np.max(x); minx2 = np.min(x); maxy2 = np.max(y2); miny2 = np.min(y2)
                #Scale(ax, max(maxx1, maxx2), min(minx1, minx2), max(maxy1, maxy2), min(miny1, miny2), sc)
    elif choice_int=="Трехмерный интеграл":
        ax = fig.add_subplot(111)
        if kp3.function(1, 1, 1)==kp3.function(2, 2, 2):
            x = np.linspace(kp3.bx(0), kp3.tx(0), int(razb / 100))
            y = np.linspace(np.min(kp3.by(x)), np.max(kp3.ty(x)), int(razb / 100))
            X, Y = np.meshgrid(x, y)
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_wireframe(X, Y, np.zeros((len(x), len(y))), color="black")
            z = np.linspace(np.min(kp3.bz(X, Y)), np.max(kp3.tz(X, Y)), int(razb / 100))
            y1 = kp3.ty(x)
            X, Z = np.meshgrid(y, z)
            maxx1 = np.max(x); minx1 = np.min(x); maxy1 = np.max(y1); miny1 = np.min(y1)
            ax.plot_wireframe(kp3.ty(X), X, Z, color="blue")
            y2 = kp3.by(x)
            ax.plot_wireframe(kp3.by(X), X, Z, color="orange")
            X, Y = np.meshgrid(x, y)
            ax.plot_wireframe(X, Y, kp.BorderZ(X, Y, kp3.tz))
            ax.plot_wireframe(X, Y, kp.BorderZ(X, Y, kp3.bz))
            maxx2 = np.max(x); minx2 = np.min(x); maxy2 = np.max(y2); miny2 = np.min(y2)
            #Scale(ax, max(maxx1, maxx2), min(minx1, minx2), max(maxy1, maxy2), min(miny1, miny2), sc)

'''''