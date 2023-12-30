import numpy as np

def function(x, y, z):
    return np.sin(x)*np.sin(y)

def tx(y, z):
    return 10
def bx(y, z):
    return 6

def ty(x, z):
    return 4
def by(x, z):
    return 0

def tz(x, y):
    return 2
def bz(x, y):
    return -2

def order():
    return['y', 'x', 'z']