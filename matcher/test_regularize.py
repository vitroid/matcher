#!/usr/bin/env python3

from math import acos, pi

import numpy as np
import itertools as it

x = np.fromstring("6.27809535  3.07847587 -1.8775237 ", sep=" ")
y = np.fromstring("2.84170353 -6.53829636 -1.56800074", sep=" ")
z = np.fromstring("2.55597332 -0.66422792  6.7812746 ", sep=" ")

def regularize(ex,ey,ez):
    error = np.dot(ex,ey)
    print(error)
    x_ort=ex-(error/2)*ey
    y_ort=ey-(error/2)*ex
    z_ort=np.cross(x_ort, y_ort)
    print(x_ort, y_ort, z_ort)
    x_new = 1./2*(3-np.dot(x_ort,x_ort))*x_ort
    y_new = 1./2*(3-np.dot(y_ort,y_ort))*y_ort
    z_new = 1./2*(3-np.dot(z_ort,z_ort))*z_ort
    return x_new, y_new, z_new

xL = np.dot(x,x)**0.5
yL = np.dot(y,y)**0.5
zL = np.dot(z,z)**0.5
alpha = acos(np.dot(y,z)/(yL*zL))*180/pi
beta  = acos(np.dot(z,x)/(zL*xL))*180/pi
gamma = acos(np.dot(x,y)/(xL*yL))*180/pi
print(xL)
print(yL)
print(zL)
print((xL+yL+zL)/3)
print(alpha)
print(beta)
print(gamma)
print((alpha+beta+gamma)/3)



x /= xL
y /= yL
z /= zL


x,y,z = regularize(x,y,z)
x,y,z = regularize(x,y,z)
print(x,y,z)
for v in (x, y, z):
    print(np.dot(v,v))
for v1, v2 in it.combinations((x, y, z), 2):
    print(np.dot(v1,v2))
