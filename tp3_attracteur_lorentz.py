#!/usr/bin/python
#-*- coding: utf-8 -*-
 
from math import *
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
# import des modules et fonctions nécessaire

t=60.
epsilon=1e-20
#epsilon=input('epsilon=')
xa,ya,za=epsilon,epsilon,epsilon
#xa,ya,za=10.0001,10.0001,20.0001
P,R,b=10.,28.,8./3.
a,t,n=0.,60.,100000
h=(t-a)/float(n)
#h=0.001
print('a, t, n =',a,t,n)
print('P, R, b =',P,R,b)
print('xa, ya, za =',xa,ya,za)
print('h=',h)
# initialisation des variables avec h le pas


def f1(x,y,t):
    return P*(y-x)
# définition de x'=f1

def f2(x,y,z,t):
    return R*x-x*z-y
# y'=f2
    
def f3(x,y,z,t):
    return x*y-b*z
# z'=f3

tt = np.empty(n+1, dtype=float)
xx = np.empty(n+1, dtype=float)
yy = np.empty(n+1, dtype=float)
zz = np.empty(n+1, dtype=float)

def euler(a,b,n,xa,ya,za,xx,yy,zz,tt):
    # xx, yy, zz sont des tableaux contenant dans l'ordre les valeurs de x,y et z(t)
    tt[0]=a
    xx[0]=xa
    yy[0]=ya
    zz[0]=za
    for i in range(1,n+1):
        tt[i]=tt[i-1]+h
        xx[i]=xx[i-1]+h*f1(xx[i-1],yy[i-1],tt[i-1])
        yy[i]=yy[i-1]+h*f2(xx[i-1],yy[i-1],zz[i-1],tt[i-1])
        zz[i]=zz[i-1]+h*f3(xx[i-1],yy[i-1],zz[i-1],tt[i-1])
    
    plt.figure('Euler')   
    plt.title('Attracteur de Lorentz par la methode de Euler')
    plt.plot(xx,yy,'b',label='y(x)')
    plt.plot(yy,zz,'g',label='z(y)')
    plt.plot(xx,zz,'r',label='z(x)')
    #plt.plot(xx,tt,'--c')
    #plt.plot(yy,tt,'--m')
    #plt.plot(zz,tt,'--y')
    plt.legend()
# algorithme de résolution d'euler



def rk(a,b,n,xa,ya,za,xx,yy,zz,tt):
    tt[0]=a
    xx[0]=xa
    yy[0]=ya
    zz[0]=za
    for i in range(1,n+1):
        tt[i]=tt[i-1]+h
          
        k11=(h/6.)*f1(xx[i-1],yy[i-1],tt[i-1])
        k12=(h/3.)*f1(xx[i-1]+k11/2.,yy[i-1]+k11/2.,tt[i-1]+h/2.)
        k13=(h/3.)*f1(xx[i-1]+k12/2.,yy[i-1]+k12/2.,tt[i-1]+h/2.)
        k14=(h/6.)*f1(xx[i-1]+k13,yy[i-1]+k13,tt[i-1]+h)
        xx[i]=xx[i-1]+(k11+k12+k13+k14)
         
        k21=(h/6.)*f2(xx[i-1],yy[i-1],zz[i-1],tt[i-1])
        k22=(h/3.)*f2(xx[i-1]+k21/2.,yy[i-1]+k21/2.,zz[i-1]+k21/2.,tt[i-1]+h/2.)
        k23=(h/3.)*f2(xx[i-1]+k22/2.,yy[i-1]+k22/2.,zz[i-1]+k22/2.,tt[i-1]+h/2.)
        k24=(h/6.)*f2(xx[i-1]+k23,yy[i-1]+k23,zz[i-1]+k23,tt[i-1]+h)
        yy[i]=yy[i-1]+(k21+k22+k23+k24)
        
        k31=(h/6.)*f3(xx[i-1],yy[i-1],zz[i-1],tt[i-1])
        k32=(h/3.)*f3(xx[i-1]+k31/2.,yy[i-1]+k31/2.,zz[i-1]+k31/2.,tt[i-1]+h/2.)
        k33=(h/3.)*f3(xx[i-1]+k32/2.,yy[i-1]+k32/2.,zz[i-1]+k32/2.,tt[i-1]+h/2.)
        k34=(h/6.)*f3(xx[i-1]+k33,yy[i-1]+k33,zz[i-1]+k33,tt[i-1]+h)
        zz[i]=zz[i-1]+(k31+k32+k33+k34)
    
    plt.figure('RK')
    plt.title('Attracteur de Lorentz par la methode RK')
    plt.plot(xx,yy,'b',label='y(x)')
    plt.plot(yy,zz,'g',label='z(y)')
    plt.plot(xx,zz,'r',label='z(x)')
    #plt.plot(xx,tt,'--c')
    #plt.plot(yy,tt,'--m')
    #plt.plot(zz,tt,'--y')
    plt.legend()
# algorithme de résolution de rk

euler(a,b,n,xa,ya,za,xx,yy,zz,tt)
rk(a,b,n,xa,ya,za,xx,yy,zz,tt)

plt.show()