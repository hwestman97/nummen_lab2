
#ATT FRÅGA
#Uppgift 2: vilken funktion
#Uppgift 3: var är matrisen?????

import numpy as np
from numpy import e, pi, log
import matplotlib.pyplot as plt
from math import exp,sqrt
import functions #byt namn!!!

G = 6.674*10**(-11) #m^3kg^(-1)s^(-2)
M = 5.974*10**(24) #kg
m = 7.348*10**(22) #kg
R = 3.844*10**8 #m
w = 2.662*10**(-6) #s^(-1)

def uppgift1():
    print(newtons_method(10,f,df,1000))

def f(r):
    f = G*M*R**2 - 2*G*M*R*r-w**2*r**3*R**2+2*w**2*R*r**4-w**2*r**5+G*M*r**2-G*m*r**2
    return f

def df(r):
    df = -2*G*M*R - 3*w**2*r**2*R**2 + 8*w**2*R*r**3-5*w**2*r**4+2*G*M*r-2*G*m*r
    return df

def newtons_method(x0, function, derivative, accuracy):
    fel = 100000
    n = 0
    while fel > accuracy:
        n += 1
        x1 = x0 -(function(x0))/(derivative(x0))
        fel = x1-x0
        x0 = x1
    print(n)
    return x0


def uppgift2():
    N = 1000
    u = []
    T = np.linspace(300,10000,N)
    for i in range(N):
        u.append(integral(T[i]))
    plt.plot(T,u)
    print('Maximum efficiency with numpy: ',round(T[np.argmax(u)],2), 'K')
    print('Maximum efficiency with golden ratio search: ',round(functions.golden_ratio_search(neg_integral,6000,7000,0.05),2), 'K')

def integral(T):
    n = 100
    h = 6.6261*10**(-34)
    c = 299792458
    lambda1 = 3.9*10**(-7)
    lambda2 = 7.5*10**(-7)
    kB = 1.38064852*10**(-23)
    a = (h*c)/(lambda2*kB*T)
    b = (h*c)/(lambda1*kB*T)
    u = (functions.quadratic_gauss(n,a,b,efficiency))
    return u

def neg_integral(T):
    u = integral(T)
    return -u

def efficiency(x):
    eta = (15/(pi**4))*(x**3/(e**x-1))
    return eta


def uppgift3():
    decay = np.array([log(50),log(33),log(31),log(14),log(20),log(13)])
    time = np.array([0.25,0.75,1.25,1.75,2.25,2.75])
    N = len(time)
    E_x = np.sum(time)/N
    E_y = np.sum(decay)/N
    E_xx = (np.sum(time**2))/N
    E_xy = (np.sum(decay*time))/N
    m = (E_xy-E_x*E_y)/(E_xx-E_x**2)
    c = (E_xx*E_y-E_x*E_xy)/(E_xx-E_x**2)
    A = e**c
    l = -m
    print('A =',A,'lambda =',l)
    x = np.linspace(0,3,100)
    y = m*x+c
    plt.plot(time,decay)
    plt.plot(x,y)
    A_m = np.matrix([[1,0.25],[1,0.75],[1,1.25],[1,1.75],[1,2.25],[1,2.75]])
    vector = np.dot(np.dot(np.linalg.inv(np.dot(A_m.transpose(),A_m)),A_m.transpose()),decay)
    print(vector)

uppgift3()

