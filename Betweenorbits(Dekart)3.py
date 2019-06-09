# -*- coding: utf-8 -*-
"""
Created on Tue May 28 22:48:24 2019

@author: FOBOS
"""

import numpy as np
from scipy.integrate import odeint,ode
import matplotlib.pyplot as plt
import math as ms
alpha0=0
const1=398440000000000
r0=384405000 
w0=0.00000265
const2=4892100640000
rm=1738000
alphaM=0.91
l=1
phi=3.14/2+ms.asin(-0.1592161887268561)-w0*842899.9+0.17
print(phi)
#alphaM=phi
#xm=r0*ms.sin(alphaM)
#ym=r0*ms.cos(alphaM)


#r=np.sqrt(y**2+x**2)
#rm=np.sqrt((xm-x)**2+(ym-y)**2) ut=((ux*sin(phi)**2+uy*cos(phi)**2)**0.5

#cos(phi)=x/r, sin(phi)=y/r
 
def fout2(t, y):# обработчик шага 
        ts.append(t)
        YY.append(list(y.copy()))
        y1, y2, y3, y4,y5 = y
        if np.sqrt((r0*ms.sin(y5)-y1)**2+(r0*ms.cos(y5)-y2)**2) <= 5000000:
            return -1

def fout1(t, y):# обработчик шага 
        ts.append(t)
        YY.append(list(y.copy()))
#        y1, y2, y3, y4 = y
        

def fun0(t,y):
    y1,y2,y3,y4=y
    r=np.sqrt((y1**2)+(y2**2))
    return [y3,y4,((-1)*const1*y1)/((np.sqrt(y1**2+y2**2))**3)+
            (F*p*ms.cos(y2/r))/(mass-mu*t),
    ((-1)*const1*y2)/((np.sqrt(y1**2+y2**2))**3)+
    (F*p*ms.cos(y1/r))/(mass-mu*t)]
    
def fun1(t,y):
    y1,y2,y3,y4=y
    return [y3,y4,((-1)*const1*y1)/((np.sqrt(y1**2+y2**2))**3)+
            (F*p*ms.sin((alpha0+w*t)*3.14/180))/(mass-mu*t),
    ((-1)*const1*y2)/((np.sqrt(y1**2+y2**2))**3)+
    (F*p*ms.cos((alpha0+w*t)*3.14/180))/(mass-mu*t)]
    
    
def fun2(t,y):
    y1,y2,y3,y4,y5=y
    return [y3,y4,((-1)*const1*y1)/((np.sqrt(y1**2+y2**2))**3)+
            (F*p*ms.sin((alpha0+w*t)*3.14/180))/(mass-mu*t)-((const2*(y1-r0*ms.sin(y5)))
            /((np.sqrt((y1-r0*ms.sin(y5))**2+(y2-r0*ms.cos(y5))**2))**3)),
    ((-1)*const1*y2)/((np.sqrt(y1**2+y2**2))**3)+
    (F*p*ms.cos((alpha0+w*t)*3.14/180))/(mass-mu*t)-((const2*(y2-r0*ms.cos(y5)))
            /((np.sqrt((y1-r0*ms.sin(y5))**2+(y2-r0*ms.cos(y5))**2))**3)),w0]    
fuel=63.7164-0.5
h=281839.5
fuel=fuel*1000+17700
print(h)
print(str((fuel-17700)/1000)+'%')
mu0=246
ut0=4130
mu1=31.4
ut1=3050
mu=mu0
ut=ut0
F=mu*ut
p=1
y0=0
x0=-6375000-h
uy0=7750.22
ux0=0
#w=(-uy0/x0)*180/3.14
w=0.062
print(w)
mass=fuel+13000+25300
YY=[]    
tmax=(fuel-17700)/mu
print(tmax)
ts=[]
ys,xs=[],[]
Y0,t0=[x0,y0,ux0,uy0],0
ODE=ode(fun1)
ODE.set_integrator('dopri5',nsteps=10000)#, max_step=0.01)
ODE.set_solout(fout1)
ODE.set_initial_value(Y0, t0) # задание начальных значений
ODE.integrate(tmax)      # решение ОДУ
print()
Y=np.array(YY)
ys,xs=list(Y[:,1]),list(Y[:,0])
print(Y[-1,2])
print(Y[-1,3])
R=(Y[len(ts)-1][0]**2+Y[len(ts)-1][1]**2)**0.5
print(R)
print()
V=(Y[len(ts)-1][2]**2+Y[len(ts)-1][3]**2)**0.5
A=(R*const1)/((-(V**2)*R+2*const1))
print(A*2)
print(r0-x0)
print(3.14*np.sqrt(A**3/const1))
print(V)
#p2
mass=43000
fuel=17700
u=ut1
mu=mu1
F=u*mu
alpha0=w*tmax
print(alpha0)
x0,y0,ux0,uy0=list(Y[-1])

tmax=270000
#    ts=np.linspace(0,tmax,5000)

p=0
ym,xm=[],[]
#    T=odeint(fun10,[x0,y0,ux0,uy0],ts)
YY=[]
ts=[]
Y0,t0=[x0,y0,ux0,uy0,alphaM],0
ODE=ode(fun2)
ODE.set_integrator('dopri5')#, max_step=0.01)
ODE.set_solout(fout2)
ODE.set_initial_value(Y0, t0) # задание начальных значений
ODE.integrate(tmax)      # решение ОДУ
print()
Y=np.array(YY)

for i in range(len(ts)):
    ym.append(r0*ms.cos(Y[i][4]))
    xm.append(r0*ms.sin(Y[i][4]))

for i in range(len(ts)):
    ys.append(Y[i][1])
    xs.append(Y[i][0])

def circle(x,r):
    return np.array(np.sqrt(r**2-x**2))
xmm=np.linspace(-r0,r0,5000)
ymm=circle(xmm,r0)
plt.figure()
plt.plot(xs,ys)
#plt.plot(xmm,-ymm)
plt.plot(xm,ym)
plt.axis('equal')
plt.grid(True)
plt.show()

#const1=const2
x,y,ux,uy,alphaM=list(Y[-1])
#phi=y/(2*A)
#print(phi)
r=np.sqrt((x-r0*ms.sin(alphaM))**2+(y-r0*ms.cos(alphaM))**2)
V=(ux**2+uy**2)**0.5
P=(ux*(x-r0*ms.sin(alphaM))+uy*(y-r0*ms.cos(alphaM)))/(V*r)
ux=V*P
print(P)
print(r)
uy=l*V*np.sqrt(1-P**2)-w0*r0
print(ux,uy)
print(np.sqrt(ux**2+uy**2))
x0=r
y0=0
alpha0=0
const1=const2
YY=[]
ts=[]
tmax=2800
Y0,t0=[x0,y0,ux,uy],0
ODE=ode(fun1)
ODE.set_integrator('dopri5')#, max_step=0.01)
ODE.set_solout(fout1)
ODE.set_initial_value(Y0, t0) # задание начальных значений
ODE.integrate(tmax)      # решение ОДУ
print()
Y=np.array(YY)
x,y,ux,uy=list(Y[-1])
print(list(Y[-1]))
w=0
alpha0=90
p=1
ts=[]
tmax=311
Y0,t0=[x,y,ux,uy],0
ODE.set_initial_value(Y0, t0)
ODE.integrate(tmax)
Y=np.array(YY)
print(list(Y[-1]))
mass+=-mu*tmax
tmax=15000
p=0
x,y,ux,uy=list(Y[-1])
Y0,t0=[x,y,ux,uy],0
ODE.set_initial_value(Y0, t0)
ODE.integrate(tmax)
Y=np.array(YY)

xs1,ys1=list(Y[:,0]),(Y[:,1])
xb1=np.linspace(-rm,rm,500)
yb1=circle(xb1,rm)
x,y,ux,uy=list(Y[-1])
print(np.sqrt(x**2+y**2))
plt.figure()
plt.plot(xs1,ys1)
plt.plot(xb1,yb1)
plt.plot(xb1,-yb1)
plt.axis('equal')
plt.grid(True)
plt.show()
u2=np.sqrt(const2/1970248)
print(u2)





"""
T=odeint(fun0,[x0,y0,ux0,uy0],ts)
for i in range(len(ts)):
ys.append(T[i][1])
xs.append(T[i][0])
print(T[len(ts)-1])
print((T[len(ts)-1][0]**2+T[len(ts)-1][1]**2)**0.5)
"""



