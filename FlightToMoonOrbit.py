# -*- coding: utf-8 -*-
"""
Created on Tue May 28 22:48:24 2019

@author: FOBOS
"""

import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import math as ms
f=np.loadtxt("BlastOffFromEarth.txt",delimiter=' ', dtype=np.float)
out=open('FromEarthtoMoonOut.txt','w')
dm, W, h, t_pr, phi =f
print(phi)

pi=ms.pi
const1=398440000000000 #гравитационный параметр Земли
r0=384405000 #радиус орбиты Луны
w0=0.00000265 #угловая скорость Луны
const2=4892100640000 #гравитационный параметр Луны
rm=1738000 #радиус Луны


alphaM=0.91-0.2827 #начальный угол луны с осью Оу в новых координатах

t0=(pi-phi+w0*t_pr+(pi/2)-alphaM)/(W-w0) #время свободного вращения до попадания в точку с координатами (х0,у0) в новых координатах
teta=w0*t0+(pi/2)-alphaM+w0*t_pr #угол поворота ,используемой системы координат, относительной системы с осью х направленной к начальному положению луны
print(teta)
print(t0)
YYY=[]
#РС
mu0=246
ut0=4130
#ЛК
mu1=31.4
ut1=3050
mu=mu0
ut=ut0
F=mu*ut
p=1 #булево состояние двигателя
y0=0
x0=-h
uy0=W*h
ux0=0
fuel=dm*100
fuel=fuel*1000+17700-710
print(h)
print(str((fuel-17700)/1000)+'%')

mass=fuel+13000+25300
YY=[]    
ts=[]
ys,xs=[],[]

def fout2(t, y):# обработчик шага(перелет) 
        ts.append(t)
        YY.append(list(y.copy()))
        y1, y2, y3, y4,y5 = y
        if np.sqrt((r0*ms.sin(y5)-y1)**2+(r0*ms.cos(y5)-y2)**2) <= 4000000:
            return -1
#Обработчики(торможение)
def fout3(t,y):
        YY.append(list(y.copy()))
        y1, y2, y3, y4,y5 = y
        y_m=y2-r0*ms.cos(y5)
        x_m=y1-r0*ms.sin(y5)
        ux_m=y3-w0*r0*ms.cos(y5)
        uy_m=y4+w0*r0*ms.sin(y5)
        YYY.append([x_m,y_m,ux_m,uy_m])
        
def fout4(t,y):
        ts.append(t)
        YY.append(list(y.copy()))
        y1, y2, y3, y4,y5 = y
        y_m=y2-r0*ms.cos(y5)
        x_m=y1-r0*ms.sin(y5)
        ux_m=y3-w0*r0*ms.cos(y5)
        uy_m=y4+w0*r0*ms.sin(y5)
        r=np.sqrt(x_m**2+y_m**2)
        arr.append(r)
        YYY.append([x_m,y_m,ux_m,uy_m])

    
def fun(t,y): #Правая часть ODE в системы Земля-ракета-Луна
    y1,y2,y3,y4,y5=y
    u=np.sqrt(y4**2+y3**2)
    return [y3,y4,((-1)*const1*y1)/((np.sqrt(y1**2+y2**2))**3)+
            (F*p*(y3/u))/(mass-mu*t)-((const2*(y1-r0*ms.sin(y5)))
            /((np.sqrt((y1-r0*ms.sin(y5))**2+(y2-r0*ms.cos(y5))**2))**3)),
    ((-1)*const1*y2)/((np.sqrt(y1**2+y2**2))**3)+
    (F*p*(y4/u))/(mass-mu*t)-((const2*(y2-r0*ms.cos(y5)))
            /((np.sqrt((y1-r0*ms.sin(y5))**2+(y2-r0*ms.cos(y5))**2))**3)),w0]

#Разгон
t1=(fuel-17700)/mu
Y0,t0=[x0,y0,ux0,uy0,alphaM],0
ODE=ode(fun)
ODE.set_integrator('dopri5',nsteps=10000)#, max_step=0.01)
ODE.set_solout(fout2)
ODE.set_initial_value(Y0, t0) # задание начальных значений
ODE.integrate(t1)      # решение ОДУ
print()
Y=np.array(YY)
ys,xs=list(Y[:,1]),list(Y[:,0])
#Сброс разгонной ступени
mass=43000
fuel=17700
u=ut1
mu=mu1
F=95750
mu=F/ut1
x0,y0,ux0,uy0,alphaM=list(Y[-1])
#Свободный полет к луне
t2=450000
p=0
ym,xm=[],[]
YY=[]
ts=[]
Y0,t0=[x0,y0,ux0,uy0,alphaM],0
ODE=ode(fun)
ODE.set_integrator('dopri5')#, max_step=0.01)
ODE.set_solout(fout2)
ODE.set_initial_value(Y0, t0) # задание начальных значений
ODE.integrate(t2)      # решение ОДУ
t2=ts[-1]
print()
Y=np.array(YY)
#Визуализация первого этапа(до торможения)
for i in range(len(ts)):
    ym.append(r0*ms.cos(Y[i][4]))
    xm.append(r0*ms.sin(Y[i][4]))

for i in range(len(ts)):
    ys.append(Y[i][1])
    xs.append(Y[i][0])

def circle(x,r):
    return np.array(np.sqrt(r**2-x**2))
plt.figure()
plt.plot(xs,ys)
plt.plot(xm,ym)
plt.axis('equal')
plt.grid(True)
plt.show()

t3=1303 #Время до начала торможения
t4=356.1 #Время торможения
YY==[]
YYY=[]
arr=[]
#свободный полет до начала торможения
x,y,ux,uy,alphaM=list(Y[-1])       
Y0,t0=[x,y,ux,uy,alphaM],0
ODE.set_solout(fout3)
ODE.set_initial_value(Y0, t0) # задание начальных значений
ODE.integrate(t3)
#Торможение
x,y,ux,uy,alphaM=list(YY[-1])
Y0,t0=[x,y,ux,uy,alphaM],0
p=-1
ODE.set_solout(fout3)
ODE.set_initial_value(Y0, t0) # задание начальных значений
ODE.integrate(t4)
mass+=-mu*(t4)
#Свободный полет для достижения периселения орбиты
p=0
t5=12915
ts=[]
x,y,ux,uy,alphaM=list(YY[-1])
Y0,t0=[x,y,ux,uy,alphaM],0
ODE.set_solout(fout4)
ODE.set_initial_value(Y0, t0) # задание начальных значений
ODE.integrate(t5)
#торможение до первой космической
x,y,ux,uy,alphaM=list(YY[-1])
Y0,t0=[x,y,ux,uy,alphaM],0
p=-1
t6=2.01
ODE.set_solout(fout3)
ODE.set_initial_value(Y0, t0) # задание начальных значений
ODE.integrate(t6)
mass+=-mu*t6
#Свободный полет по круговой орбите
p=0
t7=8000
arr=[]
x,y,ux,uy,alphaM=list(YY[-1])
Y0,t0=[x,y,ux,uy,alphaM],0
ODE.set_solout(fout4)
ODE.set_initial_value(Y0, t0) # задание начальных значений
ODE.integrate(t7)
a=max(arr)
b=min(arr)

Y=np.array(YY)

Y1=np.array(YYY)    
xs1,ys1=list(Y1[:,0]),list(Y1[:,1])
xb1=np.linspace(-rm,rm,500)
yb1=circle(xb1,rm)

#Визуализация траектории относительно Луны
plt.figure()
plt.plot(xs1,ys1)
plt.plot(xb1,yb1,color='black')
plt.plot(xb1,-yb1,color='grey')
plt.axis('equal')
plt.grid(True)
plt.show()

x_m=r0*ms.sin(Y[-1][4])
y_m=r0*ms.cos(Y[-1][4])
x,y=Y[-1][0],Y[-1][1]
print(b-rm,a-rm)
print(Y[-1])
V=np.sqrt((Y[-1][2]-w0*r0*ms.cos(Y[-1][4]))**2+(Y[-1][3]+w0*r0*ms.sin(Y[-1][4]))**2)
print(V)
print(np.sqrt(const2*2/(b+a)))
t_spent=t0+t1+t2+t3+t4+t5+t6+t7+t_pr #time spent on stages 1 and 2
print(t_spent)
H=np.sqrt(xs1[-1]**2+ys1[-1]**2)-rm
print(H)
#Выходные данные в исходной системе координат
# В исходной системе координат ось Ох направлена к центру Луны до взлета с Земли, ось Оу соноправлена с тангенсальной скорость Луны до взлета ракеты
X_m=x_m*ms.cos(teta)+y_m*ms.sin(teta)
Y_m=(-x_m*ms.sin(teta)+y_m*ms.cos(teta))*(-1)
X=x*ms.cos(teta)+y*ms.sin(teta)
Y=(-x*ms.sin(teta)+y*ms.cos(teta))*(-1)
print(X_m,Y_m,X,Y)
F_spent=(t4+t6)*mu # spent Fuel
print(F_spent)
W=(-V/(rm+H)) #положительное направление вращения-направление вращения луны
out.write(str(X)+' '+str(Y)+' '+str(X_m)+' '+str(Y_m)+' '+str(W)+' '+str(F_spent)+' '+str(t_spent))
out.close()