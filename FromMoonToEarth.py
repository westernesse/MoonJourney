import numpy as np
import matplotlib.pyplot as plt
import math as m

from scipy.integrate import ode

def fout(t, y):
    ts.append(t)
    ys.append(list(y.copy()))
    if y[0]>(R_Orbit+10**8):
        return -1
    

def func(t, y):
    
    y1, y2, y3, y4, y5 = y
    y1m, y3m = y1 - m.cos(y5)*R_Orbit, y3 - m.sin(y5)*R_Orbit
    
    ax = -1 * G * M_Earth * y1 * (y1*y1 + y3*y3)**(-1.5) - G * M_Moon * (y1m) * (y1m**2 + y3m**2)**(-1.5)
    ay = -1 * G * M_Earth * y3 * (y1*y1 + y3*y3)**(-1.5) - G * M_Moon * (y3m) * (y1m**2 + y3m**2)**(-1.5)

    return [y2, ax, y4, ay, angular_moon_speed]



#НАЧАЛЬНЫЕ УСЛОВИЯ
    
    

G = 6.6743 * 10**(-11)
M_Earth = 5.972 * 10**(24)
R_Earth = 6371000
R_Orbit = 385*10**6
M_Moon = 7.348*10**22

angular_moon_speed = 1023/R_Orbit

Moon_angle_start = 0
x_start = R_Orbit - 2000000
y_start = 2000000
Vx_start = - 1400
Vy_start = -700

time = 700000 #время



xc,yc=[],[]                                #Табличная функция для круга Земли
for i in range(0, 630):
    xc.append(R_Earth*m.cos(i/100))
    yc.append(R_Earth*m.sin(i/100))
    

    
    


y0,t0=[x_start,  Vx_start, y_start, Vy_start, Moon_angle_start], 0 # начальные условия для интегрирования

ODE = ode(func)
ODE.set_integrator('dopri5')#, max_step=0.01)
ODE.set_solout(fout)
ts, ys = [ ],[ ]
ODE.set_initial_value(y0, t0) # задание начальных значений
ODE.integrate(time)      # решение ОДУ
Y=np.array(ys)
period = len(Y)//14



#Вывод графиков для разных случаев

for k in range(1, len(Y), period):
    if (Y[k, 0]**2 + Y[k, 2]**2)**0.5 < 30000000:
        Y1 = Y[50:k]
        plt.plot(Y1[:,0],Y1[:,2],linewidth=2)#,label='k=%.1f'% k)
        plt.axis('equal')
        plt.plot(xc,yc,linewidth=1)
        xt,yt = [], []
        for i in range(0, 630):
            xt.append(200000*m.cos(i/100) + Y1[-1, 0])
            yt.append(200000*m.sin(i/100) + Y1[-1, 2])
        xj, yj = [], []
        for i in range(100):
            xj.append(Y1[-1,0]+600000*(i)/100 + 200000)
            yj.append(Y1[-1,2]+600000*((i)/100)**4 + 200000)
            xj.append(Y1[-1,0]+600000*(i)/100 + 200000)
            yj.append(Y1[-1,2]+600000*(i/100)**0.25 + 200000)
        plt.plot(xt,yt,'r', linewidth=6)
        plt.plot(xj,yj,'g', linewidth=2)
        plt.title("Flight near Earth\n "+str(k//period))
        
        plt.grid(True)
        
        plt.show()
    elif ((Y[k, 0])**2 + (Y[k, 2])**2)**0.5 > 300000000:
        if k < 20:
            Y1 = Y[0:k]
        else:
            Y1 = Y[k-20:k]
        Y2 = []
        Y3 = []
        for i in range(len(Y1)):
            Y2.append(m.cos(Y1[i,4])*R_Orbit)
            Y3.append(m.sin(Y1[i,4])*R_Orbit)
        plt.plot(Y2[:], Y3[:], 'b', linewidth = 1)
        plt.plot(Y1[:,0],Y1[:,2],linewidth=2)#,label='k=%.1f'% k)
        plt.axis('equal')
        Moon_x, Moon_y = m.cos(Y1[-1,4])*R_Orbit, m.sin(Y1[-1,4])*R_Orbit
        xn,yn = [], []
        for i in range(0, 630):
            xn.append(1737000*m.cos(i/100) + Moon_x)
            yn.append(1737000*m.sin(i/100) + Moon_y)
        plt.plot(xn,yn,linewidth=2)
        xt,yt = [], []
        for i in range(0, 630):
            xt.append(200000*m.cos(i/100) + Y1[-1, 0])
            yt.append(200000*m.sin(i/100) + Y1[-1, 2])
        xj, yj = [], []
        for i in range(100):
            xj.append(Y1[-1,0]+600000*(i)/100 + 200000)
            yj.append(Y1[-1,2]+600000*((i)/100)**4 + 200000)
            xj.append(Y1[-1,0]+600000*(i)/100 + 200000)
            yj.append(Y1[-1,2]+600000*(i/100)**0.25 + 200000)
        plt.plot(xt,yt,'r', linewidth=6)
        plt.plot(xj,yj,'g', linewidth=2)
        plt.title("Flight near Moon\n "+str(k//period))
        
        plt.grid(True)
        
        plt.show()
    else:
        Y1 = Y[:k]
        plt.plot(Y1[:,0],Y1[:,2],linewidth=2)
        Y2 = []
        Y3 = []
        for i in range(len(Y1)):
            Y2.append(m.cos(Y1[i,4])*R_Orbit)
            Y3.append(m.sin(Y1[i,4])*R_Orbit)
        plt.plot(Y2[:], Y3[:], 'b', linewidth = 1)
        plt.axis('equal')
        plt.plot(xc,yc,linewidth=1)
        Moon_x, Moon_y = m.cos(Y1[-1,4])*R_Orbit, m.sin(Y1[-1,4])*R_Orbit
        xn,yn = [], []
        for i in range(0, 630):
            xn.append(1737000*m.cos(i/100) + Moon_x)
            yn.append(1737000*m.sin(i/100) + Moon_y)
        plt.plot(xn,yn,linewidth=2)
        xt,yt = [], []
        for i in range(0, 630):
            xt.append(4000000*m.cos(i/100) + Y1[-1, 0])
            yt.append(4000000*m.sin(i/100) + Y1[-1, 2])
        xj, yj = [], []
        for i in range(100):
            xj.append(Y1[-1,0]+12000000*(i)/100 + 200000)
            yj.append(Y1[-1,2]+12000000*((i)/100)**4 + 200000)
            xj.append(Y1[-1,0]+12000000*(i)/100 + 200000)
            yj.append(Y1[-1,2]+12000000*(i/100)**0.25 + 200000)
        plt.plot(xt,yt,'r', linewidth=6)
        plt.plot(xj,yj,'g', linewidth=2)
        plt.title("General flight\n "+str(k//period))
        
        plt.grid(True)
        
        plt.show()