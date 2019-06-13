import numpy as np
import matplotlib.pyplot as plt
import math as m

from scipy.integrate import ode


with open('Escaping_from_moon.txt') as f:
    angle_start_0 = float(f.readline())
    w = float(f.readline())
    h = float(f.readline())
    T = float(f.readline())



def graph(Y, period):                                            #Функция вывода графиков
    for k in range(1, len(Y), period):                       #Период - шаг, с которым проходятся элементы списка Y с проинтегрированными значениями
        if (Y[k, 0]**2 + Y[k, 2]**2)**0.5 < 70000000:               #Одно из трех условий для вывода разных видов графиков
            
            Y1 = Y[30000:k]
            plt.plot(Y1[:,0],Y1[:,2],linewidth=2)#,label='k=%.1f'% k)
            plt.axis('equal')
            plt.plot(xc,yc,linewidth=1)
            
            xt,yt = [], []                                              #Рисую плод яблока
            for i in range(0, 630):                 
                xt.append(200000*m.cos(i/100) + Y1[-1, 0])
                yt.append(200000*m.sin(i/100) + Y1[-1, 2])
            
            xj, yj = [], []                                             #Рисую листик яблока
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
        
        elif ((Y[k, 0])**2 + (Y[k, 2])**2)**0.5 > 330000000:
            if k < 50000:
                Y1 = Y[0:k]
            else:
                Y1 = Y[k-50000:k]
            
            Y2 = []
            Y3 = []
            
            for i in range(len(Y1)):
                Y2.append(m.cos(Y1[i,4])*R_Orbit)
                Y3.append(m.sin(Y1[i,4])*R_Orbit)
            plt.plot(Y2[:], Y3[:], 'b', linewidth = 1)
            plt.plot(Y1[:,0],Y1[:,2],linewidth=2)
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




def fout1(t, y):
    ts.append(t)
    ys.append(list(y.copy()))
    if Force*t/Gas_velocity >= m1:
        return -1
    
def fout2(t, y):
    ts.append(t)
    ys.append(list(y.copy()))
    if y[0]<0 and y[2]<0 or (y[0]*y[0]+y[2]*y[2])**0.5 < 1.01*R_Earth:        #
        return -1
    
def fout3(t, y):
    ts.append(t)
    ys.append(list(y.copy()))
    
    if Force*t/Gas_velocity + m1 >= m_fuel:
        return -1
    
    
def fout4(t, y):
    ts.append(t)
    ys.append(list(y.copy()))
    if (y[0]*y[0]+y[2]*y[2])**0.5 < R_Earth + 71000:
        return -1
    

def func1(t, y):
    
    
    
    y1, y2, y3, y4, y5 = y
    y1m, y3m = y1 - m.cos(y5)*R_Orbit, y3 - m.sin(y5)*R_Orbit
    
    dm_fuel = Force*t/Gas_velocity
    a = Force/(Total_Mass - dm_fuel)
    angle = angle_start
    
    
    
    ax = -1 * G * M_Earth * y1 * (y1*y1 + y3*y3)**(-1.5) - G * M_Moon * (y1m) * (y1m**2 + y3m**2)**(-1.5) - a*m.cos(angle + 0.5*m.pi)
    ay = -1 * G * M_Earth * y3 * (y1*y1 + y3*y3)**(-1.5) - G * M_Moon * (y3m) * (y1m**2 + y3m**2)**(-1.5) - a*m.sin(angle + 0.5*m.pi)

    return [y2, ax, y4, ay, angular_moon_speed]



def func2(t, y):
    
    y1, y2, y3, y4, y5 = y
    y1m, y3m = y1 - m.cos(y5)*R_Orbit, y3 - m.sin(y5)*R_Orbit
    

    
    ax = -1 * G * M_Earth * y1 * (y1*y1 + y3*y3)**(-1.5) - G * M_Moon * (y1m) * (y1m**2 + y3m**2)**(-1.5)
    ay = -1 * G * M_Earth * y3 * (y1*y1 + y3*y3)**(-1.5) - G * M_Moon * (y3m) * (y1m**2 + y3m**2)**(-1.5)

    return [y2, ax, y4, ay, angular_moon_speed]



def func3(t, y):
    global Hmin, V, Time
    y1, y2, y3, y4, y5 = y
    y1m, y3m = y1 - m.cos(y5)*R_Orbit, y3 - m.sin(y5)*R_Orbit
    
    if (y1*y1 + y3*y3)**0.5 - R_Earth < Hmin + 1 and (y1*y1 + y3*y3)**0.5 - R_Earth > 70000:
        Hmin = (y1*y1 + y3*y3)**0.5 - R_Earth
        V = (y2*y2 + y4*y4)**0.5
        Time = y5/(2*m.pi)*27.32*24*60*60
    
    dm_fuel = Force*t/Gas_velocity + m1
    a = - Force/(Total_Mass - dm_fuel)
    angle = m.pi#*135/180
    
    ax = -1 * G * M_Earth * y1 * (y1*y1 + y3*y3)**(-1.5) - G * M_Moon * (y1m) * (y1m**2 + y3m**2)**(-1.5) + a*m.cos(angle + 0.5*m.pi)
    ay = -1 * G * M_Earth * y3 * (y1*y1 + y3*y3)**(-1.5) - G * M_Moon * (y3m) * (y1m**2 + y3m**2)**(-1.5) + a*m.sin(angle + 0.5*m.pi)

    return [y2, ax, y4, ay, angular_moon_speed]


def func4(t, y):
    global Hmin, V, Time
    
    y1, y2, y3, y4, y5 = y
    y1m, y3m = y1 - m.cos(y5)*R_Orbit, y3 - m.sin(y5)*R_Orbit
    
    if (y1*y1 + y3*y3)**0.5 - R_Earth <= Hmin + 1 and (y1*y1 + y3*y3)**0.5 - R_Earth > 70000:
        Hmin = (y1*y1 + y3*y3)**0.5 - R_Earth
        V = (y2*y2 + y4*y4)**0.5
        Time = y5/(2*m.pi)*27.32*24*60*60
    
    ax = -1 * G * M_Earth * y1 * (y1*y1 + y3*y3)**(-1.5) - G * M_Moon * (y1m) * (y1m**2 + y3m**2)**(-1.5)
    ay = -1 * G * M_Earth * y3 * (y1*y1 + y3*y3)**(-1.5) - G * M_Moon * (y3m) * (y1m**2 + y3m**2)**(-1.5)

    return [y2, ax, y4, ay, angular_moon_speed]


def funcV(t, y):
    
    dm_fuel = Force*t/Gas_velocity
    a = Force/(Total_Mass - dm_fuel)
    
    return a


#НАЧАЛЬНЫЕ УСЛОВИЯ
    
v = 0  

G = 6.6743 * 10**(-11)
M_Earth = 5.972 * 10**(24)
R_Earth = 6371000
R_Orbit = 385*10**6
R_Moon = 173800
M_Moon = 7.348*10**22
m_dry = 22500 - 17700 + 5500
m_fuel = 17700 - 11242
#m_fuel = 17700
Total_Mass = m_dry + m_fuel

angular_moon_speed = 1023/R_Orbit

h = h - R_Moon

Moon_angle_start = 0
V_On_Orbit = (G*M_Moon/(h + R_Moon))**0.5
angle_start = 28.72*m.pi/180            #28.8+-
x_start = R_Orbit + m.cos(angle_start)*(R_Moon + h)
y_start = m.sin(angle_start)*(R_Moon + h)
Vx_start = -V_On_Orbit*m.cos(angle_start + 0.5*m.pi)
Vy_start = 1023 - V_On_Orbit*m.sin(angle_start + 0.5*m.pi)

Force = 95750
Gas_velocity = 3050

Hmin = R_Orbit


time1 = 800000
time2 = 90000 #время
m1 = 4047.5  #4047.5  


Time_waiting = (h+R_Moon)/V_On_Orbit * (m.pi - angle_start_0 - angle_start)


xc,yc=[],[]                                #Табличная функция для круга Земли
for i in range(0, 630):
    xc.append(R_Earth*m.cos(i/100))
    yc.append(R_Earth*m.sin(i/100))
    

    
#этап 1, ускорение    


y0,t0=[x_start,  Vx_start, y_start, Vy_start, Moon_angle_start], 0 # начальные условия для интегрирования
print(y0)
ODE = ode(func1)
ODE.set_integrator('dopri5', nsteps=200000, max_step = 0.5)
ODE.set_solout(fout1)
ts, ys = [ ],[ ]
ODE.set_initial_value(y0, t0) # задание начальных значений
ODE.integrate(time1)      # решение ОДУ
Y=np.array(ys)
period = len(Y)//15 + 1

#graph(Y, period)        #Вывод графиков
print(len(Y))
#этап 2, свободный полет

y0,t0=[Y[-1:,0],Y[-1:,1],Y[-1:,2],Y[-1:,3],Y[-1:,4]], 0 # начальные условия для интегрирования

ODE = ode(func2)
ODE.set_integrator('dopri5', nsteps=200000, max_step = 10)
ODE.set_solout(fout2)
ODE.set_initial_value(y0, t0) # задание начальных значений
ODE.integrate(time1)      # решение ОДУ
Y=np.array(ys)
period = len(Y)//20

#graph(Y, period)      #Вывод графиков 
print(len(Y))
#Этап 3, торможение

y0,t0=[Y[-1:,0],Y[-1:,1],Y[-1:,2],Y[-1:,3],Y[-1:,4]], 0

ODE = ode(func3)
ODE.set_integrator('dopri5', nsteps=200000, max_step = 1)
ODE.set_solout(fout3)
ODE.set_initial_value(y0, t0) # задание начальных значений
ODE.integrate(time1)      # решение ОДУ
Y=np.array(ys)
period = len(Y)//5 + 1

print(len(Y))

y0,t0=[Y[-1:,0],Y[-1:,1],Y[-1:,2],Y[-1:,3],Y[-1:,4]], 0 # начальные условия для интегрирования

ODE = ode(func4)
ODE.set_integrator('dopri5', nsteps=200000, max_step = 2)
ODE.set_solout(fout4)
ODE.set_initial_value(y0, t0) # задание начальных значений
ODE.integrate(time2)      # решение ОДУ
Y=np.array(ys)
period = len(Y)//10   #Период - шаг, с которым проходятся элементы списка Y с проинтегрированными значениями
#Чем меньше период, тем больше графиков выведется, если период равен len(Y), будет один график
graph(Y, period) 
print(len(Y))

#расчет скорости отлета

y0,t0=V_On_Orbit, 0 # начальные условия для интегрирования
ODE = ode(funcV)
ODE.set_integrator('dopri5', nsteps=200000, max_step = 0.5)
ODE.set_solout(fout1)
ts, ys = [ ],[ ]
ODE.set_initial_value(y0, t0) # задание начальных значений
ODE.integrate(time1)
Y=np.array(ys)      # решение ОДУ
v = Y[-1]
print('Скорость отлета = ', v)
print('Hmin (km) = ', Hmin/1000)
print('V = ', V, Time+Time_waiting)

f = open('Approaching_Home.txt', 'w')
f.write(str(Hmin)+'\n')
f.write(str(V)+'\n')
f.write(str(T + Time + Time_waiting))
f.close()