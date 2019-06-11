import numpy as np
import matplotlib.pyplot as plt
import math as m
import random

from scipy.integrate import ode

arr = np.loadtxt("before STOP.txt", delimiter='\t', dtype=np.float)
h = arr[0] #70000 высота над поверхностью земли из файла
vx0 = arr[1] # (mu/(h+R_Earth))**0.5 скорость по орбите
vy0 = arr[2] #0
#print(h,vx0,vy0)

max_overload = 0

G=6.67e-11
totalTime=0
overload = []
g_Earth = 9.81
R_Earth = 6375*10**3 
S = 0.25*m.pi*3.9**2  
p0 = 1.584 
mu = 398600*10**9
Cx = 0.85
Cy0 = 0.134*0.85

TotalMass = 5500
F = True

result = open('landing Earth.txt', 'w')
result.close()

def Safety(F1, F2, F3, F):
    np = (F1**2 + F2**2 + F3**2)**0.5/g_Earth #в относительных единицах
    if np>=10:
        print ("Pilot is dead ;C ", np)
        F = False
    overload.append(np)
    return F


def fout1(t, y):# обработчик шага 
        ts.append(t)
        ys.append(list(y.copy()))
        y1, y2, y3, y4 = y
        if m.sqrt(y2*y2+y4*y4)<=300 and (m.sqrt(y1*y1 +y3*y3) - R_Earth)<10000:
          return -1
      
def rho(x,y):
    '''с помощью линейной апроксимации определяет плотность воздуха на необходимой нам высоте'''
    Space = [0, 1.85 * 0.00001, 1.5*0.0001, 3 * 0.0001, 1.03 * 0.001, 4 * 0.001, 7.26 * 0.001, 0.0136, 0.0251, 0.0469, 0.0889, 0.1216, 0.1665, 0.2279, 0.3119, 0.3648, 0.4135, 0.4671, 0.5258, 0.59, 0.6601, 0.7365, 0.8194, 0.9093, 1,1.1]
            # плотность для разных высот
    Space_lst = [100000, 80000, 70000, 60000, 50000, 40000, 36000, 32000, 28000, 24000, 20000, 18000, 16000, 14000, 12000, 11000, 10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000,2000,1000]    
   
    h = m.sqrt(x*x +y*y) - R_Earth
    
    i = 25
    while h > Space_lst[i]:
        i-=1
        if i < 0:
            return 0
           
    delta = h - Space_lst[i]
# разница между высотой и ближайшим значением
    delta_h = Space_lst[i-1] - Space_lst[i]
# разница между ближайшими соседями
    otn = delta/delta_h
# относительное отклонение
    p = Space[i] + ((Space[i-1] - Space[i]) * otn)
    return p

def tangage_angle(x,y):
    Space_lst = [80000, 70000, 60000, 50000, 30000, 10000]
    angle_lst = [6.5, 6.5, 6.4, 6.3, 6.2, 6.1]
    Cy_lst = np.zeros(7)
    Cy_lst[0] = 0.134*0.85*m.cos(6.5*m.pi/180)
    for i in range(1,6):
        if Cy_lst[i-1]!=Cy_lst[i]:
            Cy_lst[i]=Cy_lst[i-1]*m.cos(angle_lst[i]*m.pi/180)
    h = m.sqrt(x*x +y*y) - R_Earth
    i = 0
    while h < Space_lst[i]:
        i+=1
        if i >=5:
           break
    return angle_lst[i]*m.pi/180, Cy_lst[i]
        
# функция правых частей системы ОДУ
def f1(t, y):
         global TotalMass, max_overload, F    
         if F:
             y1, y2, y3, y4 = y
             result = open('landing Earth.txt', 'a')         
             vv=y2*y2+y4*y4
             v = m.sqrt(vv)
             r = (y1**2+y3**2)**0.5
             angle = tangage_angle(y1,y3)[0]
             Cy = tangage_angle(y1,y3)[1]
             p = rho(y1,y3)
             Q = (Cx*S*p*(v**2))/(2*TotalMass)
             N = (Cy*S*p*(v**2))/(2*TotalMass)
             Py = N*m.cos(angle)
             Pz = N*m.sin(angle)
             F = Safety(Q, Py, Pz, F)
             ax = -Q*(y2/v) - Py*(y4/v) - (mu/r**2)*y1/r
             ay = -Q*(y4/v) + Py*(y2/v) - (mu/r**2)*y3/r 
             result.write(str(r) + '\t' + str(m.sqrt(vv)) + '\t' +
                          str(m.sqrt(ax**2 +ay**2)) + '\t' + str(t) + '\n')
             result.close()
             return [y2, ax, y4,ay] 
     
x_start = 0
Vx_start = vx0
y_start = h + R_Earth
Vy_start = -vy0

xc,yc=[],[]
for i in range(90, 159):
    xc.append(0.001*R_Earth*m.cos(i/100))
    yc.append(0.001*R_Earth*m.sin(i/100))
tmax = 100000
    
y0,t0=[x_start,  Vx_start, y_start, Vy_start], 0 # начальные условия 
ODE=ode(f1)
ODE.set_integrator('dopri5')
ODE.set_solout(fout1)
ts, ys = [ ],[ ]
ODE.set_initial_value(y0, t0) # задание начальных значений
ODE.integrate(tmax)      # решение ОДУ
Y=np.array(ys)
xgraph = Y[:,0]/1000
ygraph = Y[:,2]/1000


print()
print('Конечная высота ',"%.3f" % (m.sqrt(Y[-1:,0]**2 + Y[-1:,2]**2) - R_Earth),'м')
print('Конечная скорость ', "%.3f" % (m.sqrt(Y[-1:,1]*Y[-1:,1] + Y[-1:,3]*Y[-1:,3])),'м/с')
totalTime+=ts[-1]
print('Полное время снижения в атмосфере ', "%.0f" %totalTime,'с')
max_overload = max(overload)
print ('Максимальная перегрузка ', "%.1f" %max_overload, 'g')

end = open('FinalData.txt', 'w')
end.write('Конечная высота '+ str(round(m.sqrt(Y[-1:,0]**2 + Y[-1:,2]**2) - R_Earth)) +' м \n' +
          'Конечная скорость '+ str(round(m.sqrt(Y[-1:,1]*Y[-1:,1] + Y[-1:,3]*Y[-1:,3]))) +' м/с \n'+
          'Полное время снижения в атмосфере '+ str(round(totalTime)) +' с \n' +
          'Максимальная перегрузка ' + str(round(max_overload,1)) + ' g \n' +
          'Ура! Мы вернулись на Землю!')
end.close()
plt.style.use('dark_background')
plt.title("Торможение в атмосфере \n ")
x = np.zeros(39)
y = np.zeros(39)
for i in range(39):
    if i<27:
        x[i] = i*100
        y[i] = random.randint(6375, 7000)
    else:
        x[i] = i*100 
        y[i] = random.randint(5500, 7000)       
    
plt.plot(x, y, marker ="*", c="white", linestyle=" ")
plt.plot(xgraph, ygraph, linewidth=2, color='lightcoral')
plt.axis('equal')
plt.text (1000, 5000, u'the Earth', size=40, color="bisque")
plt.plot(xc, yc, linewidth=3, color='olivedrab')
plt.xlim(0, 3500)
plt.ylim(5000, 7000)
plt.xlabel('X, км')
plt.ylabel('Y, км')
plt.grid(False)
plt.show()

def visual(data):
    for i in range(0,len(data)):
        y = data[i]
        plt.plot(i, y, c='lightcoral',marker = ".", markersize=5, linestyle="--")
        plt.pause(0.001)
    

dataa = np.zeros(len(Y[:,2]))
for i in range(len(Y[:,2])):
    dataa[i] = (m.sqrt(Y[i][0]**2 +Y[i][2]**2) - R_Earth)/1000
    

plt.figure(2)
plt.subplot(131)
plt.title("Высота над поверхностью земли \n")
plt.xlabel('Время, с')
plt.ylabel('Высота, км')
plt.plot(ts[::],dataa,c='plum',linewidth=3)

speeds = np.zeros(len(Y[:,2]))
for i in range(len(Y[:,1])):
    speeds[i] = m.sqrt(Y[i][1]**2 +Y[i][3]**2)
 
plt.subplot(133)
plt.plot(ts[::],speeds,c='plum',linewidth=3)
plt.title("Скорость КО \n")
plt.xlabel('Время, с')
plt.ylabel('Скорость, м/с')
plt.show()  