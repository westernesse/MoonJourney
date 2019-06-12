import math
from scipy.integrate import ode
import matplotlib.pyplot as plt

M_lk = 28000
M_lm = 15000
M1 = 135000
M1_fuel = 2145000 - M1   # kg
force1 = 34350 * 10 ** 3
f_vel1 = 2580
d1 = 10.1
M2 = 37600
M2_fuel = 458700 - M2
force2 = 5115 * 10 ** 3
f_vel2 = 4130
d2 = 10.1
M3 = 13000
M3_fuel = 113000 - M3
force3 = 1016 * 10 ** 3
f_vel3 = 4130
d3 = 6.6
M_total = M1 + M2 + M1_fuel + M2_fuel + M3 + M3_fuel + M_lk + M_lm
Me = 5.972 * 10 ** 24  # kg  Mass of the Earth
We = 7.29 * 10 ** (-5)
G = 6.67408 * 10 ** (-11)
R = 6375000  # m
h = 350000
vt_ = []
vr_ = []
xs = []
ys = []
xsf = []
ysf = []
h_ = []
t_ = []


def fuel1(t, y):
    global xs, ys, force1, M_total, M1_fuel, f_vel1
    xs.append(y[0])
    ys.append(y[2])
    h_.append(math.sqrt(y[0]**2 + y[2]**2) - R)
    vr_.append(- (y[1]*math.cos(math.pi - math.atan(y[2]/y[0])) + y[3]*math.cos(math.pi - math.atan(y[0]/y[2]))))
    vt_.append(math.sqrt(y[1]**2 + y[3]**2 - vr_[-1] **2 ))
    t_.append(t)
    if M1_fuel - force1/f_vel1 * t <= 0.001:
        M_total -= M1 + M1_fuel
        print("First stage detached")


def fuel2(t, y):
    global xs, ys, force2, M_total, f_vel2
    xs.append(y[0])
    ys.append(y[2])
    h_.append(math.sqrt(y[0]**2 + y[2]**2) - R)
    vr_.append(- (y[1]*math.cos(math.pi - math.atan(y[2]/y[0])) + y[3]*math.cos(math.pi - math.atan(y[0]/y[2]))))
    vt_.append(math.sqrt(y[1] ** 2 + y[3] ** 2 - vr_[-1] ** 2))
    t_.append(t)
    if M2_fuel - force2/f_vel2 * t <= 0.001:
        M_total -= M2 + M2_fuel
        print("Second stage detached")


def fuel3(t, y):
    global xs, ys, M_total, f_vel3
    xs.append(y[0])
    ys.append(y[2])
    h_.append(math.sqrt(y[0]**2 + y[2]**2) - R)
    vr_.append(- (y[1]*math.cos(math.pi - math.atan(y[2]/y[0])) + y[3]*math.cos(math.pi - math.atan(y[0]/y[2]))))
    vt_.append(math.sqrt(y[1] ** 2 + y[3] ** 2 - vr_[-1] ** 2))
    t_.append(t)
    h = math.sqrt(y[0]**2 + y[2]**2)


def fuel4(t, y):
    global xs, ys, force2, M_total, f_vel2
    xsf.append(y[0])
    ysf.append(y[2])


def angle(x, y):
    angle = 90
    h = math.sqrt(x*x + y*y) - R
    H = (10000, 12000, 14000, 16000, 18000, 20000, 22000, 24000, 70000, 75000, 80000, 90000, 100000, 110000, 185000, 198000)
    ang = (0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85)
    for i in range(16):
        if h <= H[i]:
            angle = ang[i]
            break
    return math.radians(angle)


def atm_dens(x, y):
    h1 = math.sqrt(x*x + y*y)
    if 0 < h1 <= 11000:
        D = 1.225 * (1 - h1 / 44331)**4.255876
    elif 11000 < h1 <= 20000:
        D = 1.225 * 0.297076 * 10**(-(h1-11000)/6375000)
    elif 20000 < h1 <= 32000:
        D = 1.225 * (0.978261 + h1 / 201020)**(-35.16320)
    elif 32000 < h1 <= 47000:
        D = 1.225 * (0.857003 + h1 / 57947)**(-13.20114)
    elif 47000 < h1 <= 51000:
        D = 1.225 * 0.00116533 * 10**(-(h1-47000)/7922)
    elif 51000 < h1 <= 71000:
        D = 1.225 * (0.798990 - h1 / 184809)**11.20114
    elif h1 > 71000:
        D = 0
    return (D)


def launch1(t, y):  # Запуск первой ступени
    global f_vel1, M1, M1_fuel, M_total, force1, d1
    X, X1, Y, Y1 = y
    ang = angle(X, Y) + math.atan(-Y/X)
    S = math.pi * d1 ** 2 / 4  # Параметры для силы сопротивления
    C = 0.1
    v2 = X1**2 + Y1**2
    D = atm_dens(X, Y)
    dM_fuel = force1/f_vel1 * t
    a = force1/(M_total - dM_fuel)
    res = C*S*v2*D / 2     # сила сопротивления

    ax = - X * G * Me / ((X*X + Y*Y) ** 1.5) + (a - res) * math.cos(ang)
    ay = - Y * G * Me / ((X*X + Y*Y) ** 1.5) - (a - res) * math.sin(ang)

    return [X1, ax, Y1, ay]


def launch2(t, y):  # Запуск второй ступени
    global f_vel2, M2, M2_fuel, M_total, force2, d2
    X, X1, Y, Y1 = y
    ang = angle(X, Y) + math.atan(-Y/X)
    S = math.pi * d2 ** 2 / 4
    C = 0.1
    v2 = X1**2 + Y1**2
    D = atm_dens(X, Y)

    dM_fuel = force2/f_vel2 * t
    a = force2/(M_total - dM_fuel)

    res = C*S*v2*D / 2
    ax = - X * G * Me / ((X*X + Y*Y) ** 1.5) + (a - res) * math.cos(ang)
    ay = - Y * G * Me / ((X*X + Y*Y) ** 1.5) - (a - res) * math.sin(ang)

    return [X1, ax, Y1, ay]


def launch3(t, y):  # Запуск третьей ступени
    global f_vel3, M3, M3_fuel, M_total, force3, d3
    X, X1, Y, Y1 = y
    ang = angle(X, Y) + math.atan(-Y/X)
    S = math.pi * d3 ** 2 / 4  # Параметры для силы сопротивления
    C = 0.1
    v2 = X1**2 + Y1**2
    D = atm_dens(X, Y)

    dM_fuel = force3/f_vel3 * t
    a = force3/(M_total - dM_fuel)
    res = C*S*v2*D / 2     # сила сопротивления

    ax = - X * G * Me / ((X*X + Y*Y) ** 1.5) + (a - res) * math.cos(ang)
    ay = - Y * G * Me / ((X*X + Y*Y) ** 1.5) - (a - res) * math.sin(ang)

    return [X1, ax, Y1, ay]


def free_flight(t, y):
    X, X1, Y, Y1 = y
    ax = - X * G * Me / ((X * X + Y * Y) ** 1.5)
    ay = - Y * G * Me / ((X * X + Y * Y) ** 1.5)
    return [X1, ax, Y1, ay]


y0 = [R, 0, 0, -We * R]  # Начальные условия
t0 = 0

f = ode(launch1)
f.set_integrator('dopri5', nsteps=100000)
f.set_initial_value(y0, t0)
f.set_solout(fuel1) # функция, которая вызывается после каждой итерации

y1 = f.integrate(M1_fuel*f_vel1 / force1)  # решение системы для первой ступени; аргумент - время, возвращает список: x, x1, y, y1
t1 = M1_fuel*f_vel1 / force1


f = ode(launch2)
f.set_integrator('dopri5', nsteps=100000)
f.set_initial_value(y1, t1)
f.set_solout(fuel2) # функция, которая вызывается после каждой итерации

y2 = f.integrate(M2_fuel*f_vel2 / force2)  # решение для второй ступени
t2 = M2_fuel*f_vel2 / force2


f = ode(launch3)
f.set_integrator('dopri5', nsteps=100000)
f.set_initial_value(y2, t2)
f.set_solout(fuel3) # функция, которая вызывается после каждой итерации

# y3 = f.integrate(t2 + (M3_fuel * 0.35) * f_vel3 / force3) 572.8
y3 = f.integrate(487.5)
t3 = 487.5

print("free flight")
print("Height: ",math.sqrt(y3[0] ** 2 + y3[2] ** 2) - R, ", Velocity: ", math.sqrt(y3[1] ** 2 + y3[3] ** 2), ", Fuel left: ", (M3_fuel - force3/f_vel3 * (t3 - t2)) / M3_fuel, sep="")

f = ode(free_flight)
f.set_integrator('dopri5')
f.set_initial_value(y3, t3)
f.set_solout(fuel4)
f.integrate(10000)


x, vx, y, vy = y3
dm = (M3_fuel - force3/f_vel3 * (t3 - t2)) / M3_fuel
h = math.sqrt(x**2 + y**2)
nx, ny = math.cos(math.atan(-y/x)), math.sin(math.atan(-y/x))
w = math.sqrt((ny * vx + nx * vy)**2) / h
phi = math.atan(-y/x)
f = open("BlastOffFromEarth.txt", 'w')
print(dm, w, h, t3, phi, file=f)
f.close()

xc, yc = [], [] # рисует окружность Земли для графика
for i in range(0, 630):
    xc.append(R*math.cos(i/100))
    yc.append(R*math.sin(i/100))

with plt.style.context("dark_background"):
    graph = plt.figure()
    graph.set_facecolor("black")
    graph = plt.gca()
    graph.set_facecolor("black")
    plt.plot(xc, yc, color='olive', linewidth="2")
    plt.axis('equal')
    plt.plot(xs, ys, linewidth='4', color="#FFD6F9")
    plt.plot(xsf, ysf, color='#E8DAD5', alpha=0.8, linewidth='0.5')
    plt.title("Выход на Низкую Околоземную орбиту", fontsize="15")
    plt.show()
    graph1 = plt.figure()
    graph1.add_subplot(311)
    plt.plot(t_, h_)
    plt.title("Высота РН", fontsize="15")
    plt.ylabel("Высота, м")
    graph1.add_subplot(312)
    plt.plot(t_, vr_)
    plt.title("Радиальная скорость", fontsize="15")
    plt.ylabel("Скорость, м/с")
    graph1.add_subplot(313)
    plt.plot(t_, vt_)
    plt.title("Тангенциальная скорость", fontsize="15")
    plt.xlabel("Время, с")
    plt.ylabel("Скорость, м/с")
    plt.show()
