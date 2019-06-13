from math import *
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import matplotlib.transforms as trns
import pylab
from numpy import *
import random
import os

G = 6.67e-11  # Гравитационная постоянная
M_Moon = 7.36e22  # Масса Луны
R_Moon = 1738000  # Радиус луны
m_ship = 6835 #(пс+вс)лм
m_lm = 2315
m_fuel_1 = 8165
m_fuel_2 = 2355
F_1 = 45040
F_2 = 15600
h_0_Moon = 2
u_fuel = 3050
dt = 0.01

X, Y, R, V, A, T = [], [], [], [], [], []

file = open('FromEarthToMoonOut.txt').readlines()
prev = array([[float(i) for i in file[j].split()] for j in range((len(file)))])
x_Moon, y_Moon, x_ship, y_ship, FullT = ((prev[0][0])/1000), ((prev[0][1])/1000), ((prev[0][2])/1000), ((prev[0][3])/1000), prev[0][6]

def change_coordinates (x_Moon, y_Moon, x_ship, y_ship):
    New_x =(x_ship - x_Moon)*(y_Moon/sqrt(x_Moon**2 + y_Moon**2)) - (y_ship - y_Moon)*(x_Moon/sqrt(x_Moon**2 + y_Moon**2))
    New_y =(x_ship - x_Moon)*(x_Moon/sqrt(x_Moon**2 + y_Moon**2)) + (y_ship - y_Moon)*(y_Moon/sqrt(x_Moon**2 + y_Moon**2))
    R_orbit = sqrt(New_x**2 + New_y**2)*1000
    V_orbit = sqrt(G * M_Moon / R_orbit)
    return New_x, New_y, R_orbit, V_orbit

new_x, new_y, R_orbit, V_orbit = change_coordinates (x_Moon, y_Moon, x_ship, y_ship)

def gg(x, y, G, M_Moon):  # гравитация
    gx = - (G*M_Moon*x)/((x**2 + y**2)**1.5)
    gy = - (G*M_Moon*y)/((x**2 + y**2)**1.5)
    return gx, gy

def current_position(x, y, r):  # местоположение РН на орбите
    if x <= 0:
        alpha = acos(y/r)
    else:
        alpha = acos(-y/r) + pi
    position = r*alpha
    return position

def newData(x, y, vx, vy, ax, ay):
    x = x + vx*dt + ax*(dt**2)/2
    y = y + vy*dt + ay*(dt**2)/2
    vx = vx + ax*dt
    vy = vy + ay*dt
    return x, y, vx, vy

def free(x, y, vx, vy):  # свободный полет
    while x <= 0:
        ax, ay = gg(x, y, G, M_Moon)
        x, y, vx, vy = newData(x, y, vx, vy, ax, ay)
    return x, y, vx, vy

def desired_x(x=0, y=R_orbit, vx=-V_orbit, vy=0):  # место начала торможения
    height = y
    speed = vx
    x_i = x
    y_i = y
    vy_i = vy
    while height > R_Moon + ((y_i - R_Moon) / 2):  # ((y_probe - rad) * 0.95):
        speed += 0.1
        x = x_i
        y = y_i
        vx = speed
        vy = vy_i
        x, y, vx, vy = free(x, y, vx, vy)
        height = abs(y)
    desiredx = (-speed**2 + V_orbit**2)/(2*F_1/(m_ship + m_fuel_1))
    return desiredx

def after_impulse(x_ship_orbit, y_ship_orbit, mfuel = m_fuel_1): # скорость и координаты перед началом свободного полёта
    global X, Y, R, V, A, T
    x = desired_x()
    y = sqrt(R_orbit**2 - x**2)
    decel_startco = current_position(x, y, R_orbit)
    rocket_co = current_position(x_ship_orbit, y_ship_orbit, R_orbit)
    if decel_startco >= rocket_co:
        t_1stage = (decel_startco - rocket_co) / V_orbit
    else:
        t_1stage = (decel_startco + 2 * pi * R_orbit - rocket_co) / V_orbit
    time = 0
    vx = - (V_orbit*y)/R_orbit
    vy = (V_orbit*x)/R_orbit
    ax = F_1/(m_ship + m_fuel_1)
    while x > 0:
        gx, ay = gg(x, y,  G, M_Moon)
        afuel = ax - gx
        x, y, vx, vy = newData(x, y, vx, vy, ax, ay)
        mfuel -= ((m_ship + mfuel)*afuel/u_fuel)*dt
        time += dt
        X.append(x/1000), Y.append(y/1000), R.append(sqrt(x ** 2 + y ** 2)/1000), V.append(vx ** 2 + vy ** 2)
        A.append(sqrt(ax ** 2 + ay ** 2)), T.append(time)
    return mfuel, vx, vy, x, y, t_1stage, time, decel_startco

mfuel_aftdecel, vx, vy, x, y, t_1stage, time, decel_startco = after_impulse(new_x, new_y)

print('Конец торможения, начало свободного полёта')
print("x =", round(x, 4), "y =", round(-y, 4), "Vx =", round(vx, 4), "Vy =", round(-vy, 4), "масса топлива =",
      round(mfuel_aftdecel, 4))

def writing_coordinates(x, y, vx, vy, mfuel, time, f_max=F_1):  # свободный полёт
    global X, Y, R, V, A, T
    n = 1
    while y > 0:  # Половина траектории
        ax, ay = gg(x, y, G, M_Moon)
        x, y, vx, vy = newData(x, y, vx, vy, ax, ay)
        time += dt
        n += 1
        if n % 10 == 0:
            X.append(x/1000), Y.append(y/1000), R.append(sqrt(x ** 2 + y ** 2)/1000), V.append(vx ** 2 + vy ** 2)
            A.append(sqrt(ax ** 2 + ay ** 2)), T.append(time)
    ax, ay = gg(x, y, G, M_Moon)
    while (vx ** 2) < abs(2 * x * (- (f_max/(m_ship + mfuel)) + ax)):  # скорость при которой остаётся горизонтальная составляющая V
        ax, ay = gg(x, y, G, M_Moon)
        x, y, vx, vy = newData(x, y, vx, vy, ax, ay)
        time += dt
        n += 1
        if n % 10 == 0:
            X.append(x/1000), Y.append(y/1000), R.append(sqrt(x ** 2 + y ** 2)/1000), V.append(vx ** 2 + vy ** 2)
            A.append(sqrt(ax ** 2 + ay ** 2)), T.append(time)
    return x, y, vx, vy, time


x_1, y_1, vx_1, vy_1, time = writing_coordinates(x, y, vx, vy, mfuel_aftdecel, time)
print("Конец свободного полёта, начало торможения")
print("x =", round(x_1, 4), "y =", round(-y_1, 4), "Vx =", round(vx_1, 4), "Vy =", round(-vy_1, 4))

def stop(x, y, vx, vy, mfuel, time, f_max=F_1):  # торможение
    global X, Y, R, V, A, T
    ax_G, ay_G = gg(x, y,  G, M_Moon)
    ax = - (F_1/ (mfuel + m_ship)) + ax_G
    while vx > 0:
        ax_G, ay_G = gg(x, y,  G, M_Moon)
        ax_decel = - ax + ax_G
        ay_decel = - sqrt(((f_max / (m_ship + mfuel)) ** 2) - (ax_decel ** 2))
        ay = ay_G + ay_decel
        if mfuel < 7800:
            if vy > 100:
                mfuel -= (f_max /u_fuel) * dt
            else:
                ay = ay_G
                mfuel -= (ax_decel * (m_ship + mfuel) / u_fuel) * dt
        else:
            mfuel -= (f_max / u_fuel) * dt
        time += dt
        x, y, vx, vy = newData(x, y, vx, vy, ax, ay)
        X.append(x/1000), Y.append(y/1000), R.append(sqrt(x ** 2 + y ** 2)/1000), V.append(vx ** 2 + vy ** 2)
        A.append(sqrt(ax ** 2 + ay ** 2)), T.append(time)
    return x, y, vx, vy, mfuel, time

print("Конец торможения, начало вертикальной посадки")
x_2, y_2, vx_2, vy_2, mfuel_2, time = stop(x_1, y_1, vx_1, vy_1, mfuel_aftdecel,time)
print("x =", round(x_2, 4), "y =", round(-y_2, 4), "Vx =", round(vx_2, 4), "Vy =", round(-vy_2, 4),
      "масса топлива =", round(mfuel_2, 4))


def vertical_stop(x, y, vx, vy, mfuel, time):  # вертикальная посадка
    global X, Y, R, V, A, T
    while vy > 0:
        ax_G, ay_G = gg(x, y,  G, M_Moon)
        ay = - (vy ** 2 / (2 * abs(abs(y) - R_Moon)))
        x, y, vx, vy = newData(x, y, vx, vy, ax_G, ay)
        mfuel -= ((ay_G - ay) * (m_ship + mfuel) / u_fuel) * dt
        time += dt
        X.append(x/1000), Y.append(y/1000), R.append(sqrt(x ** 2 + y ** 2)/1000), V.append(vx ** 2 + vy ** 2)
        A.append(sqrt(ax_G ** 2 + ay ** 2)), T.append(time)
    return x, y, vx, vy, mfuel, time

print('Успешная посадка ⚐')
x_3, y_3, vx_3, vy_3, mfuel_3, t_2stage = vertical_stop(x_2, y_2, vx_2, vy_2, mfuel_2, time)
print("Место прилунения: x =", round(x_3, 4), "y =", round(-y_3, 4), "Скорость при посадке: v_x =", round(
    vx_3, 4),  "v_y =", round(-vy_3, 4), "Масса топлива =", round(mfuel_3, 4))

t_landing = t_1stage + t_2stage
position_of_lk = (decel_startco + V_orbit * t_2stage)% (2 * pi * R_orbit)
y_lk = R_orbit * cos(position_of_lk / R_orbit)
x_lk = - R_orbit * sin(position_of_lk / R_orbit)  # местоположение ракеты после посадки

print("Местоположения ракеты носителя на орбите: x =", round(x_lk, 4), "y =", round(-y_lk, 4), "Полное время =", round(
    t_landing + FullT, 4))


plt.style.use('dark_background')
fig, ax = plt.subplots()
xc,yc=[],[]
for i in range(0, 630):
    xc.append(R_Moon*0.001*cos(i/100))
    yc.append(R_Moon*0.001*sin(i/100))
plt.plot(xc,yc,linewidth=2, c = 'bisque')

x = zeros(39)
y = zeros(39)
for i in range(39):
    x[i] = random.randint(-2000, 2000)
    y[i] = random.randint(-2000, 2000)
plt.plot(x, y, marker="*", c="lightsteelblue", linestyle=" ")

plt.axis('equal')
plt.plot(list(X), list(Y), marker = "*", c="lightcoral", markersize=0.1)
plt.xlabel("Координата x, км")
plt.ylabel("Координата y, км")
ax.set_title("Траектория ракеты вблизи луны")
pylab.xlim(-2000, 2000)
pylab.ylim(-2000, 2000)
plt.grid(False)
plt.savefig("Moon1.png", dpi=300)
plt.show()

fig, ax = plt.subplots()

x = zeros(39)
y = zeros(39)
for i in range(39):
    x[i] = random.randint(0, 4000)
    y[i] = random.randint(0, 3000000)
plt.plot(x, y, marker="*", c="powderblue", linestyle=" ")

plt.plot(list(T), list(V), marker = "*", c = "olive", markersize=0.1)
plt.ylabel('Скорость, м/с ')
plt.xlabel("Время, с")
ax.set_title("Скорость ЛМ")
plt.grid(False)
plt.savefig("Moon2.png", dpi=300)
plt.show()

fig, ax = plt.subplots()

x = zeros(39)
y = zeros(39)
for i in range(39):
    x[i] = random.randint(0, 4000)
    y[i] = random.randint(1738, 1800)
plt.plot(x, y, marker="*", c="lightsteelblue", linestyle=" ")

plt.plot(list(T), list(R), marker = "*", c = "mediumvioletred", markersize=0.1)
plt.ylabel("Расстояние до поверхности луны, км")
plt.xlabel("Время, c")
ax.set_title("Расстояние от центра луны")
plt.grid(False)
plt.savefig("Moon3.png", dpi=300)
plt.show()

x, y, Vx, Vy, t_3stage = x_3, -y_3, 0, 0, 0
X, Y, R, V, T = [], [], [], [], []

def ready_steady(angle):
    dt = 1
    global x, y, Vx, Vy, m_fuel_2, m_lm, t_3stage, X, Y, R, V, T
    gx = (G * M_Moon * x) / ((x ** 2 + y ** 2) ** 1.5)
    gy = (G * M_Moon * y) / ((x ** 2 + y ** 2) ** 1.5)
    m0 = m_lm + m_fuel_2
    m = m0 - (F_2/u_fuel) * dt
    Vr = u_fuel*log(m/m0)
    x = x + Vx * dt + 0.5 * gx * dt ** 2 - cos(radians(angle)) * ((-m0 / (F_2/u_fuel) + dt) * Vr - u_fuel*dt)
    y = y + Vy * dt - 0.5 * gy * dt ** 2 - sin(radians(angle)) * ((-m0 / (F_2/u_fuel) + dt) * Vr - u_fuel*dt)
    Vx = Vx + gx * dt - cos(radians(angle)) * Vr
    Vy = Vy - gy * dt - sin(radians(angle)) * Vr
    m_fuel_2 = m_fuel_2 - (F_2/u_fuel) * dt
    t_3stage += dt
    X.append(x / 1000), Y.append(y / 1000), R.append(sqrt(x ** 2 + y ** 2) / 1000), V.append(sqrt(Vx ** 2 + Vy ** 2))
    T.append(t_3stage)



def go(x_0, y_0):
    global x, y, t_3stage, X, Y, R, V, T, m_fuel_2
    t_3stage, i = 0, 0
    x = x_0
    y = y_0
    while sqrt(x ** 2 + y ** 2)  <= 1740000:
        i += 1
        ready_steady(90 - i)
    while sqrt(x ** 2 + y ** 2) <= 1748000:
        ready_steady(38)
    i = 0
    while sqrt(x ** 2 + y ** 2) <= 1764100:
        i += 0.25
        ready_steady(38-i)
    while sqrt(x ** 2 + y ** 2) <= 1788500:
        ready_steady((acos(y/sqrt(x ** 2 + y ** 2)))*(180/pi))
    i=0
    while i<35:
        i+=1
        ready_steady(270)
    i=0
    while i < 1:
        i += 1
        ready_steady(0)
    V4 = sqrt(Vx ** 2 + Vy ** 2)
    H4 = sqrt(x ** 2 + y ** 2)
    print("Место стыковки: x =", round(x, 4), "y =", round(y, 4), "Скорость на орбите: v =", round(
    V4, 4), "Высота орбиты =", round(H4, 4), "Масса топлива =", round(m_fuel_2, 4))
    Esc_moon = open("Escaping from moon.txt", 'w')
    Esc_moon.write(str(atan(x/y)) + "\n" + str(sqrt(Vx ** 2 + Vy ** 2)) + "\n" + str(sqrt(x ** 2 + y ** 2))
                   + "\n" + str(t_landing + 3202 + FullT))
    Esc_moon.close()
    drawing(X, Y, R, V, T)
    return (V4, H4, m_fuel_2)


def drawing(X, Y, R, V, T ):
    plt.style.use('dark_background')
    fig, ax = plt.subplots()
    xc, yc = [], []
    for i in range(0, 630):
        xc.append(R_Moon * 0.001 * cos(i / 100))
        yc.append(R_Moon * 0.001 * sin(i / 100))
    plt.plot(xc, yc, linewidth=2, c='bisque')

    x = zeros(39)
    y = zeros(39)
    for i in range(39):
        x[i] = random.randint(-75, 365)
        y[i] = random.randint(1650, 1790)
    plt.plot(x, y, marker="*", c="lightsteelblue", linestyle=" ")

    plt.plot(list(X), list(Y), marker="*", c="mediumaquamarine", markersize=0.1)
    plt.xlabel("Координата x, км")
    plt.ylabel("Координата y, км")
    ax.set_title("Траектория ракеты вблизи луны")
    pylab.xlim(-75, 365)
    pylab.ylim(1650, 1790)
    plt.grid(False)
    plt.rcParams['examples.directory'] = os.getcwd()
    image_data = cbook.get_sample_data('rocket.png')
    image = plt.imread(image_data)
    im = ax.imshow(image, origin='lower', extent=[327, 335, 1749, 1772])
    trans_data = trns.Affine2D().rotate_deg_around(333, 1761, 75) + ax.transData
    im.set_transform(trans_data)
    plt.savefig("Moon4.png", dpi=300)
    plt.show()

    plt.style.use('dark_background')
    fig, ax = plt.subplots()

    x = zeros(39)
    y = zeros(39)
    for i in range(39):
        x[i] = random.randint(0, 500)
        y[i] = random.randint(0, 2000)
    plt.plot(x, y, marker="*", c="lightsteelblue", linestyle=" ")

    plt.plot(list(T), list(V), marker="*", c="lightcoral", markersize=0.1)
    plt.xlabel("Время, с")
    plt.ylabel("Скорость, м/с")
    ax.set_title("Скорость ЛМ")
    pylab.xlim(0, 500)
    pylab.ylim(0, 2000)
    plt.grid(False)
    plt.savefig("Moon5.png", dpi=300)
    plt.show()

    plt.style.use('dark_background')
    fig, ax = plt.subplots()

    x = zeros(39)
    y = zeros(39)
    for i in range(39):
        x[i] = random.randint(0, 500)
        y[i] = random.randint(1738, 1800)
    plt.plot(x, y, marker="*", c="lightsteelblue", linestyle=" ")

    plt.plot(list(T), list(R), marker="*", c="darkturquoise", markersize=0.1)
    plt.xlabel("Время, с")
    plt.ylabel("Высота, км")
    ax.set_title("Расстояние от центра луны")
    pylab.xlim(0, 500)
    pylab.ylim(1738, 1800)
    plt.grid(False)
    plt.savefig("Moon6.png", dpi=300)
    plt.show()

t_landing = t_1stage + t_2stage
position_of_lk = (decel_startco + V_orbit * t_2stage) % (2 * pi * R_orbit)
y_lk = R_orbit * cos(position_of_lk / R_orbit)
x_lk = - R_orbit * sin(position_of_lk / R_orbit)



def vzlet():
    position_of_lk1 = (decel_startco + V_orbit * (
    ((R_orbit * (asin(-319716.9245 / R_orbit) + 3*pi) - decel_startco) / V_orbit))) % (2 * pi * R_orbit)
    t_blast = (R_orbit * (asin(-319716.9245 / R_orbit) + 3*pi) - decel_startco) / V_orbit
    y_lk1 = R_orbit * cos(position_of_lk1 / R_orbit)
    x_lk1 = R_orbit * sin(position_of_lk1 / R_orbit)
    go(0, R_Moon)
    print("Чилим на луне", round((t_blast - t_3stage)/60),"минут ♫ ℂℍℐℒℒ ♫")
    print("Местоположения ракеты носителя на орбите. Момент стыковки: x =", round(x_lk1, 4), "y =", round(-y_lk1, 4),
          "Полное время =", round(t_landing+t_blast+FullT))


vzlet()
print ("Ура! Летим на землю!🌍")
