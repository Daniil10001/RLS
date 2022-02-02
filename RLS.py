import os
global tx, bx
from collections import namedtuple
from math import *
from tkinter import *
from tkinter import messagebox as mb
from tkinter.filedialog import askopenfilename
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

root = Tk()
root.title("Рельсотрон Python")
s = ['Сила тока A', 'Ширина проводника м', 'Масса кг', 'Магнитное поле Тс',
     "Количество проводников", "Арэдиномический коэффицент",
     "Площадь миделя м*м",
     "Угол", "Длина ствола м", 'Время с']
root.rowconfigure(list(range(len(s))), minsize=0, weight=1)
root.columnconfigure([0, 1, 2], minsize=0, weight=1)


def opeen():
    s = askopenfilename()
    print(s, os.path.normpath(s))
    try:
        txt = open(os.path.normpath(s), 'r').read()
    except:
        mb.showerror('Error','Файл нельзя открыть')
        return 0
    a = []
    a = txt.split(';')
    print(a, len(bx))
    if len(a) == len(bx):
        for i in range(len(bx)):
            bx[i].delete(0, 'end')
            bx[i].insert(0, a[i])
    else:
        mb.showerror('Alert','Неправильный файл')


def writee():
    from tkinter import filedialog as fd
    a = []
    for i in range(len(bx)):
        a.append(bx[i].get())
    try:
        f = fd.asksaveasfile(mode='w', defaultextension=".txt", filetypes=(("TXT files", "*.txt"),)).write(';'.join(a))
    except:
        pass


mainmenu = Menu(root)
root.config(menu=mainmenu)

filemenu = Menu(mainmenu, tearoff=0)
filemenu.add_command(label="Открыть", command=opeen)
filemenu.add_command(label="Сохранить", command=writee)
mainmenu.add_cascade(label="Файл",
                     menu=filemenu)


# @numba.jit(nopython=False)
def fnc(ghyt):
    i, l, m, bf, p, cf, sp, af, l2, t = ghyt  # получение данных
    f = i * l * bf * p  # сила ампера
    v0 = (l2 * 2 * (f / m)) ** 0.5  # начальная скорость
    print(v0, f, f / p, p)
    print(1)
    params = {'legend.fontsize': 14, 'figure.figsize': (10, 7), 'axes.labelsize': 14,
              'axes.titlesize': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14}
    pylab.rcParams.update(params)
    print(2)
    a = np.array([-6.3759, -7.3012, -1.1817])
    b = np.array([-0.4754, -0.0096, -0.0068, -0.0120, 0.0042]);
    c = np.array([0.1803, 0.0872, -0.0153, 0.0145, 0]);
    mu = 398600.4415e9
    Re = 6371000.0
    params = namedtuple("params", "CD CL mass")
    params.mass = m  # 4000.0
    params.CD = cf  # 0.02  # Коэффицент лобового сопротивления
    params.CL = 0  # Коэффицент подъёмной силы
    params.Sm = sp  # 1  # Площадь миделя
    h0 = 0  # Начальная высота [м]
    # v0 = 1000  # Начальная скорость [м/c]
    theta0 = radians(af)  # 1.22173  # начальный угол наклона траектории [радиан]

    def rho(h):
        x = h / 50.0 - 1
        sa = a[0] + a[1] * x + a[2] * x * x
        sbc = sum((b[i] * np.cos((i + 1) * np.pi * x) + c[i] * np.sin((i + 1) * np.pi * x) for i in range(5)))
        return np.exp(sa + sbc)

    def g_acc(h):
        return mu / (Re + h) ** 2

    def dydt(t, y, p):
        r = y[0]
        v = y[1]
        theta = y[2]
        h = (r - Re) * 0.001
        q = rho(h) * v * v / 2
        g = g_acc(r + Re)
        dv = - q * p.CD * p.Sm / p.mass - g * np.sin(theta)
        dtheta = (q * p.CL * p.Sm / p.mass - (g - v * v / r) * np.cos(theta)) / v
        dr = v * np.sin(theta)
        dx = v * np.cos(theta)
        return (dr, dv, dtheta, dx)

    def event_h_eq_0(t, y):
        return y[0] - Re

    event_h_eq_0.direction = -1
    event_h_eq_0.terminal = True
    print(100)
    sol = integrate.solve_ivp(lambda t, y: dydt(t, y, params), [0, t], [Re + h0, v0, theta0, 0], method='BDF',
                              events=event_h_eq_0, rtol=1e-8)
    print(1000)
    fig3, ax3 = plt.subplots(figsize=(6, 6))

    af = sol.y[3]
    bf = sol.y[0]
    originX = originY = 0
    theta1 = np.arange(0, 360, 0.001)
    ax3.plot(np.cos(theta1 * np.pi / 180) * (6371000.0) + originX, np.sin(theta1 * np.pi / 180) * (6371000.0) + originY,
             'o', ms=1)
    theta = (af % 6371000.0) * 360 / 6371000.0
    tenDegreePtsX = np.array([np.cos(theta[i] * np.pi / 180) * (bf[i]) + originX for i in range(len(bf))])
    tenDegreePtsY = np.array([np.sin(theta[i] * np.pi / 180) * (bf[i]) + originY for i in range(len(bf))])
    ax3.plot(tenDegreePtsX, tenDegreePtsY, ms=1)
    ax3.set_title('Траектория полета')
    mx = max(max(tenDegreePtsX), max(tenDegreePtsY))
    ax3.set_xlim(-mx * 1.2, mx * 1.2)
    ax3.set_ylim(-mx * 1.2, mx * 1.2)

    fig2, ax2 = plt.subplots(figsize=(8, 6))
    fig, ax = plt.subplots(figsize=(9, 6))
    fig1, ax1 = plt.subplots(figsize=(9, 6))
    ax1.set_title('Зависимость высоты м от пути м')
    ax.set_title('Зависимость высоты м от времени с')
    ax2.set_title('Зависимость скорости м/с от времени с')
    ax.plot(sol.t, (sol.y[0] - Re))
    ax1.plot(sol.y[3], (sol.y[0] - Re))  #
    ax.minorticks_on()
    ax1.minorticks_on()
    ax2.minorticks_on()
    ax.grid(which='major', color='k', linewidth=0.5)
    ax1.grid(which='major', color='k', linewidth=0.5)
    ax.grid(which='minor', color='k', linestyle=':')
    ax1.grid(which='minor', color='k', linestyle=':')
    ax2.grid(which='major', color='k', linewidth=0.5)
    ax2.grid(which='minor', color='k', linestyle=':')
    ax.set_xlabel("t, c", fontsize=14)
    ax.set_ylabel("h, м", fontsize=14)
    ax2.set_xlabel("t, c", fontsize=14)
    ax2.set_ylabel("v, м/с", fontsize=14)
    ax1.set_xlabel("l, м", fontsize=14)
    ax1.set_ylabel("h, м", fontsize=14)
    ax2.plot(sol.t, sol.y[1])
    plt.show(block=False)


def raschot():
    global tx, bx
    st = []
    for i in range(len(bx)):
        if bx[i].get() == '':
            st.append(s[i])
    if len(st) > 0:
        print(st)
        mb.showerror(title="Alert!", message="Не все данные введены!\n" + ', '.join(st))
        return
    i = float(bx[0].get())
    l = float(bx[1].get())
    m = float(bx[2].get())
    b = float(bx[3].get())
    p = float(bx[4].get())
    c = float(bx[5].get())
    sp = float(bx[6].get())
    a = float(bx[7].get())
    l2 = float(bx[8].get())
    t = float(bx[9].get())
    fnc((i, l, m, b, p, c, sp, a, l2, t))
    print(t)


tx = [Label(master=root, text=s[i], foreground="#000", font="16") for i in range(len(s))]
bx = [Entry(master=root, foreground="#000", font="16") for i in range(len(s))]
for i in range(len(s)):
    tx[i].grid(column=1, row=i)
    bx[i].grid(column=2, row=i)
# tx = Label(text="Clicks 0", background="#555", foreground="#000",font="16")
# tx.place(relx=.0, rely=.0, anchor="nw", height=30, width=150, bordermode=OUTSIDE)
# pr = Entry(text="Clicks 0", background="#555", foreground="#000",font="16")
# pr.place(relx=.5, rely=.0, anchor="nw", height=30, width=150, bordermode=OUTSIDE)
# btn1 = Button(master=root,text="Траектория\nдо земли", background="#555", foreground="#ccc",
#              padx="20", pady="8", font="16")
# btn1.grid(column=0, row=0)
btn2 = Button(master=root, text="Траектория\nпо времени", background="#555", foreground="#ccc",
              padx="20", pady="8", font="16")
btn2.grid(column=0, row=0)

btn = Button(master=root, text="Рассчитать", background="#555", foreground="#ccc",
             padx="20", pady="8", font="16", command=raschot)
btn.grid(column=1, row=len(s))

root.mainloop()
