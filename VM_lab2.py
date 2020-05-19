
#import tkinter as tk
#root = tk.Tk()

#text_widget = tk.Text(root)
#scrollbar = tk.Scrollbar(root, orient="vertical", command=text_widget.yview)
#text_widget.configure(yscrollcommand=scrollbar.set)

#scrollbar.pack(side="right", fill="y")
#text_widget.pack(side="left", fill="both", expand=True)

#root.mainloop()



import matplotlib
matplotlib.use('TkAgg')
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from tkinter import *
import tkinter.ttk as ttk
import tkinter as tk
import time, threading
import math
from matplotlib import pyplot as plt




def run_through_method(A, B, C, F):
    size = len(F)
    i = 1;
    a = [1] * (size - 1);
    b = [1] * (size - 1);
    x = [1] * size
    a[0] = -C[0] / B[0]
    b[0] = F[0] / B[0]
    while (i < size  - 1):
        a[i] = -C[i] / (A[i] * a[i - 1] + B[i])
        b[i] = (F[i] - A[i] * b[i - 1]) / (A[i] * a[i - 1] + B[i])
        i = i + 1

    x[size - 1] = (F[size - 1] - A[size - 1] * b[size - 2]) / (B[size - 1] + A[size - 1] * a[size - 2])
    i = size - 2
    while (i > -1):
        x[i] = a[i] * x[i + 1] + b[i]
        i = i - 1
    return x



def decision(mpb, window, l, T, h_t, h_x, b0, b1, b2, f0, f1):
    if (h_t / (h_x * h_x) >= 0.25):
        tk.messagebox.showinfo(title="Error", message="Шаг по x и шаг по t не удовлетворяют условию устойчивости:\nτ/(h*h) < 1/4")
        return
    time1 = time.time()
    max = int(2 * T / h_t)
    mpb["maximum"] = max # для строки прогресса
    mpb["value"] = 0
    window.update_idletasks()
    n = int(l / h_x) + 1
    m = int(T / h_t)
    # коэффициенты для метода прогонки
    A = [0] * n
    B = [0] * n
    C = [0] * n
    f_ = [0] * n
    # начальные функции
    b = [0] * n
    f = [0] * n
    x = [0] * n
    # часть 1
    j = 0
    while (j < m):
        i = 0
        while (i < n):
            if (j == 0):
                x[i] = i * h_x
                f[i] = 1 / l + f0 * math.cos(math.pi * x[i] / l) + f1 * math.cos(2 * math.pi * x[i] / l)
                f_[i] = f[i]
                b[i] = b0 + b1 * math.cos(math.pi * x[i] / l) + b2 * math.cos(2 * math.pi * x[i] / l)
            else:
                f_[i] = w[i]
            B[i] = 1 + 2 * h_t / (h_x * h_x) - h_t * b[i] 
            C[i] = -h_t / (h_x * h_x)
            A[i] = -h_t / (h_x * h_x)
            i = i + 1
        # подсчет коэффициентов для 0 и n-1
        A[0] = 0
        B[0] = 1 + h_t / (h_x * h_x) - h_t * b[0]
        C[0] = -h_t / (h_x * h_x)
        A[n - 1] = -h_t / (h_x * h_x)
        B[n - 1] = 1 + h_t / (h_x * h_x) - h_t * b[n - 1]
        C[n - 1] = 0

        w = run_through_method(A, B, C, f_)
        mpb["value"] = mpb["value"] + 1
        window.update_idletasks()
        j = j + 1
    # подсчет интеграла
    w_I = w[0] + w[n-1]
    j = 1
    while(j < n - 1):
        if (j%2 == 1):
            w_I = w_I + 4 * w[j]
        else:
            w_I = w_I + 2 * w[j]
        j = j + 1
    w_I = w_I * h_x / 3
    for i in range(n):
        w[i] = w[i] / w_I
    j = 0
    y = [0] * n
    # часть 2
    I = 0
    while (j < m):
        i = 0
        while (i < n):
            if (j == 0):
                f[i] = 1 / l + f0 * math.cos(math.pi * x[i] / l) + f1 * math.cos(2 * math.pi * x[i] / l)
                f_[i] = f[i]
                b[i] = b0 + b1 * math.cos(math.pi * x[i] / l) + b2 * math.cos(2 * math.pi * x[i] / l)
            else:
                    f_[i] = y[i]
            B[i] = 1 + 2 * h_t / (h_x * h_x) - h_t * b[i] + h_t * I
            C[i] = -h_t / (h_x * h_x)
            A[i] = -h_t / (h_x * h_x)
            i = i + 1
        A[0] = 0
        B[0] = 1 + h_t / (h_x * h_x) - h_t * b[0] + h_t* I
        C[0] = -h_t / (h_x * h_x)
        A[n - 1] = -h_t / (h_x * h_x)
        B[n - 1] = 1 + h_t / (h_x * h_x) - h_t * b[n - 1] + h_t * I
        C[n - 1] = 0
        y = run_through_method(A, B, C, f_)
        # подсчет интеграла
        I = y[0] * b[0]+ y[n - 1] * b[n - 1]
        k = 1
        while (k < n - 1):
            if (k%2 == 1):
                I = I + 4 * y[k] * b[k]
            else:
                I = I + 2 * y[k] * b[k]
            k = k + 1
        I = I * h_x / 3
        mpb["value"] = mpb["value"] + 1
        window.update_idletasks()
        j = j + 1

    time2 = time.time()
    time_ = time2 - time1
    str_ = "Время выполнения: " + str(round(time_, 2))
    label = Label(text = str_, font = "Arial 10")
    label.grid(row=3, column=6)
    window.update_idletasks()

    global green, _x, a, canvas
    green = w
    _x = x
    fig = Figure(figsize=(7,5)) # размер графика
    a = fig.add_subplot(111)
    a.plot(x, f, label = "Начальная температура", color='blue')
    a.plot(x, y, label = "Конечная температура", color='red')
    a.set_ylabel("T", fontsize=14)
    a.set_xlabel("X", fontsize=14)
    a.legend()
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.get_tk_widget().place(relx=0.2, rely=0.2) #сдвиг по окну
    canvas.draw()
    id = fig.canvas.mpl_connect('key_press_event', continue_) # реакция на нажатие пробела



def continue_(event):
    global green, _x, a, canvas
    a.plot(_x, green, label = "Конечная температура (часть A)", color='green')
    a.legend()
    canvas.draw()




window= Tk()
window.title("Численное решение уравнения теплопроводности")
window.geometry("900x650")

label1 = Label(text = "Длина стержня ", font = "Arial 9")
label2 = Label(text = "Время воздействия ", font = "Arial 9")
label3 = Label(text = "Шаг по t ", font = "Arial 9")
label4 = Label(text = "Шаг по x ", font = "Arial 9")
label6 = Label(text = "b₀ = ", font = "Arial 11")
label7 = Label(text = "b₁ = ", font = "Arial 11")
label8 = Label(text = "b₂ = ", font = "Arial 11")
label9 = Label(text = "ϕ₀ = ", font = "Arial 11")
label10 = Label(text = "ϕ₁ = ", font = "Arial 11")
label5 = Label(text = "Строка прогресса ", font = "Arial 10")

label1.grid(row=0, column=0)
label2.grid(row=1, column=0)
label3.grid(row=1, column=2)
label4.grid(row=0, column=2)
label6.grid(row=3, column=0)
label7.grid(row=3, column=2)
label8.grid(row=3, column=4)
label9.grid(row=2, column=0)
label10.grid(row=2, column=2)
label5.grid(row=1, column=6)

_l=DoubleVar()
_t=DoubleVar()
_h1=DoubleVar()
_h2=DoubleVar()
_b0=DoubleVar()
_b1=DoubleVar()
_b2=DoubleVar()
_f0=DoubleVar()
_f1=DoubleVar()

entry1 = Entry(textvariable = _l, font = "Arial 9")
entry2 = Entry(textvariable = _t, font = "Arial 9")
entry3 = Entry(textvariable = _h1, font = "Arial 9")
entry4 = Entry(textvariable = _h2, font = "Arial 9")
entry5 = Entry(textvariable = _b0, font = "Arial 9")
entry6 = Entry(textvariable = _b1, font = "Arial 9")
entry7 = Entry(textvariable = _b2, font = "Arial 9")
entry8 = Entry(textvariable = _f0, font = "Arial 9")
entry9 = Entry(textvariable = _f1, font = "Arial 9")

entry1.grid(row=0, column=1, padx=5, pady=5)
entry2.grid(row=1, column=1, padx=5, pady=5)
entry3.grid(row=1, column=3, padx=5, pady=5)
entry4.grid(row=0, column=3, padx=5, pady=5)
entry5.grid(row=3, column=1, padx=5, pady=5)
entry6.grid(row=3, column=3, padx=5, pady=5)
entry7.grid(row=3, column=5, padx=5, pady=5)
entry8.grid(row=2, column=1, padx=5, pady=5)
entry9.grid(row=2, column=3, padx=5, pady=5)


entry1.delete(0, END)
entry1.insert(0, "20.0")

entry2.delete(0, END)
entry2.insert(0, "10.0")

entry3.delete(0, END)
entry3.insert(0, "0.01")

entry4.delete(0, END)
entry4.insert(0, "0.2")

entry5.delete(0, END)
entry5.insert(0, "0.001")

entry6.delete(0, END)
entry6.insert(0, "0.03")

entry7.delete(0, END)
entry7.insert(0, "0.0")

entry8.delete(0, END)
entry8.insert(0, "0.0")

entry9.delete(0, END)
entry9.insert(0, "0.0")


mpb = ttk.Progressbar(window, orient ="horizontal",length = 180, mode ="determinate")
mpb.grid(row = 2, column=6)


global green, _x, a, canvas



btn = Button(text= "Расчет", font = "Arial 10", command =lambda: decision(mpb, window, float(_l.get()), float(_t.get()), float(_h1.get()), float(_h2.get()), float(_b0.get()), float(_b1.get()), float(_b2.get()), float(_f0.get()), float(_f1.get())))
btn.grid(row=0,column=6)

window.mainloop()
