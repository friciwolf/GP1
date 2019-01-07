#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 08:35:24 2019

@author: Mate
"""

import numpy as np
import matplotlib.pyplot as plt
import praktikum.analyse as anal
import praktikum.cassy as cassy

def c_open(file):
    """
    Zum Einlesen der CASSY-Messungen
    Parameter: file - Messdatei
    returns: f, U, I 
    """
    data = cassy.CassyDaten(file)
    f = data.messung(1).datenreihe("f_B1").werte
    U = data.messung(1).datenreihe("U_A1").werte
    return np.array(f), np.array(U)

def c_open2(file):
    """
    Zum Einlesen der CASSY-Messungen
    Parameter: file - Messdatei
    returns: f, U, I 
    """
    data = cassy.CassyDaten(file)
    f = data.messung(1).datenreihe("f_A1").werte
    U = data.messung(1).datenreihe("U_B1").werte
    return np.array(f), np.array(U)

def find_max(x,y):
    ymax = 0.0
    index = 0
    for i in range(len(y)):
        if y[i] > ymax:
            index = i
            ymax = y[i]
    return x[index], y[index]

L=0.421
s_L=0.0005
s_L_sys = (0.3+0.2*2)*10**-3
f, U = c_open("V2.1B_richtig.lab")
f2, U2 = c_open2("V2.1B_nochwas.lab")
#f3, U3 = c_open("V2.1B_zuweit.lab")
f = np.append(f, f2)
U = np.append(U, U2)
f_max_0 = np.array([])#Betrachte 5 Peaks, gesucht werden sie mit 4 verschiedenen Methoden, Ergebnis wird danach gemittelt
fmax_ablesen = [410.168, 811.134, 1608.93, 2007, 2415.98]
fmax_peakfinder = np.array([])
fmax_maxima = np.array([])
for i in range(5):
    f_bereich = np.array([])
    U_bereich = np.array([])
    for fr in fmax_ablesen:
        for i_fre in range(len(f)):
            if np.abs(f[i_fre]-fmax_ablesen[i])<=20:
                f_bereich = np.append(f_bereich, f[i_fre])
                U_bereich = np.append(U_bereich, U[i_fre])
    fmax_maxima = np.append(fmax_maxima, find_max(f_bereich, U_bereich)[0])
    fmax_peakfinder = np.append(fmax_peakfinder, anal.peak(f_bereich, U_bereich, min(f_bereich), max(f_bereich)))
for fr in fmax_peakfinder:
    plt.axvline(fr, color="red", ls = "dashed")
for fr in fmax_maxima:
    plt.axvline(fr, color="blue", ls = "dashed")  
for fr in fmax_ablesen:
    plt.axvline(fr, color="green", ls = "dashed")
fmax_matlab = fmax_ablesen #TODO: finde peaks mit Matlab
plt.plot(f, U, linestyle = "None", marker="x", markersize=5)
plt.title("Resonanzfrequenzen")
plt.xlabel("Frequenz / Hz")
plt.ylabel("Spannung / V")
plt.show()
plt.close()
s_f_0 = np.array([])
n_0 = np.array([1, 2, 4, 5, 6])
for i in range(5):
    s_f_0 = np.append(s_f, np.std([fmax_ablesen[i], fmax_peakfinder[i], fmax_maxima[i], fmax_matlab[i]]))
    f_max_0 = np.append(f_max, np.average([fmax_ablesen[i], fmax_peakfinder[i], fmax_maxima[i], fmax_matlab[i]]))
n = n_0[:-1] #ignoriere den letzten Punkt
s_f=s_f_0[:-1]
f_max=f_max_0[:-1]
a,ea, b, eb,  chiq, corr = anal.lineare_regression(n,f_max, s_f)
plt.errorbar(n, f_max, s_f, ls="None", markersize=5, marker="x")
plt.errorbar(n_0[4], f_max_0[4], s_f_0[4], fmt="o", marker="x", markersize=5, color="red")
plt.plot(np.arange(0.5, 6.5, 0.01), a*np.arange(0.5, 6.5, 0.01)+b, label="f = "+str(a)+"*n+"+str(b))
plt.title("Resonanzfrequenzen")
plt.xlabel("Ordnung n")
plt.ylabel("Frequenz / Hz")
plt.legend()
plt.show()

plt.errorbar(n, f_max-n*a-b, s_f, fmt="o",capsize=5, markersize=5)
plt.errorbar(n_0[4], f_max_0[4]-6*a-b, s_f_0[4], fmt="o",capsize=5, markersize=5, color="red")
plt.axhline(0)
plt.title("Residumplot")
plt.xlabel("Ordnung n")
plt.ylabel("Frequenz / Hz")
plt.show()
print("Die Schallgeschwindigkeit ergibt sich daraus auf:")
print("v =", a*2*L, "+-", a*2*L*np.sqrt((ea/a)**2+(s_L/L)**2))
#für die systematische Fehler - Verschiebemethode
a1,ea, b1, eb,  chiq1, corr = anal.lineare_regression(n,f_max-s_f, s_f)
a2,ea, b2, eb,  chiq2, corr = anal.lineare_regression(n,f_max+s_f, s_f)
s_a_sys = (np.abs(a-a1)+np.abs(a-a2))*0.5
print("+-", a*2*L*np.sqrt((s_a_sys/a)**2+(s_L_sys/L)**2))