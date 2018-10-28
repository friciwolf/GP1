#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 27 16:08:54 2018

@author: Mate
<<<<<<< HEAD
=======
<<<<<<< HEAD
=======
>>>>>>> 8791f053e71a070377bb6fcd68c1266330aa5963

hello world - asdgasd
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy
from praktikum import analyse as anal
import praktikum.cassy as cassy

#Fehler Beachten+Daraufaddieren!!!!!
R=8.314
N=10000 #Anzahl der Elemente die zur Berechnung von L betrachtet werden = max. obere Temperatur

#Kalibrierungsdaten
#Aus dem Resultat von Kalib.py
def T_real(T_gemessen):
    return 0.9727951049804687*(T_gemessen-50) + 1.126251220703125 + 50

t, p, T = [], [], [] #Temperatur in °C!!!
data = cassy.CassyDaten("Dampfdruckkurve.lab")
t = data.messung(1).datenreihe("t").werte
T = data.messung(1).datenreihe("&J_A11").werte
p = data.messung(1).datenreihe("p_B1").werte

#Mitbetrachtung der Kalibrierung:
T = T_real(T)
T = T+273.15
#Plots

plt.plot(t, T-273.15)
plt.title("Abkühlung des Kolbens")
plt.ylabel("Temperatur in K")
plt.xlabel("Zeit in s")
plt.show()

plt.plot(t, p)
plt.title("Dampfdruckänderung im Kolben")
plt.ylabel("Druck in Pa")
plt.xlabel("Zeit in s")
plt.show()

L_Werte=[[],[]]
T_Werte=[[],[]]
for i in range(1,11):
    xi=np.power(T[:1000*i], -1)-np.power(T[0], -1)
    yi=(-np.log(p)[:1000*i]+np.log(p[0]))*R
    L,eL,b,eb,chiq,corr = anal.lineare_regression_xy(xi, yi, np.ones(len(xi)), np.ones(len(yi)))
    L_Werte[0].append(L)
    L_Werte[1].append(eL)
    T_Werte[0].append(T[1000*i])
    T_Werte[1].append(1)

L,eL,b,eb,chiq,corr = anal.lineare_regression_xy(np.array(T_Werte[0]), np.array(L_Werte[0]), np.array(T_Werte[1]), np.array(L_Werte[1]))
plt.plot(T_Werte[0], L_Werte[0])
x=np.arange(min(T_Werte[0])-2,max(T_Werte[0])+2,0.01)
plt.plot(x, L*x+b)
plt.title("Lamda in Abhängigkeit von T")
plt.show()
print("Tabelle von Werten")
print("L || T")
for i in range(1,11):
    print(str(T_Werte[0][i-1]) + " || " + str(L*T_Werte[0][i-1]+b))

x=np.power(T[:N], -1)-np.power(T[0], -1)
y=(-np.log(p)[:N]+np.log(p[0]))*R
plt.plot(x, y)
plt.xlabel("1/T-1/T_0")
plt.ylabel("-R*ln(p/p_0)")

L,eL,b,eb,chiq,corr = anal.lineare_regression_xy(x, y, np.ones(len(x)), np.ones(len(y)))
plt.plot(x, L*x+b, color='green')
plt.title(str(T[N]))
plt.show()
print('L = (%g +- %g) J/mol,   b = (%g +- %g) J/(mol*K),  chi2/dof = %g / %g  corr = %g' % (L, eL, b, eb, chiq, len(x)-2, corr))
