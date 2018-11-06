#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 27 16:08:54 2018

@author: Mate, Patrick
"""
"""
Ausgabe Kalib.py:
    
Damit ergibt sich die Korrekturgleichung:
T_gem = 0.9728*(T_th-50 °C) + 1.13°C + 50°C
oder nach T_th umgestellt:
T_th = (T_gem-1.13-50°C)/0.9728 + 50°C

Somit haben wir T_th(T100) = 100.0 °C +-0.06 °C +-0.52 °C
und T_th(T0) = -0.0 °C +-0.06 °C +-0.38 °C
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
from praktikum import analyse as ana
import praktikum.cassy as cassy
import Kalib

def T_real(T_gem):
    return  (T_gem-Kalib.b-50)/Kalib.a +50
    
#Fehler Beachten+Daraufaddieren!!!!!
R=sc.R
N=10000 #Anzahl der Elemente die zur Berechnung von L betrachtet werden = max. obere Temperatur

t, p, T = [], [], [] #Temperatur in °C!!!
data = cassy.CassyDaten("Dampfdruckkurve.lab")
t = data.messung(1).datenreihe("t").werte
T = data.messung(1).datenreihe("&J_A11").werte
p = data.messung(1).datenreihe("p_B1").werte

#Mitbetrachtung der Kalibrierung:
T = T_real(T) #TODO: Fehler
T = T+273.15
#Plots

plt.plot(t, T, color='red')
plt.title("Abkühlung des Kolbens")
plt.ylabel("Temperatur in K")
plt.xlabel("Zeit in s")
plt.savefig("Images/Dampfdruck_T.pdf")
plt.figure()

plt.plot(t, p)
plt.title("Dampfdruckänderung im Kolben")
plt.ylabel("Druck in hPa")
plt.xlabel("Zeit in s")
plt.savefig("Images/Dampfdruck_p.pdf")
plt.figure()

L_Werte=[[],[]]
T_Werte=[[],[]]
for i in range(1,11):
    xi=np.power(T[:1000*i], -1)-np.power(T[0], -1)
    yi=(-np.log(p)[:1000*i]+np.log(p[0]))*R
    L,eL,b,eb,chiq,corr = ana.lineare_regression_xy(xi, yi, np.ones(len(xi)), np.ones(len(yi)))
    L_Werte[0].append(L)
    L_Werte[1].append(eL)
    T_Werte[0].append(T[1000*i])
    T_Werte[1].append(1) #TODO: Fehler

L,eL,b,eb,chiq,corr = ana.lineare_regression_xy(np.array(T_Werte[0]), np.array(L_Werte[0]), np.array(T_Werte[1]), np.array(L_Werte[1]))
plt.plot(T_Werte[0], L_Werte[0])
x=np.arange(min(T_Werte[0])-2,max(T_Werte[0])+2,0.01)
plt.plot(x, L*x+b)
plt.title("Verdampfungsenthalpie in Abhängigkeit von T")
plt.xlabel('Temperatur in K')
plt.ylabel('Enthalpie in J/mol')
plt.figure()
print("\nTabelle von Werten")
print("    T     ||        L   ")
for i in range(1,11):
    print('{:.5f} || {}'.format(T_Werte[0][i-1],L*T_Werte[0][i-1]+b))

x=np.power(T[:N], -1)-np.power(T[0], -1)
y=(-np.log(p)[:N]+np.log(p[0]))*R
plt.plot(x, y)
plt.xlabel("1/T-1/T_0")
plt.ylabel("-R*ln(p/p_0)")

L,eL,b,eb,chiq,corr = ana.lineare_regression_xy(x, y, np.ones(len(x)), np.ones(len(y)))
plt.plot(x, L*x+b, color='green')
plt.title(str(T[N]))
plt.show()
print('\nL = (%g +- %g) J/mol,   b = (%g +- %g) J/(mol*K),  chi2/dof = %g / %g  corr = %g' % (L, eL, b, eb, chiq, len(x)-2, corr))
