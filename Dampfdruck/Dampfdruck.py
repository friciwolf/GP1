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
import Rauschmessung

def T_th(T_gem):
    return (T_gem-Kalib.b-50)/Kalib.a +50
    
#Fehler Beachten+Daraufaddieren!!!!!
R=sc.R
N=10000 #Anzahl der Elemente die zur Berechnung von L betrachtet werden = max. obere Temperatur

t, p, T = [], [], [] #Temperatur in °C!!!
data = cassy.CassyDaten("Dampfdruckkurve.lab")
t = data.messung(1).datenreihe("t").werte
T = data.messung(1).datenreihe("&J_A11").werte
p = data.messung(1).datenreihe("p_B1").werte


#-----------------------------------------------------------------------#
#Auswertung der Dampfdruckkurve und Bestimmung der Verdampfungsenthalpie#
#-----------------------------------------------------------------------#

#------------------------------------------------------------
#Ausführen der Kalibrierung
T_real,sTstat,sTsys=Kalib.T_th(T,Kalib.a,Kalib.b,Rauschmessung.sT,Kalib.ea,Kalib.eb,Kalib.sigmaT_Gerät_gr_70,Kalib.easys,Kalib.ebsys)
T_real = T_real+273.15 #zu Kelvin
#TODO: Fehler

#Plot T
plt.plot(t, T_real, color='red')
plt.title("Abkühlung des Kolbens")
plt.ylabel("Temperatur in K")
plt.xlabel("Zeit in s")
plt.savefig("Images/Dampfdruck_T.pdf")
plt.figure()

#Plot p
plt.plot(t, p)
plt.title("Dampfdruckänderung im Kolben")
plt.ylabel("Druck in hPa")
plt.xlabel("Zeit in s")
plt.savefig("Images/Dampfdruck_p.pdf")
plt.figure()

#------------------------------------------------------------
#Berechnung der Verdampfungsenthalpie
L_Werte=[[],[],[], []] #L_sys, s_sys, L_stat, s_stat
T_Werte=[[],[],[]]
for i in range(1,11):
    xi=np.power(T_real[:1000*i], -1)-np.power(T_real[0], -1)
    yi=(-np.log(p)[:1000*i]+np.log(p[0]))*R
    sxi_sys = np.sqrt(np.power(T_real[:1000*i], -4)*sTsys[:1000*i]**2+ np.power(T_real[0], -4)*sTsys[0])
    syi_sys = R*Kalib.sigmap*np.sqrt(np.power(p[:1000*i], -2)+np.power(p[0], -2))
    Lsys,eLsys,bsys,ebsys,chiqsys,corr = ana.lineare_regression_xy(xi, yi, sxi_sys, syi_sys)
    sxi_stat = np.sqrt(np.power(T_real[:1000*i], -4)*sTstat[:1000*i]**2+ np.power(T_real[0], -4)*sTstat[0])
    syi_stat = R*Kalib.sigmap*np.sqrt(np.power(p[:1000*i], -2)+np.power(p[0], -2))
    Lstat,eLstat,bstat,ebstat,chiqstat,corr = ana.lineare_regression_xy(xi, yi, sxi_stat, syi_stat)
    L_Werte[0].append(Lsys)
    L_Werte[1].append(eLsys)
    L_Werte[2].append(Lstat)
    L_Werte[3].append(eLstat)
    T_Werte[0].append(T_real[1000*i])
    T_Werte[1].append(sTstat[1000*i])
    T_Werte[2].append(sTsys[1000*i])

#wichte L nach den Fehlern
L_gew = (L_Werte[0]/np.power(L_Werte[1], 2) + L_Werte[2]/np.power(L_Werte[3], 2))/(np.power(L_Werte[1], -2) + np.power(L_Werte[3], -2))

#Plot der Verdampfungsenthalpie
L,eL,b,eb,chiq,corr = ana.lineare_regression_xy(np.array(T_Werte[0]), np.array(L_gew), np.array(T_Werte[1]), np.array(L_Werte[1]))
plt.plot(T_Werte[0], L_Werte[0])
plt.errorbar(T_Werte[0], L_gew, L_Werte[1], T_Werte[2], ecolor="red", fmt='or', markersize=4, capsize=3)
x=np.arange(min(T_Werte[0])-2,max(T_Werte[0])+2,0.01)
plt.plot(x, L*x+b)
plt.title("Verdampfungsenthalpie in Abhängigkeit von T")
plt.xlabel('Temperatur in K')
plt.ylabel('Enthalpie in J/mol')
plt.savefig('Images/Dampfdruck_L.pdf')
plt.figure()

#Wertetabelle
print("\nTabelle von Werten")
print("    T     ||        L   ")
for i in range(1,11):
    print('{:.5f} || {}'.format(T_Werte[0][i-1],L*T_Werte[0][i-1]+b))

#Für die Clausius-Clapeyron-Ausgleichsgerade
x=np.power(T_real[:N], -1)-np.power(T_real[0], -1)
sx = np.sqrt(np.power(T_real[:N], -4)*sTsys[:N]**2+ np.power(T_real[0], -4)*sTsys[0])
y=(-np.log(p)[:N]+np.log(p[0]))*R
sy = R*Kalib.sigmap*np.sqrt(np.power(p[:N], -2)+np.power(p[0], -2))
plt.plot(x, y)
plt.xlabel("$1/T-1/T_0$")
plt.ylabel("$-R*\ln(p/p_0)$")

L,eL,b,eb,chiq,corr = ana.lineare_regression_xy(x, y, sx, sy)
plt.plot(x, L*x+b, color='green')
plt.title("Clausius-Clapeyron-Ausgleichsgerade")
plt.savefig("Images/Dampfdruck_ln.pdf")
print('\nL = (%g +- %g) J/mol,   b = (%g +- %g) J/(mol*K),  chi2/dof = %g / %g  corr = %g' % (L, eL, b, eb, chiq, len(x)-2, corr))
plt.show()