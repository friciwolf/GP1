#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 16:29:50 2019

@author: Mate
"""

import numpy as np
import scipy as sc
import praktikum.cassy as cassy
import praktikum.analyse as anal
import matplotlib.pyplot as plt

def c_open(file, N_Messung):
    """
    Zum Einlesen der CASSY-Messungen
    Parameter: file - Messdatei
    returns: t, U, I 
    """
    data = cassy.CassyDaten(file)
    s = []
    Dt_A1 = []
    R = []
    for i in range(N_Messung):
        s.append(data.messung(i+1).datenreihe("s").werte)
        Dt_A1.append(data.messung(i+1).datenreihe("&Dt_A1").werte)
        R.append(data.messung(i+1).datenreihe("R_B1").werte)
    return np.array(s), np.array(Dt_A1), np.array(R)

##############################
#Kalibrierung
##############################
data = cassy.CassyDaten("Kalibrierung_Wegaufnehmer.lab")
s_k = []
R_1 = []
s_k = data.messung(1).datenreihe("s").werte
R_k = data.messung(1).datenreihe("R_B1").werte
s_1 = np.array(s_k)-np.mean(s_k)
R_1 = np.array(R_k)-np.mean(R_k)
plt.errorbar(s_1, R_1, yerr = np.ones(len(s_1))*0.005/np.sqrt(12), xerr=np.ones(len(s_1))*0.1/np.sqrt(12), fmt ="x")
a, ea, b, eb, chiq, corr = anal.lineare_regression_xy(np.array(s_1), np.array(R_1), np.ones(len(s_1))*0.1/np.sqrt(12),np.ones(len(s_1))*0.005/np.sqrt(12))
plt.plot(np.arange(-13, 13, 0.01), a*np.arange(-13, 13, 0.01)+b, label=str(a)+"*x+"+str(b))
print("Parameter a=", a, "b=" ,b)
print("Fehler auf a: s_a=", ea, "; s_b=", eb)
print("chiq/Ndf = ", chiq/(len(s_1)-2))
plt.xlabel("s/cm")
plt.ylabel("R/ $\Omega$")
plt.title("Kalibriereung des Wegaufnehmers")
plt.show()

plt.errorbar(s_1,R_1-a*s_1-b, np.sqrt(((np.ones(len(s_1))*0.005/np.sqrt(12)))**2+(a*np.ones(len(s_1))*0.1/np.sqrt(12))**2), fmt="x", capsize=5)
plt.axhline(0)
plt.title("Residuenplot zur Kalibrierung des Wegaufnehmers")
plt.xlabel("s/cm")
plt.ylabel("R / $\Omega$")
plt.show()

##############################
#Einlesen der Rohdaten
##############################
s, D, R = c_open("V2.1A.lab", 6)
sig_R = 0.005/np.sqrt(12)
#plt.plot(s,D)
#plt.show()
s2 = []
D2 = []
R2 = []
sig_D2 = []
for i in range(len(s)):
    s2.append(np.average(s[i])) #falls benötigt
    sig_D2.append(np.std(D[i]))
    D2.append(np.average(D[i]))
    R2.append(np.average(R[i]))
sig_D2 = np.array(sig_D2)
s2 = (R2-b-np.mean(R_k)+np.mean(s_k)*a)/a #Umrechnung der Widerstände in Länge
sig_s2 = sig_R*np.ones(len(R2))/a
sig_sys_s2 = np.sqrt((R2/a**2)*ea**2) #sys. Fehler auf R?
plt.plot(D2,s2, label="Rohdaten")
plt.ylabel("s / cm")
plt.xlabel("t / s")
plt.title("Messung der Schallgeschwindigkeit")
a,ea, b, eb,  chiq, corr = anal.lineare_regression_xy(np.array(D2[:-1]),np.array(s2[:-1]), sig_D2[:-1], sig_s2[:-1])
print(a, "*t+", b)
plt.plot(np.arange(min(D2), max(D2), 0.0001), np.arange(min(D2), max(D2), 0.0001)*a+b, label=str(a)+"*t+"+str(b))
plt.errorbar(D2, np.array(s2), np.sqrt((np.ones(len(s2))*0.1/np.sqrt(12))**2+(a*sig_D2)**2), fmt="o",capsize=5, markersize=5)
plt.legend()
plt.show()

plt.errorbar(D2[:-1], np.array(s2[:-1])-a*np.array(D2[:-1])-b, np.sqrt((np.ones(len(s2[:-1]))*0.1/np.sqrt(12))**2+(a*sig_D2[:-1])**2), fmt="o",capsize=5, markersize=5)
#plt.errorbar(D2[-1], np.array(s2[-1])-a*np.array(D2[-1])-b, np.sqrt((0.1/np.sqrt(12))**2+(a*sig_D2[-1])**2), fmt="o",capsize=5, markersize=5, color="r") #nehme inh raus, der Punkt ist anhand der Daten zu unsicher
plt.axhline(0)
plt.title("Residumplot")    
plt.xlabel("t / s")
plt.ylabel("s /cm")
plt.show()
print("chiq/Ndf = ", chiq/(len(D2[:-1])-2))
print("a = ", a, "+-", ea)
#für den sys. Fehler (Verschiebemethode):
a,ea, b, eb,  chiq, corr = anal.lineare_regression_xy(np.array(D2[:-1]),np.array(s2[:-1])-sig_sys_s2[:-1], sig_D2[:-1], sig_sys_s2[:-1])
a1 = a
a,ea, b, eb,  chiq, corr = anal.lineare_regression_xy(np.array(D2[:-1]),np.array(s2[:-1])+sig_sys_s2[:-1], sig_D2[:-1], sig_sys_s2[:-1])
a2 = a
sa = 0.5*(np.abs(a-a1)+np.abs(a-a2))
print("+-", sa)
print("(Schallgeschwindigkeit)")