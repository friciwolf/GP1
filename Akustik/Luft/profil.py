#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 11:01:47 2019

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
    R = data.messung(1).datenreihe("R_B1").werte
    U = data.messung(1).datenreihe("U_A1").werte
    return np.array(R), np.array(U)

def finde_Knoten(R, U):
    R_kn = []
    I = 30
    for i in range(I, len(R)-I):
        minumum = True
        for i2 in range(i-I, i+I):
            if U[i2]<U[i]:
                minumum=False
        if (minumum and R[i] not in R_kn): R_kn.append(R[i])
    return np.array(R_kn)

f = 1608.56165094
s_f = 0.67901543
s_R  = 0.005/np.sqrt(12)


data = cassy.CassyDaten("Kalibrierung_Wegaufnehmer.lab")
s_k = []
R_1 = []
s_k = data.messung(1).datenreihe("s").werte
R_k = data.messung(1).datenreihe("R_B1").werte
s_1 = np.array(s_k)-np.mean(s_k)
R_1 = np.array(R_k)-np.mean(R_k)
a, ea, b, eb, chiq, corr = anal.lineare_regression_xy(np.array(s_1), np.array(R_1), np.ones(len(s_1))*0.1/np.sqrt(12),np.ones(len(s_1))*0.005/np.sqrt(12))

def R_to_s(R):
    return (R-b-np.mean(R_k)+np.mean(s_k)*a)/a

R, U = c_open("V2.1C.lab")
#finde die Knoten
R_peak = finde_Knoten(R, U)
s_peak = R_to_s(R_peak)
s_s_peak = s_R/a * np.ones(len(s_peak))
s_sys_s = np.sqrt((R_peak/a**2)*ea**2) #ist der sys. Fehler auf R bekannt?
for r in R_peak:
    plt.axvline(r)
plt.plot(R, U, ls="None", marker="x", markersize=2)
plt.xlabel("R / $\Omega$")
plt.ylabel("U / V")
plt.title("Schwingungsprofil des Rohres")
plt.show()
print("c =",(s_peak[0]-s_peak[1])*2*f, "+-", np.sqrt((2*f)**2*(s_s_peak[0]**2+s_s_peak[1]**2)+((s_peak[0]-s_peak[1])*2)**2*s_f**2), "+-", np.sqrt((2*f)**2*(s_sys_s[0]**2+s_sys_s[1]**2)))