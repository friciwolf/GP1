#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 13:27:50 2018

@author: Mate
"""
import matplotlib.pyplot as plt
import praktikum.cassy as cassy
from scipy.signal import find_peaks_cwt
import numpy as np
from pylab import *
import scipy.constants

#---------------------------------------------------
#Definitionen
#---------------------------------------------------
omega = [] #Array der Omega-Werte der Messungen mit der Stange
lp = 0.67817 #Pendellänge in m
rp = 0.04 #Radius des Pendelkörpers in m
def lok_max(t, U):
    """
    Berechnet die lokalen Maxima einer Schwingung
    """
    Tmax, m = [], []
    for i in range(4, len(U)-4):
        maximum = True
        for j in range(i-4, i+4):
            if U[j]>U[i]:
                maximum = False
                break
        if maximum:
            max_wahr = True
            for k in range(len(m)):
                #if np.round(m[k], 1)==np.round(U[i],2):
                if np.abs(Tmax[k]-t[i])<0.5:
                    max_wahr = False
            if max_wahr:
                Tmax.append(t[i])
                m.append(U[i])
    return Tmax, m

def c_open(file):
    """
    Zum Einlesen der CASSY-Messungen
    Parameter: file - Messdatei
    returns: t, U
    """
    data = cassy.CassyDaten(file)
    t = data.messung(1).datenreihe("t").werte
    U = data.messung(1).datenreihe("U_A1").werte
    return t, U

#---------------------------------------------------
#Alleinige Schwingung der Stange 1
#---------------------------------------------------
t, U = c_open("Messung_Stange1.lab")
peaks = lok_max(t,U)

plt.title("Schwingung der Stange 1 ohne Masse")
plt.ylabel("U / V")
plt.xlabel("t / s")
plt.plot(t, U)
scatter(peaks[0], peaks[1], color="r")
plt.show()

#T und omega:
T = 0.0
for i in range(len(peaks[0])-1):
    T += (peaks[0][i+1]-peaks[0][i])/(len(peaks[0])-1)
print("Schwingung der Stange 1:")
print("Die Periodendauer, die Frequenz und die Kreisfrequenz betragen somit:")
print("T=",T, "s,", "f =", 1/T, "Hz,", "w =", scipy.pi*2/T, "Hz")

#---------------------------------------------------
#Schwingung der Stange 1 mit der Masse darauf
#---------------------------------------------------
t, U = c_open("Stange_Masse1.lab")
peaks = lok_max(t,U)

plt.title("Schwingung der Stange 1 mit der Masse")
plt.ylabel("U / V")
plt.xlabel("t / s")
plt.plot(t, U)
scatter(peaks[0], peaks[1], color="r")
plt.show()

#T und omega:
T = 0.0
for i in range(len(peaks[0])-1):
    T += (peaks[0][i+1]-peaks[0][i])/(len(peaks[0])-1)
print("Schwingung der Stange 1 mit der Masse:")
print("Die Periodendauer, die Frequenz und die Kreisfrequenz betragen somit:")
print("T=",T, "s,", "f =", 1/T, "Hz,", "w =", scipy.pi*2/T, "Hz")
omega.append(scipy.pi*2/T)
#---------------------------------------------------
#Schwingung der Stange 2 mit der Masse darauf - beste Messung
#---------------------------------------------------
t, U = c_open("Stasnge_Masse2_beste.lab")
peaks = lok_max(t,U)
plt.plot(t, U)
scatter(peaks[0], peaks[1], color="r")
plt.title("Schwingung der Stange 1 mit der Masse")
plt.ylabel("U / V")
plt.xlabel("t / s")
plt.show()

#T und omega:
T = 0.0
for i in range(len(peaks[0])-1):
    T += (peaks[0][i+1]-peaks[0][i])/(len(peaks[0])-1)
print("Schwingung der Stange 2 mit der Masse: - beste Resultat")
print("Die Periodendauer, die Frequenz und die Kreisfrequenz betragen somit:")
print("T=",T, "s,", "f =", 1/T, "Hz,", "w =", scipy.pi*2/T, "Hz")
omega.append(scipy.pi*2/T)

omega_average = np.mean(omega)
g = omega_average**2 * lp * (1+0.5*(rp/lp)**2)
print("g =", g)