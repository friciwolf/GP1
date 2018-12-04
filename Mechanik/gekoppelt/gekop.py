#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 15:03:13 2018

@author: Mate
"""

import matplotlib.pyplot as plt
import praktikum.cassy as cassy
import praktikum.analyse as anal
from scipy.signal import find_peaks_cwt
import numpy as np
from pylab import *
import scipy.constants

def calc_T(peaks):
    T = 0.0
    for i in range(len(peaks[0])-1):
        T += (peaks[0][i+1]-peaks[0][i])/(len(peaks[0])-1)
    return T

def lok_max(t, U):
    """
    Berechnet die lokalen Maxima einer Schwingung
    """
    Tmax, m = [], []
    I = 20 #Umgeungsintervall
    for i in range(I, len(U)-I):
        maximum = True
        for j in range(i-I, i+I):
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
    returns: t, U1, U2
    """
    data = cassy.CassyDaten(file)
    t = data.messung(1).datenreihe("t").werte
    U1 = data.messung(1).datenreihe("U_A1").werte
    U2 = -data.messung(1).datenreihe("U_B1").werte
    return t, U1, U2

#Omegas
ws = []
wsf = []
W = []
kappa = []
M = 1.0733
ls=1.0

#---------------------------------------------------
#Gegensinnige Schwingung bei lF=0.527 m
#---------------------------------------------------
t, U1, U2 = c_open("527cm/Gegensinnig.lab")
U1 = np.array(U1) - np.average(U1)
U2 = np.array(U2) - np.average(U2)

peaks1 = lok_max(t,U1)
peaks2 = lok_max(t,U2)

subplot(2,1,1)
plt.title("Gegensinnige Schwingung l_F=52.7cm")
plt.ylabel("U1 / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("U2 / V")
plt.xlabel("t / s")
plt.plot(t, U2)
scatter(peaks2[0], peaks2[1], color="r")
plt.savefig("images/527_gegen.png")
plt.show()
plt.close()

T1 = calc_T(peaks1)
T2 = calc_T(peaks2)
T = (T1+T2)*0.5
wsf.append(scipy.pi*2/T)

#---------------------------------------------------
#Gleichsinnige Schwingung bei lF=0.527 m
#---------------------------------------------------
t, U1, U2 = c_open("527cm/Gleichsinnig.lab")
U1 = np.array(U1) - np.average(U1)
U2 = np.array(U2) - np.average(U2)

peaks1 = lok_max(t,U1)
peaks2 = lok_max(t,U2)

subplot(2,1,1)
plt.title("Gleichsinnige Schwingung l_F=52.7cm")
plt.ylabel("U1 / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("U2 / V")
plt.xlabel("t / s")
plt.plot(t, U2)
scatter(peaks2[0], peaks2[1], color="r")
plt.savefig("images/527_gleich.png")
plt.show()
plt.close()

T1 = calc_T(peaks1)
T2 = calc_T(peaks2)
T = (T1+T2)*0.5
ws.append(scipy.pi*2/T)


#---------------------------------------------------
#Gegensinnige Schwingung bei lF=0.781 m
#---------------------------------------------------
t, U1, U2 = c_open("781cm/Gegensinnig.lab")
U1 = np.array(U1) - np.average(U1)
U2 = np.array(U2) - np.average(U2)

peaks1 = lok_max(t,U1)
peaks2 = lok_max(t,U2)

subplot(2,1,1)
plt.title("Gegensinnige Schwingung l_F=78.1cm")
plt.ylabel("U1 / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("U2 / V")
plt.xlabel("t / s")
plt.plot(t, U2)
scatter(peaks2[0], peaks2[1], color="r")
plt.savefig("images/781_gegen.png")
plt.show()
plt.close()

T1 = calc_T(peaks1)
T2 = calc_T(peaks2)
T = (T1+T2)*0.5
wsf.append(scipy.pi*2/T)

#---------------------------------------------------
#Gleichsinnige Schwingung bei lF=0.781 m
#---------------------------------------------------
t, U1, U2 = c_open("781cm/Gleichsinnig.lab")
U1 = np.array(U1) - np.average(U1)
U2 = np.array(U2) - np.average(U2)

peaks1 = lok_max(t,U1)
peaks2 = lok_max(t,U2)

subplot(2,1,1)
plt.title("Gleichsinnige Schwingung l_F=78.1cm")
plt.ylabel("U1 / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("U2 / V")
plt.xlabel("t / s")
plt.plot(t, U2)
scatter(peaks2[0], peaks2[1], color="r")
plt.savefig("images/781_gleich.png")
plt.show()
plt.close()

T1 = calc_T(peaks1)
T2 = calc_T(peaks2)
T = (T1+T2)*0.5
ws.append(scipy.pi*2/T)


#---------------------------------------------------
#Gegensinnige Schwingung bei lF=0.2786m
#---------------------------------------------------
t, U1, U2 = c_open("2786cm/Gegensinnig1.lab")
U1 = np.array(U1) - np.average(U1)
U2 = np.array(U2) - np.average(U2)

peaks1 = lok_max(t,U1)
peaks2 = lok_max(t,U2)

subplot(2,1,1)
plt.title("Gegensinnige Schwingung l_F=27.86cm")
plt.ylabel("U1 / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("U2 / V")
plt.xlabel("t / s")
plt.plot(t, U2)
scatter(peaks2[0], peaks2[1], color="r")
plt.savefig("images/2786_gegen.png")
plt.show()
plt.close()

T1 = calc_T(peaks1)
T2 = calc_T(peaks2)
T = (T1+T2)*0.5
wsf.append(scipy.pi*2/T)

#---------------------------------------------------
#Gleichsinnige Schwingung bei lF=0.2786m - 1A
#---------------------------------------------------
t, U1, U2 = c_open("2786cm/Gleichsinnig1A.lab")
U1 = np.array(U1) - np.average(U1)
U2 = np.array(U2) - np.average(U2)

peaks1 = lok_max(t,U1)
peaks2 = lok_max(t,U2)

subplot(2,1,1)
plt.title("Gleichsinnige Schwingung l_F=27.86cm (1A)")
plt.ylabel("U1 / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("U2 / V")
plt.xlabel("t / s")
plt.plot(t, U2)
scatter(peaks2[0], peaks2[1], color="r")
plt.savefig("images/2786_gleich1A.png")
plt.show()
plt.close()

T1 = calc_T(peaks1)
T2 = calc_T(peaks2)
T11 = (T1+T2)*0.5

#---------------------------------------------------
#Gleichsinnige Schwingung bei lF=0.2786m - 1B
#---------------------------------------------------
t, U1, U2 = c_open("2786cm/Gleichsinnig1B.lab")

U1 = U1[:-1]
U2 = U2[:-1]
t = t[:-1]

U1 = np.array(U1) - np.average(U1)
U2 = np.array(U2) - np.average(U2)

peaks1 = lok_max(t,U1)
peaks2 = lok_max(t,U2)

subplot(2,1,1)
plt.title("Gleichsinnige Schwingung l_F=27.86cm (1B)")
plt.ylabel("U1 / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("U2 / V")
plt.xlabel("t / s")
plt.plot(t, U2)
scatter(peaks2[0], peaks2[1], color="r")
plt.savefig("images/2786_gleich1B.png")
plt.show()
plt.close()

T1 = calc_T(peaks1)
T2 = calc_T(peaks2)
T22 = (T1+T2)*0.5
T = (T11+T22)*0.5
ws.append(scipy.pi*2/T)


#---------------------------------------------------
#Schwebung bei lF=0.527 m 
#---------------------------------------------------

t, U1, U2 = c_open("527cm/Schwebung.lab")
U1 = np.array(U1) - np.average(U1)
U2 = np.array(U2) - np.average(U2)

peaks1 = lok_max(t,U1)
peaks2 = lok_max(t,U2)

subplot(2,1,1)
plt.title("Schwebung bei l_F=52.7cm")
plt.ylabel("U1 / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("U2 / V")
plt.xlabel("t / s")
plt.plot(t, U2)
scatter(peaks2[0], peaks2[1], color="r")
plt.savefig("images/527_schwe.png")
plt.show()
plt.close()

omega_fft1, A1 = anal.fourier_fft(t, U1)
omega_fft2, A2 = anal.fourier_fft(t, U2)
plt.plot(omega_fft1, A1)
plt.plot(omega_fft2, A2)
plt.show()

#---------------------------------------------------
#Schwebung bei lF=0.781 m 
#---------------------------------------------------

t, U1, U2 = c_open("781cm/Schwebung.lab")
U1 = np.array(U1) - np.average(U1)
U2 = np.array(U2) - np.average(U2)

peaks1 = lok_max(t,U1)
peaks2 = lok_max(t,U2)

subplot(2,1,1)
plt.title("Schwebung bei l_F=78.1cm")
plt.ylabel("U1 / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("U2 / V")
plt.xlabel("t / s")
plt.plot(t, U2)
scatter(peaks2[0], peaks2[1], color="r")
plt.savefig("images/781_schwe.png")
plt.show()
plt.close()

omega_fft1, A1 = anal.fourier_fft(t, U1)
omega_fft2, A2 = anal.fourier_fft(t, U2)
plt.plot(omega_fft1, A1)
plt.plot(omega_fft2, A2)
plt.show()

#---------------------------------------------------
#Schwebung bei lF=0.2786 m 
#---------------------------------------------------

t, U1, U2 = c_open("2786cm/Schwebung.lab")
U1 = np.array(U1) - np.average(U1)
U2 = np.array(U2) - np.average(U2)

peaks1 = lok_max(t,U1)
peaks2 = lok_max(t,U2)

subplot(2,1,1)
plt.title("Schwebung bei l_F=27.86cm")
plt.ylabel("U1 / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("U2 / V")
plt.xlabel("t / s")
plt.plot(t, U2)
scatter(peaks2[0], peaks2[1], color="r")
plt.savefig("images/_schwe.png")
plt.show()
plt.close()

omega_fft1, A1 = anal.fourier_fft(t, U1)
omega_fft2, A2 = anal.fourier_fft(t, U2)
plt.plot(omega_fft1, A1)
plt.plot(omega_fft2, A2)
plt.show()

#---------------------------------------------------
#Kappa Berechnen
#---------------------------------------------------
wsf = np.array(wsf)
ws = np.array(ws)
kappa = np.array((wsf**2-ws**2)/(wsf**2+ws**2))
#TODO: Fehler; hier 0.1 nur zufällig!
m,em,b,eb,chiq,corr = anal.lineare_regression(np.array([0.527,0.781,0.2786])**-2, np.array(kappa)**(-1), 0.1*np.ones(len(kappa)))
print(b)
print(M*ls*9.81/m)
print(0.5*9.81*(0.0654/(0.38-0.27)+0.0678/(0.387-0.270)))