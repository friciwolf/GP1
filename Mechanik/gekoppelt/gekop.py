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
lF = np.array([0.527,0.781,0.2786])
st = 0.02
sTs = []
sTsf = []

#---------------------------------------------------
#Gegensinnige Schwingung bei lF=0.527 m
#---------------------------------------------------
t, U1, U2 = c_open("527cm/Gegensinnig.lab")
U1 = np.array(U1) - np.average(U1)
U2 = np.array(U2) - np.average(U2)

peaks1 = lok_max(t,U1)
peaks2 = lok_max(t,U2)

subplot(2,1,1)
plt.title("Gegensinnige Schwingung  $l_F$=52.7cm")
plt.ylabel("$U_1$ / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("$U_2$ / V")
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

sT1 = np.sqrt(2)*st/len(peaks1)
sT2 = np.sqrt(2)*st/len(peaks2)
sTsf.append(0.5*np.sqrt(sT1**2+sT2**2))

#---------------------------------------------------
#Gleichsinnige Schwingung bei lF=0.527 m
#---------------------------------------------------
t, U1, U2 = c_open("527cm/Gleichsinnig.lab")
U1 = np.array(U1) - np.average(U1)
U2 = np.array(U2) - np.average(U2)

peaks1 = lok_max(t,U1)
peaks2 = lok_max(t,U2)

subplot(2,1,1)
plt.title("Gleichsinnige Schwingung  $l_F$=52.7cm")
plt.ylabel("$U_1$ / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("$U_2$ / V")
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

sT1 = np.sqrt(2)*st/len(peaks1)
sT2 = np.sqrt(2)*st/len(peaks2)
sTs.append(0.5*np.sqrt(sT1**2+sT2**2))


#---------------------------------------------------
#Gegensinnige Schwingung bei lF=0.781 m
#---------------------------------------------------
t, U1, U2 = c_open("781cm/Gegensinnig.lab")
U1 = np.array(U1) - np.average(U1)
U2 = np.array(U2) - np.average(U2)

peaks1 = lok_max(t,U1)
peaks2 = lok_max(t,U2)

subplot(2,1,1)
plt.title("Gegensinnige Schwingung  $l_F$=78.1cm")
plt.ylabel("$U_1$ / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("$U_2$ / V")
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

sT1 = np.sqrt(2)*st/len(peaks1)
sT2 = np.sqrt(2)*st/len(peaks2)
sTsf.append(0.5*np.sqrt(sT1**2+sT2**2))

#---------------------------------------------------
#Gleichsinnige Schwingung bei lF=0.781 m
#---------------------------------------------------
t, U1, U2 = c_open("781cm/Gleichsinnig.lab")
U1 = np.array(U1) - np.average(U1)
U2 = np.array(U2) - np.average(U2)

peaks1 = lok_max(t,U1)
peaks2 = lok_max(t,U2)

subplot(2,1,1)
plt.title("Gleichsinnige Schwingung  $l_F$=78.1cm")
plt.ylabel("$U_1$ / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("$U_2$ / V")
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

sT1 = np.sqrt(2)*st/len(peaks1)
sT2 = np.sqrt(2)*st/len(peaks2)
sTs.append(0.5*np.sqrt(sT1**2+sT2**2))

#---------------------------------------------------
#Gegensinnige Schwingung bei lF=0.2786m
#---------------------------------------------------
t, U1, U2 = c_open("2786cm/Gegensinnig1.lab")
U1 = np.array(U1) - np.average(U1)
U2 = np.array(U2) - np.average(U2)

peaks1 = lok_max(t,U1)
peaks2 = lok_max(t,U2)

subplot(2,1,1)
plt.title("Gegensinnige Schwingung  $l_F$=27.86cm")
plt.ylabel("$U_1$ / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("$U_2$ / V")
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

sT1 = np.sqrt(2)*st/len(peaks1)
sT2 = np.sqrt(2)*st/len(peaks2)
sTsf.append(0.5*np.sqrt(sT1**2+sT2**2))

#---------------------------------------------------
#Gleichsinnige Schwingung bei lF=0.2786m - 1A
#---------------------------------------------------
t, U1, U2 = c_open("2786cm/Gleichsinnig1A.lab")
U1 = np.array(U1) - np.average(U1)
U2 = np.array(U2) - np.average(U2)

peaks1 = lok_max(t,U1)
peaks2 = lok_max(t,U2)

subplot(2,1,1)
plt.title("Gleichsinnige Schwingung  $l_F$=27.86cm (1A)")
plt.ylabel("$U_1$ / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("$U_2$ / V")
plt.xlabel("t / s")
plt.plot(t, U2)
scatter(peaks2[0], peaks2[1], color="r")
plt.savefig("images/2786_gleich1A.png")
plt.show()
plt.close()

T1 = calc_T(peaks1)
T2 = calc_T(peaks2)
T11 = (T1+T2)*0.5

sT1 = np.sqrt(2)*st/len(peaks1)
sT2 = np.sqrt(2)*st/len(peaks2)
sT11 = 0.5*np.sqrt(sT1**2+sT2**2)

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
plt.title("Gleichsinnige Schwingung  $l_F$=27.86cm (1B)")
plt.ylabel("$U_1$ / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("$U_2$ / V")
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

sT1 = np.sqrt(2)*st/len(peaks1)
sT2 = np.sqrt(2)*st/len(peaks2)
sT22 = 0.5*np.sqrt(sT1**2+sT2**2)
sTs.append(0.5*np.sqrt(sT11**2+sT22**2))


#---------------------------------------------------
#Schwebung bei lF=0.527 m 
#---------------------------------------------------

t, U1, U2 = c_open("527cm/Schwebung.lab")
U1 = np.array(U1) - np.average(U1)
U2 = np.array(U2) - np.average(U2)

peaks1 = lok_max(t,U1)
peaks2 = lok_max(t,U2)

subplot(2,1,1)
plt.title("Schwebung bei  $l_F$=52.7cm")
plt.ylabel("$U_1$ / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("$U_2$ / V")
plt.xlabel("t / s")
plt.plot(t, U2)
scatter(peaks2[0], peaks2[1], color="r")
plt.savefig("images/527_schwe.png")
plt.show()
plt.close()

omega_fft1, A1 = anal.fourier_fft(t, U1)
omega_fft2, A2 = anal.fourier_fft(t, U2)
plt.title("Fouriertransformierte der Schwebungen bei  $l_F$=52.7cm")
plt.ylabel("Häufigkeit")
plt.xlabel("f / Hz")
plt.xlim(0.4,0.8)
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
plt.title("Schwebung bei  $l_F$=78.1cm")
plt.ylabel("$U_1$ / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("$U_2$ / V")
plt.xlabel("t / s")
plt.plot(t, U2)
scatter(peaks2[0], peaks2[1], color="r")
plt.savefig("images/781_schwe.png")
plt.show()
plt.close()

omega_fft1, A1 = anal.fourier_fft(t, U1)
omega_fft2, A2 = anal.fourier_fft(t, U2)
plt.title("Fouriertransformierte der Schwebungen bei  $l_F$=78.1cm")
plt.ylabel("Häufigkeit")
plt.xlabel("f / Hz")
plt.xlim(0.4,0.8)
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
plt.title("Schwebung bei  $l_F$=27.86cm")
plt.ylabel("$U_1$ / V")
plt.xlabel("t / s")
plt.plot(t, U1)
scatter(peaks1[0], peaks1[1], color="r")

subplot(2,1,2)
plt.ylabel("$U_2$ / V")
plt.xlabel("t / s")
plt.plot(t, U2)
scatter(peaks2[0], peaks2[1], color="r")
plt.savefig("images/2786_schwe.png")
plt.show()
plt.close()

omega_fft1, A1 = anal.fourier_fft(t, U1)
omega_fft2, A2 = anal.fourier_fft(t, U2)
plt.title("Fouriertransformierte der Schwebungen bei  $l_F$=27.86cm")
plt.ylabel("Häufigkeit")
plt.xlabel("f / Hz")
plt.xlim(0.4,0.8)
plt.plot(omega_fft1, A1)
plt.plot(omega_fft2, A2)
plt.show()

#---------------------------------------------------
#Kappa Berechnen - über die Periodendauer T
#---------------------------------------------------
wsf = np.array(wsf)
ws = np.array(ws)
print(wsf)
print(ws)
sTs = np.array(sTs)
sTsf = np.array(sTsf)
kappa = np.array((wsf**2-ws**2)/(wsf**2+ws**2))
print(kappa)

skappa = 4*np.power(ws*wsf, 3)*(scipy.pi*2)**(-1)*np.power(ws**2+wsf**2, -2)*np.sqrt((sTs/ws)**2+(sTsf/wsf)**2)

print(skappa)
m,em,b,eb,chiq,corr = anal.lineare_regression(lF**-2, np.array(kappa)**(-1)-1, skappa)
print(b)
print(eb)
plt.errorbar(lF**-2, kappa**-1, skappa/(kappa)**2,fmt="x",markersize=4, capsize=3)
plt.plot(np.arange(min(lF**-2)-0.2, max(lF**-2)+0.2, 0.1), np.arange(min(lF**-2)-0.2, max(lF**-2)+0.2, 0.1)*m+b)
plt.title("Die lineare Regression bzgl. $\kappa^{-1}$")
plt.xlabel("$l_F^{-2}$")
plt.ylabel("$\kappa^{-1}$ / $m^{-2}$")
plt.savefig("images/D_linreg.png")
plt.show()

plt.title("Residuenplot")
plt.errorbar(lF**-2, kappa**-1-lF**-2*m-b, skappa/(kappa)**2,fmt="x",markersize=4, capsize=3)
plt.axhline(0)
plt.xlabel("$l_F^{-2}$")
plt.ylabel("$\kappa^{-1}$ / $m^{-2}$")
plt.savefig("images/D_residum.png")
plt.show()

print(M*ls*9.81/m)

#---------------------------------------------------
#Kappa Berechnen - über die FFT
#---------------------------------------------------
print("fft:")
omegasf_fft = np.array([0.6043956044, 0.6715506716, 0.5524356])*2*scipy.pi
omegas_fft = np.array([0.5311355311, 0.5311355311, 0.5310706873])*2*scipy.pi
print(wsf)
print(omegasf_fft)
print(ws)
print(omegas_fft)

somegasf_fft = 0.0061*np.ones(len(omegasf_fft))
somegas_fft = 0.0061*np.ones(len(omegasf_fft))
print("kappa")
print(kappa)
kappa = np.array((omegasf_fft**2-omegas_fft**2)/(omegasf_fft**2+omegas_fft**2))
print(kappa)
print(skappa)
skappa = 4*omegasf_fft*omegas_fft/(omegasf_fft**2+omegas_fft**2)**2 * np.sqrt((omegasf_fft*somegas_fft)**2+(omegas_fft*somegasf_fft)**2)

print(skappa)
m,em,b,eb,chiq,corr = anal.lineare_regression(lF**-2, np.array(kappa)**(-1)-1, skappa)
print(b)
print(eb)
plt.errorbar(lF**-2, kappa**-1, skappa/(kappa)**2,fmt="x",markersize=4, capsize=3)
plt.plot(np.arange(min(lF**-2)-0.2, max(lF**-2)+0.2, 0.1), np.arange(min(lF**-2)-0.2, max(lF**-2)+0.2, 0.1)*m+b)
plt.title("Die lineare Regression bzgl. $\kappa^{-1}$")
plt.xlabel("$l_F^{-2}$")
plt.ylabel("$\kappa^{-1}$ / $m^{-2}$")
plt.savefig("images/Dfft_linreg.png")
plt.show()

plt.title("Residuenplot")
plt.errorbar(lF**-2, kappa**-1-lF**-2*m-b, skappa/(kappa)**2,fmt="x",markersize=4, capsize=3)
plt.axhline(0)
plt.xlabel("$l_F^{-2}$")
plt.ylabel("$\kappa^{-1}$ / $m^{-2}$")
plt.savefig("images/Dfft_residum.png")
plt.show()

print(M*ls*9.81/m)


L0 = 0.270
sL = 0.001
sL0sys = 0.0003+0.0002*0.270
l1 = 0.380
sL1sys = 0.0003+0.0002*0.380
l2 = 0.387
sL2sys = 0.0003+0.0002*0.387
m1 = 0.0654
m2 = 0.0678
sm = 0.00005

D1 = 9.81*(0.0654/(0.38-0.27))
sD1 = D1 * np.sqrt((sm/m)**2+(sL/(l1-L0))**2+(sL/(l1-L0))**2)
sD1s = D1 * np.sqrt((sm/m)**2+(sL0sys/(l1-L0))**2+(sL1sys/(l1-L0))**2)
D2 = 9.81*(0.0678/(0.387-0.270))
sD2 = D2 * np.sqrt((sm/m)**2+(sL/(l2-L0))**2+(sL/(l2-L0))**2)
sD2s = D2 * np.sqrt((sm/m)**2+(sL0sys/(l2-L0))**2+(sL2sys/(l2-L0))**2)
print((D1/sD1**2+D2/sD2**2)/(sD1**-2+sD2**-2))
print(1/(sD1**-2+sD2**-2))
print(1/(sD1s**-2+sD2s**-2))