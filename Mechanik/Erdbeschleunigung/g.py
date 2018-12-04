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
somega = [] #Array der Omega-Fehler der Messungen mit der Stange
lp = 0.67817 #Pendellänge in m
slp = 0.001 #Ableseunsicherheit des Maßbands - größte Unsicherheit
rp = 0.04 #Radius des Pendelkörpers in m
srp = 0.01 /np.sqrt(12)
dt = 0.02
sslp = 0.7 *pow(10, -3) #Kalibrierungssicherheit (systematisch)
ssrp = 0.01/np.sqrt(12)
def lok_max(t, U):
    """
    Berechnet die lokalen Maxima einer Schwingung
    """
    Tmax, m = [], []
    I = 7 #Umgeungsintervall
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
    returns: t, U
    """
    data = cassy.CassyDaten(file)
    t = data.messung(1).datenreihe("t").werte
    U = data.messung(1).datenreihe("U_A1").werte
    return t, U

#---------------------------------------------------
#Frequenzvergleich der Stangen
#---------------------------------------------------
data = cassy.CassyDaten("Frequenzvergleich Stangen.lab")
t = data.messung(1).datenreihe("t").werte
U1 = -data.messung(1).datenreihe("U_A1").werte
U2 = data.messung(1).datenreihe("U_B1").werte

U1 = U1[:-1]
U1 = np.array(U1)-np.average(U1)
U2 = np.array(U2)-np.average(U2)

plt.title("Schwingungen der Stangen 1 und 2 mit der Masse auf 2")
plt.ylabel("U / V")
plt.xlabel("t / s")
plt.plot(t[:-1], U1)
plt.plot(t, U2)
peaks1 = lok_max(t,U1)
peaks2 = lok_max(t,U2)
scatter(peaks1[0], peaks1[1], color="b")
scatter(peaks2[0], peaks2[1], color="orange")

plt.savefig("Images/Frequenzvergleich.pdf")
plt.show()

#T und omega:
T1 = 0.0
dT1 = 0.0
T2 = 0.0
dT2 = 0.0
for i in range(len(peaks1[0])-1):
    T1 += (peaks1[0][i+1]-peaks1[0][i])/(len(peaks1[0])-1)
for i in range(len(peaks2[0])-1):
    T2 += (peaks2[0][i+1]-peaks2[0][i])/(len(peaks2[0])-1)
dT1 = np.sqrt(2) * dt /len(peaks1[0]) #TODO: check
dT2 = np.sqrt(2) * dt /len(peaks2[0]) #TODO: check
print("Schwingung der Stangen:")
print("Die Periodendauer, die Frequenz und die Kreisfrequenz betragen somit:")
print("T1=",T1, "+-", dT1, "s,", "f1 =", 1/T1,"+-",1/(T1**2)*dT1, "Hz,", "w1 =", scipy.pi*2/T1,"+-" ,scipy.pi*2/(T1**2)*dT1 ,"Hz")
print("T2=",T2, "+-", dT2, "s,", "f2 =", 1/T2,"+-",1/(T2**2)*dT2, "Hz,", "w2 =", scipy.pi*2/T2,"+-" ,scipy.pi*2/(T2**2)*dT2 ,"Hz")
#print("T2=",T2, "s,", "f2 =", 1/T2, "Hz,", "w2 =", scipy.pi*2/T2, "Hz")
print("relativer Unterschied:", (T1-T2)/T2)


#---------------------------------------------------
#Alleinige Schwingung der Stange 1
#---------------------------------------------------
t, U = c_open("Messung_Stange1.lab")

U = np.array(U) - np.average(U)

subplot(2,1,2)
plt.ylabel("U / V")
plt.xlabel("t / s")
plt.plot(t[1000:2300], U[1000:2300])
peaks = lok_max(t[1000:2300],U[1000:2300])
scatter(peaks[0], peaks[1], color="r")

subplot(2,1,1)
plt.title("Schwingung der Stange 1 ohne Masse")
plt.ylabel("U / V")
plt.xlabel("t / s")
plt.plot(t, U)
peaks = lok_max(t,U)
scatter(peaks[0], peaks[1], color="r")

plt.savefig("Images/Stange1_oM.pdf")
plt.show()

#T und omega:
T = 0.0
dT = 0.0
for i in range(len(peaks[0])-1):
    T += (peaks[0][i+1]-peaks[0][i])/(len(peaks[0])-1)
dT = np.sqrt(2) * dt /len(peaks1[0])
print("Schwingung der Stange 1:")
print("Die Periodendauer, die Frequenz und die Kreisfrequenz betragen somit:")
#print("T=",T, "s,", "f =", 1/T, "Hz,", "w =", scipy.pi*2/T, "Hz")
print("T=",T, "+-", dT, "s,", "f =", 1/T,"+-",1/(T**2)*dT, "Hz,", "w1 =", scipy.pi*2/T,"+-" ,scipy.pi*2/(T**2)*dT ,"Hz")

#---------------------------------------------------
#Schwingung der Stange 1 mit der Masse darauf
#---------------------------------------------------
t, U = c_open("Stange_Masse1.lab")
t = t[:-1]
U = U[:-1]
U = np.array(U) - np.average(U)

subplot(2,1,2)
plt.ylabel("U / V")
plt.xlabel("t / s")
plt.plot(t[400:1200], U[400:1200])
peaks = lok_max(t[400:1200],U[400:1200])
scatter(peaks[0], peaks[1], color="r")

subplot(2,1,1)
plt.plot(t, U)
peaks = lok_max(t,U)
scatter(peaks[0], peaks[1], color="r")
plt.title("Schwingung der Stange 1 mit der Masse")
plt.ylabel("U / V")
plt.xlabel("t / s")

plt.savefig("Images/Stange1_mM.pdf")
plt.show()

#T und omega:
T = 0.0
dT = 0.0
for i in range(len(peaks[0])-1):
    T += (peaks[0][i+1]-peaks[0][i])/(len(peaks[0])-1)
dT = np.sqrt(2) * dt /len(peaks1[0])
print("Schwingung der Stange 1 mit der Masse:")
print("Die Periodendauer, die Frequenz und die Kreisfrequenz betragen somit:")
#print("T=",T, "s,", "f =", 1/T, "Hz,", "w =", scipy.pi*2/T, "Hz")
print("T=",T, "+-", dT, "s,", "f =", 1/T,"+-",1/(T**2)*dT, "Hz,", "w1 =", scipy.pi*2/T,"+-" ,scipy.pi*2/(T**2)*dT ,"Hz")
omega.append(scipy.pi*2/T)
somega.append(scipy.pi*2/(T**2)*dT)
#---------------------------------------------------
#Schwingung der Stange 2 mit der Masse darauf - beste Messung
#---------------------------------------------------
t, U = c_open("Stasnge_Masse2_beste.lab")
U = np.array(U) - np.average(U)

subplot(2,1,2)
plt.ylabel("U / V")
plt.xlabel("t / s")
plt.plot(t[400:1200], U[400:1200])
peaks = lok_max(t[400:1200],U[400:1200])
scatter(peaks[0], peaks[1], color="r")

subplot(2,1,1)
plt.plot(t, U)
peaks = lok_max(t,U)
scatter(peaks[0], peaks[1], color="r")
plt.title("Schwingung der Stange 1 mit der Masse")
plt.ylabel("U / V")
plt.xlabel("t / s")

plt.savefig("Images/Stange1_mM-b.pdf")
plt.show()

#T und omega:
T = 0.0
sT = 0.0
for i in range(len(peaks[0])-1):
    T += (peaks[0][i+1]-peaks[0][i])/(len(peaks[0])-1)
dT = np.sqrt(2) * dt /len(peaks1[0])
print("Schwingung der Stange 2 mit der Masse: - beste Resultat")
print("Die Periodendauer, die Frequenz und die Kreisfrequenz betragen somit:")
#print("T=",T, "s,", "f =", 1/T, "Hz,", "w =", scipy.pi*2/T, "Hz")
print("T=",T, "+-", dT, "s,", "f =", 1/T,"+-",1/(T**2)*dT, "Hz,", "w1 =", scipy.pi*2/T,"+-" ,scipy.pi*2/(T**2)*dT ,"Hz")
omega.append(scipy.pi*2/T)
somega.append(scipy.pi*2/(T**2)*dT)

#TODO gew. Mittel
omega_aver = (omega[0]/(somega[0])**2+omega[1]/(somega[1])**2)/(np.power(somega[0], -2)+np.power(somega[1], -2))
somega_aver = 1.0/(np.power(somega[0], -2)+np.power(somega[1], -2))
g = omega_aver**2 * lp * (1+0.5*(rp/lp)**2)
print("somega", somega_aver*2*g/omega_aver)
print("sl", (omega_aver**2-0.5*(rp/lp)**2)*slp)
print("sr", srp * omega_aver**2*(rp/lp))
sg = np.sqrt((somega_aver*2*g/omega_aver)**2+((omega_aver**2-0.5*(rp/lp)**2)*slp)**2+(srp * omega_aver**2*(rp/lp))**2)
ssg = np.sqrt(((omega_aver**2-0.5*(rp/lp)**2)*sslp)**2+(ssrp * omega_aver**2*(rp/lp))**2)
print("g =", g, "+-", sg, "+-", ssg)