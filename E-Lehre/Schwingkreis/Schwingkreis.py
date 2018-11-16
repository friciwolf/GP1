#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 08:35:53 2018

@author: Mate
"""
import matplotlib.pyplot as plt
import numpy as np
import praktikum.analyse as anal
import praktikum.cassy as cassy

def f_open(loc):
    """
    Zum Einlesen der Dateien der Oszilloskop
    Parameter: loc - Ort der Datei
    returns: t, U
    """
    file = open(loc)
    t, U = [], []
    for l in file:
        data = l.split(",") # 3<=>t; 4<=>U
        t.append(float(data[3]))
        U.append(float(data[4]))
    return t, U

def c_open(file):
    """
    Zum Einlesen der CASSY-Messungen
    Parameter: file - Messdatei
    returns: t, U, I 
    """
    data = cassy.CassyDaten(file)
    t = data.messung(1).datenreihe("t").werte
    I = data.messung(1).datenreihe("I_A2").werte
    U = data.messung(1).datenreihe("U_B2").werte 
    return t, U, I

def lok_max(t, U):
    """
    Berechnet die lokalen Maxima einer Schwingung
    """
    Tmax, m = [], []
    for i in range(10, len(U)-11):
        if max(np.array(U[i-10:i+10])) == U[i]:
            Tmax.append(t[i])
            m.append(U[i])
    return Tmax, m

#TODO: Entscheidung ob Einhüllende Sinn fürs Plot macht
def delta(t, U):
    """
    Berechnet das Delta einer Schwingung
    """
    delta = 0.0
    tmaxima, maxima = lok_max(t, U)
    delta = np.log(maxima[1]/maxima[0])/(tmaxima[1]-tmaxima[0])
    return delta

#Plot der Messungen mit dem Oszilloskop
for i in range(4):
    t1,t2, U1, U2 = [], [], [], []
    t1, U1 = f_open("ALL000"+str(i)+"/U1.dat")
    t2, U2 = f_open("ALL000"+str(i)+"/U2.dat")
    plt.title("Entladung des Kondensators über den Schwingkreis \n R = ??? $\Omega$")
    plt.plot(t1, U1, label="U1")
    plt.plot(t2, U2, label="U2")
    plt.savefig("Images/O_"+str(i)+".pdf")
    plt.legend()
    plt.show()
    plt.close()

#Plot der Messungen mit den Widerständen R:
R = [5.1, 47, 100, 1000, 208, 173, 137, 100.6]
for i in range(len(R)): 
    f = str(R[i])+"_Messung_"+str(i+1)+".lab"
    t, U, I = c_open(f)
    plt.title("Entladung des Kondensators über den Schwingkreis \n R = " +str(R[i]) + " $\Omega$")
    d = delta(t, U)
    plt.plot(t, U[0]*np.exp(d*np.array(t)), label = "Einhüllende")
    plt.plot(t, I, label="I")
    plt.plot(t, U, label="U")
    plt.savefig("Images/C_"+str(R[i])+".pdf")
    plt.legend()
    plt.show()
    plt.close()

#Grenzfall:
f = "Grenzfall_Messung_9.lab"
t, U, I = c_open(f)
plt.title("Entladung des Kondensators über den Schwingkreis - Grenzfall bei \n R = 111,3 $\Omega$")
plt.plot(t, I, label="I")
plt.plot(t, U, label="U")
plt.savefig("Images/grenzfall.pdf")
plt.legend()
plt.show()
plt.close()