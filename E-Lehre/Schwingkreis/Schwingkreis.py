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
    t = data.messung(1).datenreihe("t").werte, 
    I = data.messung(1).datenreihe("I_A2").werte
    U = data.messung(1).datenreihe("U_B2").werte, 
    return t, U, I
#Hier kommt noch eine Schleife zur Erstellung der Plots...
t1,t2, U1, U2 = [], [], [], []
t1, U1 = f_open("ALL0000/U1.dat")
t2, U2 = f_open("ALL0000/U2.dat")

plt.title("Entladung des Kondensators durch den Schwingkreis \n $R = 10 \Omega$")
plt.plot(t1, U1, label="U1")
plt.plot(t2, np.array(U2)*10, label="U2 - scale factor: 10")
plt.legend()
plt.show()

f = "100_Messung_3.lab"
t, U, I = c_open(f)
plt.title("Entladung des Kondensators über den Schwingkreis \n R = " + f.split("_")[0]+ " $\Omega$")
plt.plot(t[0], I, label="I")
plt.plot(t[0], U[0], label="U")
plt.legend()
plt.show()