#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 18 11:01:17 2018

@author: Mate
"""

import numpy as np
import praktikum.cassy as cassy

def gew_mittel(x,ex):
    S1 = 0.0
    S2 = 0.0
    errx = 0.0
    for i in range(len(x)):
        S1 += x[i]*np.power(ex[i],-2)
        S2 += np.power(ex[i],-2)
        errx += np.power(ex[i],-2)
    return [S1/S2, np.power(errx,-1)]

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
    for i in range(8, len(U)-8):
        maximum = True
        for j in range(i-8, i+8):
            if U[j]>U[i]:
                maximum = False
                break
        if maximum:
            max_wahr = True
            for k in range(len(m)):
                if np.round(m[k], 1)==np.round(U[i],1):
                    max_wahr = False
            if max_wahr:
                Tmax.append(t[i])
                m.append(U[i])
    return Tmax, m