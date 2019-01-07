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

R, U = c_open("V2.1C.lab")
plt.plot(R, U, ls="None", marker="x", markersize=2)
plt.show()