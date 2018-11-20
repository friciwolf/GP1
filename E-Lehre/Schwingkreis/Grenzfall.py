#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 17:33:47 2018

@author: Mate
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 08:35:53 2018

@author: Mate
"""
import matplotlib.pyplot as plt
import praktikum.analyse as anal
import numpy as np
from defs import *

#Plot der Messungen mit den Widerständen R:
R = [5.1, 47, 100, 1000, 208, 173, 137, 100.6]
offsetkorr = [1601, 401, 171, 900, 343, 220, 150, 60]

for i in range(2, len(R)):
    if R[i]!=100.6 and R[i]!=100:
        f = str(R[i])+"_Messung_"+str(i+1)+".lab"
        t, U, I = c_open(f)
        if R[i] == 173: #Abschneiden der Werte der Wiederaufladung
            U=U[:1501]
            t=t[:1501]
            I=I[:1501]
        #Korreketur des Offsets:
        dU = np.average(U[offsetkorr[i]:])
        U=U-dU
        #Lineare Regression für Delta
        t_sliced = t[:offsetkorr[i]]
        U_sliced = U[:offsetkorr[i]]
        #plot
        plt.title("Entladung des Kondensators über den Schwingkreis")
        plt.xlabel("t /s")
        plt.ylabel("$U_C$ /V")
        plt.plot(t[:min(offsetkorr)], U[:min(offsetkorr)], label=str(R[i])+ " Ω")

#Grenzfall:
f = "Grenzfall_Messung_9.lab"
t, U, I = c_open(f)
plt.plot(t[:min(offsetkorr)], U[:min(offsetkorr)],"--", label="Beobachteter Grenzfall - 111.3 Ω")
plt.legend()
plt.savefig("Images/Grenzfall_plots.png")
plt.show()
plt.close()