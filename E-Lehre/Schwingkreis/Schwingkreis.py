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

def s_U_sys(U): return 0.01+0.005*U

def delta(t, U):
    """
    Berechnet das Delta einer Schwingung
    """
    d, ed = [], []
    tmaxima, maxima = lok_max(t, U)
    if len(tmaxima) <= 3: return [1.0, 1.0]
    for i in range(len(maxima)-1):
        d.append(np.log(maxima[i]/maxima[i+1])/(tmaxima[i+1]-tmaxima[i]))
        ed.append(s_U_Ger * np.sqrt(np.power(maxima[i],-2)+np.power(maxima[i+1],-2))/((tmaxima[i+1]-tmaxima[i])))
    return gew_mittel(d, ed)

#Plot der Messungen mit dem Oszilloskop
#for i in range(4):
#    t1,t2, U1, U2 = [], [], [], []
#    t1, U1 = f_open("ALL000"+str(i)+"/U1.dat")
#    t2, U2 = f_open("ALL000"+str(i)+"/U2.dat")
#    plt.title("Entladung des Kondensators über den Schwingkreis \n R = ??? $\Omega$")
#    plt.plot(t1, U1, label="U1")
#    plt.plot(t2, U2, label="U2")
#    plt.savefig("Images/O_"+str(i)+".pdf")
#    plt.legend()
#    plt.show()
#    plt.close()

#Plot der Messungen mit den Widerständen R:
R = [5.1, 47, 100, 1000, 208, 173, 137, 100.6]
offsetkorr = [1601, 401, 171, 900, 343, 220, 150, 60]
R_plt = []
Delta_werte = []
Delta_werte_err = []
s_U_Ger = 10.0*np.power(2.0, -11) #Auflösungsvermögen in V

for i in range(len(R)):
    f = str(R[i])+"_Messung_"+str(i+1)+".lab"
    t, U, I = c_open(f)
    if R[i] == 173: #Abschneiden der Werte der Wiederaufladung
        U=U[:1501]
        t=t[:1501]
        I=I[:1501]
    #Korreketur des Offsets:
    dU = np.average(U[offsetkorr[i]:])
    U=U-dU

    #Schwingfall
    if i==0 or i==1:
        #Berechnung von Delta aus dem Amplitudenverhältnis
        d, ed = delta(t, np.abs(U))
        maxima = lok_max(t, np.abs(U))
        emax = []
        emax = s_U_Ger * np.ones(len(maxima[0]))
        if len(maxima[0])>3: plt.plot(maxima[0], maxima[1], marker="x")
        if d!=1.0: plt.plot(t, U[0]*np.exp(-d*np.array(t)), label = "Einhüllende")
        print("Die Dämpfungskonstante Delta berechnet aus der Amplitude: \n d=" + str(np.round(d,4)) + "+-" + str(np.round(ed,4)))
        
        #Berechnung von Delta aus der linearen Regression
        #TODO: Fehler zu groß!
        d2,ed2,b,eb,chiq,corr = anal.lineare_regression(np.array(maxima[0]), np.log(maxima[1]), np.log(emax))
        print("Die Dämpfungskonstante Delta berechnet aus der linearen Regression: \n d=" + str(-d2) + "+-" + str(ed2))
        #plot
        plt.title("Entladung des Kondensators über den Schwingkreis \n R = " +str(R[i]) + " $\Omega$")
        plt.errorbar(maxima[0],maxima[1],yerr=emax, fmt="x")
        plt.axvline(x=t[offsetkorr[i]], color="red", linestyle = "--")
        plt.plot(t, I*20, label="20$\cdot$I")
        plt.plot(t, U, label="U")
        plt.savefig("Images/C_"+str(R[i])+".pdf")
        plt.legend()
        plt.show()
        plt.close()
        
        plt.title("Entladung des Kondensators über den Schwingkreis - log \n R = " +str(R[i]) + " $\Omega$")
        plt.errorbar(maxima[0],np.log(maxima[1]),yerr=np.log(emax), fmt="x",markersize=4, capsize=3)
        plt.plot(t[:offsetkorr[i]], -d*t[:offsetkorr[i]]+b, label="Regressionsgerade")
        #plt.savefig("Images/C_"+str(R[i])+".pdf")
        plt.legend()
        plt.show()
        plt.close()
        
    #i=2: Grenzfall
    #Kriechfall
    if i>=3 and i!=7:
        #Lineare Regression für Delta
        t_sliced = t[:offsetkorr[i]]
        U_sliced = U[:offsetkorr[i]]
        d,ed,b,eb,chiq,corr = anal.lineare_regression(np.array(t_sliced), np.log(U_sliced), np.log(s_U_Ger*np.ones(len(U_sliced))))
        print("Die Dämpfungskonstante Delta berechnet aus der linearen Regression: \n d=" + str(d) + "+-" + str(ed))
        #plot
        plt.title("Entladung des Kondensators über den Schwingkreis \n R = " +str(R[i]) + " $\Omega$")
        plt.axvline(x=t[offsetkorr[i]], color="red", linestyle = "--")
        plt.plot(t, I*20, label="20$\cdot$I")
        plt.plot(t, U, label="U")
        plt.plot(t_sliced, np.exp(d*t_sliced+b), label="Einhüllende")
        plt.savefig("Images/C_"+str(R[i])+".pdf")
        plt.legend()
        plt.show()
        plt.close()
        
        plt.title("Entladung des Kondensators über den Schwingkreis - log \n R = " +str(R[i]) + " $\Omega$")
        plt.errorbar(t_sliced,np.log(U_sliced),yerr=np.log(s_U_Ger)*np.ones(len(U_sliced)), fmt="x",markersize=4, capsize=3)
        plt.plot(t[:int(offsetkorr[i]*1.1)], d*t[:int(offsetkorr[i]*1.1)]+b, label="Regressionsgerade")
        #plt.savefig("Images/C_"+str(R[i])+".pdf")
        plt.legend()
        plt.show()
        plt.close()
    
    R_plt.append(R[i])
    Delta_werte.append(np.abs(d))
    Delta_werte_err.append(np.abs(ed))
    
i2L,ei2L,R_Rest,eR_Rest,chiq,corr = anal.lineare_regression(np.array(R), np.array(Delta_werte), np.array(Delta_werte_err))
plt.errorbar(R, Delta_werte, Delta_werte_err, fmt="o",markersize=4, capsize=3)
plt.plot(R, np.array(R)*i2L+R_Rest)
plt.title("Die Veränderung der Dämpfungskonstante als Funktion von R")
plt.show()
plt.close()
print((i2L*2)**-1, " WTF??")

print("===="*10)
print("Also was wir haben:")
R_real=5.12
print(2*9*10**-3*Delta_werte[0]-R_real)
    
##Grenzfall:
#f = "Grenzfall_Messung_9.lab"
#t, U, I = c_open(f)
#plt.title("Entladung des Kondensators über den Schwingkreis \n Grenzfall bei \n R = 111,3 $\Omega$")
#plt.plot(t, I, label="I")
#plt.plot(t, U, label="U")
#plt.savefig("Images/grenzfall.pdf")
#plt.legend()
#plt.show()
#plt.close()