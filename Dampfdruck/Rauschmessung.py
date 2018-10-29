#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 27 17:24:44 2018

@author: Mate
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy
import praktikum.cassy as cassy

def gauss(x, m, s):
    return(np.power(scipy.pi*2*s**2, -0.5)*np.exp(-(x-m)**2/(2*s**2)))

t, p, T = [], [], [] #Temperatur in °C!
data = cassy.CassyDaten("Rauschmessung.lab")
t = data.messung(1).datenreihe("t").werte
T = data.messung(1).datenreihe("&J_A11").werte
p = data.messung(1).datenreihe("p_B1").werte

print("Temperaturverteilung")
print("m=" + str(np.average(T))+" ; s="+str(np.std(T))+", err="+str(np.std(T)))
plt.hist(T,7, normed=True)
x=np.arange(np.average(T)-0.25, np.average(T)+0.25, 0.001)
plt.plot(x, gauss(x, np.average(T), np.std(T)))
plt.title(u"Rauschmessung bei Raumtemperatur\n $\sigma$=" + str(np.std(T))+" ; $\mu$="+str(np.average(T)))
plt.xlabel(u"T/°C")
plt.ylabel('Relatives Vorkommen')
plt.savefig("Images/RauschmessungRT_T_histo.jpg")
plt.figure()

plt.plot(t, T)
plt.title("Rauschmessung bei Raumtemperatur")
plt.ylabel(u"T/°C")
plt.xlabel('Zeit/s')
plt.plot(t,np.array([np.mean(T)]*len(t)),color='green')
plt.savefig("Images/RauschmessungRT_T.jpg")
plt.figure()

print("==="*20)
print("Druckverteilung") #vllt Poisson???
print("m=" + str(np.average(p))+", s="+str(np.std(p))+", err="+str(np.std(p)/np.sqrt(len(p))))

x=np.arange(min(p), max(p), 0.01)
plt.hist(p, 3, normed=True)#, bins=np.arange(984, 986, 0.25))
plt.plot(x, gauss(x, np.average(p), np.std(p)))
plt.title(u"Rauschmessung Druck \n s=" + str(np.std(p))+", m="+str(np.average(p)))
plt.xlabel(u"Druck/hPa")
plt.ylabel('Relatives Vorkommen')
plt.savefig("Images/RauschmessungRT_p_histo.jpg")
plt.figure()

plt.plot(t, p)
plt.ylabel(u"Druck/hPa")
plt.xlabel('Zeit/s')
plt.title("Rauschmessung Druck")
plt.savefig("Images/RauschmessungRT_p.jpg")
plt.plot(t,np.array([np.mean(p)]*len(t)),color='green')
