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

t, p, T = [], [], [] #Temperatur in °C!!!
data = cassy.CassyDaten("Rauschmessung.lab")
t = data.messung(1).datenreihe("t").werte
T = data.messung(1).datenreihe("&J_A11").werte
p = data.messung(1).datenreihe("p_B1").werte

print("Temperaturverteilung") #vllt Poisson???
print("m=" + str(np.average(T))+", s="+str(np.std(T))+", err="+str(np.std(T)/np.sqrt(len(T))))
plt.hist(T,7, normed=True)
x=np.arange(min(T), max(T), 0.001)
plt.plot(x, gauss(x, np.average(T), np.std(T)))
plt.title(u"Temperaturverteilung+Gauß\n s=" + str(np.std(T))+", m="+str(np.average(T)))
plt.xlabel(u"T/°C")
plt.ylabel('rel')
plt.figure()

plt.plot(t, T)
plt.title("Rauschmessung Raumtemperatur")
plt.ylabel(u"T/°C")
plt.xlabel('time/s')
plt.plot(t,np.array([np.mean(T)]*len(t)),color='green')
plt.figure()

print("==="*20)
print("Druckverteilung") #vllt Poisson???
print("m=" + str(np.average(p))+", s="+str(np.std(p))+", err="+str(np.std(p)/np.sqrt(len(p))))

x=np.arange(min(p), max(p), 0.01)
plt.plot(x, gauss(x, np.average(p), np.std(p)))
plt.hist(p, 3, normed=True)#, bins=np.arange(984, 986, 0.25))
plt.title(u"Druckverteilung+Gauß \n s=" + str(np.std(p))+", m="+str(np.average(p)))
plt.xlabel(u"p/hPa")
plt.ylabel('rel')
plt.figure()

plt.plot(t, p)
plt.ylabel(u"p/hPa")
plt.xlabel('time/s')
plt.title("Rauschmessung Atmosphärendruck")
plt.plot(t,np.array([np.mean(p)]*len(t)),color='green')
