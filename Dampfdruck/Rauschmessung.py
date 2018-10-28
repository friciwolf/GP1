#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 27 17:24:44 2018

@author: Mate
Hi!
"""
####
import matplotlib.pyplot as plt
import numpy as np
import scipy

def gauss(x, m, s):
    return(np.power(scipy.pi*2*s**2, -0.5)*np.exp(-(x-m)**2/(2*s**2)))

f = open("pT-Gauss.dat")
t, p, T = [], [], [] #Temperatur in °C!!!
for line in f:
    data = line.split("\t")
    t.append(float(data[1]))
    T.append(float(data[2]))
    p.append(float(data[4]))
    
print("Temperaturverteilung") #vllt Poisson???
print("s=" + str(np.std(T))+", m="+str(np.average(T)))
plt.hist(T)
x=np.arange(min(T), max(T), 0.001)
plt.plot(x, gauss(x, np.average(T), np.std(T))*100) #TODO: Falsche Verteilung??? Faktor 100 daneben!!
plt.hist(T, normed=True, bins=np.arange(min(T), max(T), 0.5))
plt.title("Temperaturverteilung - falsch, Gauss mit Faktor 100!\n s=" + str(np.std(T))+", m="+str(np.average(T)))
plt.show()

plt.plot(t, T)
plt.title("Rauschmessung Raumtemperatur")
plt.show()

print("==="*20)
print("Druckverteilung") #vllt Poisson???
print("s=" + str(np.std(p))+", m="+str(np.average(p)))

x=np.arange(min(p), max(p), 0.01)
plt.plot(x, gauss(x, np.average(p), np.std(p)))
plt.hist(p, normed=True, bins=np.arange(984, 986, 0.25))
plt.title("Druckverteilung \n s=" + str(np.std(p))+", m="+str(np.average(p)))
plt.show()

plt.plot(t, p)
plt.title("Rauschmessung Atmosphärendruck")
plt.show()
