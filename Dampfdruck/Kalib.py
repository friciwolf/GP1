#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 27 19:27:20 2018

@author: Mate
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy
import praktikum.cassy as cassy
import praktikum.analyse as ana

def gauss(x, m, s):
    return(np.power(scipy.pi*2*s**2, -0.5)*np.exp(-(x-m)**2/(2*s**2)))

t, p, T = [], [], [] #Temperatur in °C!!!
data = cassy.CassyDaten("Kalibrierung 1 + Dichtigkeit.lab")
t = data.messung(1).datenreihe("t").werte
T = data.messung(1).datenreihe("&J_A11").werte
p = data.messung(1).datenreihe("p_B1").werte
#Fehler ausrechnen!!
sigmat=np.std(t)
sigmaT=np.std(T)
sigmap=np.std(p)
errt=sigmat/np.sqrt(len(t))
errT=sigmaT/np.sqrt(len(T))
errp=sigmap/np.sqrt(len(p))
####Temperaturkalibrierung:

#T=0°C
N_Tslice=2000
t_sliced=t[N_Tslice:]
T=T[N_Tslice:]
print("Temperaturkalibrierung")
T0=np.average(T)
sigmaT0 = np.std(T)
print("s=" + str(sigmaT0)+", m="+str(T0))
plt.title("Temperatur bei T=0°C\n s=" + str(sigmaT0)+", m="+str(T0))
plt.plot(np.arange(min(t_sliced), max(t_sliced), 0.1),np.ones(len(t_sliced))*T0, color="green")
plt.plot(t_sliced, T)
plt.xlabel('time/s')
plt.ylabel(u'T/°C')
plt.figure()

x=np.arange(min(T), max(T), 0.01)
plt.hist(T, normed=True)
plt.plot(x, gauss(x, np.average(T), np.std(T)))
plt.title("Streuung der gemessenen Temperaturwerte bei T=0°C")
plt.ylabel('rel')
plt.xlabel(u'T/°C')
plt.figure() #=> mehr oder weniger Gauß -> man kann den Durchschnitt nehmen...

plt.plot(t, p)
m,em,b,eb,chiq,corr = ana.lineare_regression_xy(t, p, sigmat*np.ones(len(t)), sigmap*np.ones(len(p)))
plt.plot(t, m*t+b, color='green')
plt.title("Dichtigkeitsmessung aka Druckanstieg im evakuierten Kolben")
plt.xlabel('time/s')
plt.ylabel(u'p/hPa')
plt.figure()
print('m = (%g +- %g) hPa/s,   b = (%g +- %g) hPa,  chi2/dof = %g / %g  corr = %g' % (m, em, b, eb, chiq, len(t)-2, corr))
print('Verlust: ({} +-{}) mbar/min'.format(60*m,60*em)) #so genehmigt

#T=100°C:
data = cassy.CassyDaten("Temperatur, Druck bei 100.lab")
t = data.messung(1).datenreihe("t").werte
p = data.messung(1).datenreihe("p_B1").werte
T = data.messung(1).datenreihe("&J_A11").werte

N_slice = [240, 1500] #Betrachte dieses Werteintervall, weil außerhalb von diesem komische Sachen passieren...
t = t[240:1500]
T = T[240:1500]
p = p[240:1500]
T100 = np.average(T)
sigmaT100 = np.std(T)

x = np.arange(min(t),max(t), 0.01)
plt.plot(t, p)
plt.plot(x, np.ones(len(x))*np.average(p), c="green")
plt.title("Dampfdruck bei T=100°C\n m = " + str(np.average(p)))
plt.xlabel('time/s')
plt.ylabel(u'p/hPa')
plt.figure()

plt.plot(t, T)
x = np.arange(min(t),max(t), 0.01)
plt.plot(x, np.ones(len(x))*T100, color="green")
plt.title("Temperatur bei T=100°C\n m ="+ str(T100))
plt.xlabel('time/s')
plt.ylabel(u'T/°C')
plt.figure()  

x = np.arange(min(T),max(T), 0.01)
plt.hist(T, normed=True)
plt.plot(x, gauss(x, T100, sigmaT100))
plt.title("Streuung der gemessenen Temperaturwerte bei T=100°C")
plt.ylabel('rel')
plt.xlabel(u'T/°C')
#=> merkwürdig... sehr merkwürdig, die hier ist nicht gaussförmig genug...

print("==="*20)
print("Resultate:")
print("T bei T=0°C = " + str(T0)+" +-"+str(sigmaT0)+ " oder +-" + str(0.04443156295496722) + " (aus der Rauschmessung)"+ " oder +-" +" DIGITALISIERUNGSFEHLER")
print("T bei T=100°C = " + str(T100)+" +-"+str(sigmaT100)+ " oder +-" + str(0.04443156295496722) + " (aus der Rauschmessung)"+ " oder +-" +" DIGITALISIERUNGSFEHLER")
#TODO: Fehler einfügen
a,ea,b,eb,chiq,corr = ana.lineare_regression_xy(np.array([0, 100])-50, np.array([T0, T100])-50, np.ones(2), np.ones(2))
print("Damit ergibt sich die Korrekturgleichung:")
print("T_real = "+ str(a)+"*(T_gemessen-50) + " + str(b)+" + 50")
