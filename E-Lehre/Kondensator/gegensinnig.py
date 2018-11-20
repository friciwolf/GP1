# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 14:19:52 2018

@author: Christine Falter
"""

from __future__ import print_function
from praktikum import analyse
from praktikum import cassy
from scipy.signal import find_peaks_cwt
import numpy as np
from pylab import *


close('all')

def calc_mean_distance(ind,v):
    '''
    Berechnet den Mittelwert der Abstände einer Peakgruppe für die mittlere Periodendauer
    ind:    array der Peakindexe
    v:      array der Zeit
    return: Mittelwert T
    '''
    dis = 0
    i = 0
    while i < len(ind)-1:
        dis = dis + v[ind[i+1]]-v[ind[i]]
        i = i+1
        
    dis = dis/(len(ind)-1)
    return dis

def round_good(m,s,err): 
    '''
    Passt Mittelwert- und Standardabweichungsnachkommastellen an den Fehler an (2 signifikante Stellen des Fehlers)
    m=mean
    s=std
    err=std/sqrt(len)
    '''
    i=2
    test=list('{:.20f}'.format(err))
    while test[i]=='0':
        i+=1
    return np.round([m,s,err],i)

#Daten einlesen
data = cassy.CassyDaten('gegen.koppl.eisen.lab')
U1 = data.messung(1).datenreihe('U_B1').werte
U2 = data.messung(1).datenreihe('U_B2').werte
t = data.messung(1).datenreihe('t').werte

#Rohdatenplot
subplot(2,1,1)
plot(t,U1)
ylabel('$U$ / V')
xlabel('t/s')
xlim(0,0.04)
title('erster Schwingkreis')
grid()

subplot(2,1,2)
plot(t,U2)
ylabel('$U$ / V')
xlabel('t/s')
xlim(0,0.04)
title('zweiter Schwingkreis')
grid()

subplots_adjust(hspace = 0.5)
#fig.savefig('Images/Rohdaten_gegensinnig.pdf')
plt.figure()

#FFT + Plot
w1_array,A1 = analyse.fourier_fft(t,U1)
w1 = analyse.peakfinder_schwerpunkt(w1_array,A1)#Hauptfrequenz1
print("Frequenz 1 aus fft: ",w1)

plt.plot(w1_array, A1/max(A1))
plt.axvline(x=w1, color="darkred", linestyle = "--") 
plt.text(w1,0.5,'Max: {} Hz'.format(round(w1,3)))
plt.title('FFT Schwingung 1')
plt.xlim(0,4000)
plt.xlabel('Frequenz / Hz')
plt.ylabel('rel. Häufigkeiten')
plt.figure()


w2_array,A2 = analyse.fourier_fft(t,U2)
w2 = analyse.peakfinder_schwerpunkt(w2_array,A2) #Hauptfrequenz2
print("Frequenz 2 aus fft: ",w2)


plt.plot(w2_array, A2/max(A2))
plt.axvline(x=w2, color="darkred", linestyle = "--") 
plt.text(w1,0.5,'Max: {} Hz'.format(round(w2,3)))
plt.title('FFT Schwingung 2')
plt.xlim(0,4000)
plt.xlabel('Frequenz / Hz')
plt.ylabel('rel. Häufigkeiten')

#Hauptfrequenz
w_mean = (w1+w2)/2
w_std = np.sqrt((w1-w_mean)**2+(w2-w_mean)**2)
x,w_mean,w_std=round_good(0.0,w_mean,w_std)
print("mittlere Frequenz aus fft",w_mean,"+-",w_std)

fig = figure()

#Peakbestimmung + Plots
t_err = 2*(t[1]-t[0]) #Fehler für die Peakbestimmung etwa 2 Messpunkte

ind_max = find_peaks_cwt(U1,np.arange(1,20))
ind_max = ind_max[:16]

subplot(2,1,1)
scatter(t[:3000],U1[:3000],marker = ".")
scatter(t[ind_max],U1[ind_max],marker = "x",color = 'red')
ylabel('$U$ / V')
xlabel('t/s')
title('erster Schwingkreis')
grid()

dis = calc_mean_distance(ind_max,t)
freq1 = 1/dis
freq1_err = (freq1/dis)*t_err
print("\nFrequenz 1 aus Peaks: ",freq1, "+-",freq1_err)

ind_max = find_peaks_cwt(U2,np.arange(1,20))
ind_max = ind_max[:16]

subplot(2,1,2)
scatter(t[:3000],U2[:3000],marker = ".")
scatter(t[ind_max],U2[ind_max],marker = "x",color = 'red')
ylabel('$U$ / V')
xlabel('t/s')
title('zweiter Schwingkreis')
grid()

dis = calc_mean_distance(ind_max,t)
freq2 = 1/dis
freq2_err = (freq2/dis)*t_err
print("Frequenz 2 aus Peaks: ",freq2, "+-",freq2_err)

subplots_adjust(hspace = 0.5)
#fig.savefig('Images/Peakanalyse_gegensinnig.pdf')

#gewichteter Mittelwert/Fehler
x,freq_mean,freq_std=round_good(0.0,*analyse.gewichtetes_mittel(np.array([freq1,freq2]),np.array([freq1_err,freq2_err])))
print("mittlere Frequenz aus Peaks",freq_mean,"+-",freq_std)

