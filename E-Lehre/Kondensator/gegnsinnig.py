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
    dis = 0
    i = 0
    while i < len(ind)-1:
        dis = dis + v[ind[i+1]]-v[ind[i]]
        i = i+1
        
    dis = dis/(len(ind)-1)
    return dis

data = cassy.CassyDaten('Gegensinnig/gegen.koppl.eisen.lab')
U1 = data.messung(1).datenreihe('U_B1').werte
U2 = data.messung(1).datenreihe('U_B2').werte
t = data.messung(1).datenreihe('t').werte

w1,A1 = analyse.fourier_fft(t,U1)
w1 = analyse.peakfinder_schwerpunkt(w1,A1)
print("Frequenz 1 aus fft: ",w1)

w2,A2 = analyse.fourier_fft(t,U2)
w2 = analyse.peakfinder_schwerpunkt(w2,A2)
print("Frequenz 2 aus fft: ",w2)

w_mean = (w1+w2)/2
w_std = np.sqrt((w1-w_mean)**2+(w2-w_mean)**2)
print("mittlere Frequenz aus fft",w_mean,"+-",w_std)

fig = figure()

subplot(2,1,1)
plot(t,U1)
ylabel('$U$ / V')
xlabel('t/s')
title('erster Schwingkreis')
grid()

subplot(2,1,2)
plot(t,U2)
ylabel('$U$ / V')
xlabel('t/s')
title('zweiter Schwingkreis')
grid()

subplots_adjust(hspace = 0.5)

fig.savefig('Auswertung/Rohdaten_gegensinnig.pdf')

fig = figure()

#Fehler für die Peakbestimmung etwa 2 Messpunkte
t_err = 2*(t[1]-t[0])

subplot(2,1,1)
ind_max = find_peaks_cwt(U1,np.arange(1,20))
ind_max = ind_max[:16]

scatter(t[:3000],U1[:3000],marker = ".")
scatter(t[ind_max],U1[ind_max],marker = "x",color = 'red')
dis = calc_mean_distance(ind_max,t)
freq1 = 1/dis
freq1_err = (freq1/dis)*t_err
print("Frequenz 1 aus Peaks: ",freq1, "+-",freq1_err)
ylabel('$U$ / V')
xlabel('t/s')
title('erster Schwingkreis')
grid()

subplot(2,1,2)
ind_max = find_peaks_cwt(U2,np.arange(1,20))
ind_max = ind_max[:16]

scatter(t[:3000],U2[:3000],marker = ".")
scatter(t[ind_max],U2[ind_max],marker = "x",color = 'red')
dis = calc_mean_distance(ind_max,t)
freq2 = 1/dis
freq2_err = (freq2/dis)*t_err
print("Frequenz 2 aus Peaks: ",freq2, "+-",freq2_err)
ylabel('$U$ / V')
xlabel('t/s')
title('zweiter Schwingkreis')
grid()

subplots_adjust(hspace = 0.5)

#gewichteter Mittelwert/Fehler
freq_mean = (freq1/freq1_err**2+freq2/freq2_err**2)/(1/freq1_err**2+1/freq2_err**2)
freq_std = np.sqrt(1/(1/freq1_err**2+1/freq2_err**2))
print("mittlere Frequenz aus Peaks",freq_mean,"+-",freq_std)

fig.savefig('Auswertung/Peakanalyse_gegensinnig.pdf')

