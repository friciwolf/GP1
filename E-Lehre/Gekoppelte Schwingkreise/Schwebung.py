# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 00:43:38 2018

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

data = cassy.CassyDaten('Schwebung/Schwebung_1_250n.lab')
U1 = data.messung(1).datenreihe('U_B1').werte
t = data.messung(1).datenreihe('t').werte

ind_max = find_peaks_cwt(U1,np.arange(1,20))
ind_max = ind_max[:20]
ind_max = ind_max[5:]


w_array,A = analyse.fourier_fft(t,U1)
ind_w = find_peaks_cwt(A,np.arange(1,20))
axvline(t[ind_w[0]], color="RED", linestyle='dashed', linewidth=1)


figure()
plot(w_array,A)
xlim(0,4000)

figure()
#Rohdaten
scatter(t,U1,marker = ".")
ylabel('$U$ / V')
xlabel('t/s')
axvline(t[len(t)-4000], color="RED", linestyle='dashed', linewidth=1)
grid()

#off_set anpassung
figure()
hist(U1[4000:])
title("Offset-Korrektur")
off_set = np.mean(U1[4000:])
print(off_set)

U1 = U1-off_set

t_err = 2*(t[1]-t[0])
#Peakanalyse
figure()
scatter(t[:3000],U1[:3000],marker = ".")
scatter(t[ind_max],U1[ind_max],marker = "x",color = 'red')
dis = calc_mean_distance(ind_max,t)


freq1 = 1/dis
freq1_err = (freq1/dis)*t_err
print("Frequenz 1 aus Peaks: ",freq1, "+-",freq1_err)
ylabel('$U$ / V')
xlabel('t/s')
grid()

