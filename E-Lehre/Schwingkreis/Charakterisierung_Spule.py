# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 11:58:52 2018

@author: Christine Falter
"""

from __future__ import print_function
from praktikum import analyse
from praktikum import cassy
from scipy.signal import find_peaks_cwt
import numpy as np
from pylab import *

L_fin = []
C_fin = []
R_L_fin = []
d_fin = []

close('all')

'''
Dieses Programm soll den Schwingkreis analysieren 
sowie die verwendete Spule charakterisieren. Die Kabelwerte
könnten aus den gekoppelten Schwingungen bestimmt werden und 
müssten von R_Rest noch abgezogen werden um R_L zu erhalten

delta aus Einhüllende
freq aus Peaks
'''

def c_open(file):
    """
    Zum Einlesen der CASSY-Messungen
    Parameter: file - Messdatei
    returns: t, U, I 
    """

    data = cassy.CassyDaten(file)
    t = data.messung(1).datenreihe("t").werte
    I = data.messung(1).datenreihe("I_A2").werte
    U = data.messung(1).datenreihe("U_B2").werte 

    return t, U, I


def calc_mean_distance(ind,v):
    dis = 0
    i = 0
    while i < len(ind)-1:
        dis = dis + v[ind[i+1]]-v[ind[i]]
        i = i+1
        
    dis = dis/(len(ind)-1)
    return dis

def analyse_Schwingung(file,R):
    t, U, I = c_open(file)
    fig = figure()
    
    #Frequenz berechnen
    ind_max = find_peaks_cwt(U,np.arange(1,50))
    ind_max = ind_max[:12]
    ind_max = ind_max[1:]
    dis = calc_mean_distance(ind_max,t)
    
    hist(U[1000:])
    title("Offset-Korrektur")
    off_set = np.mean(U[1000:])
    print(off_set)
     
    U=U-off_set
    fig = figure()
    
    t_err = 2*(t[1]-t[0])
    freq = 1/dis
    freq_err = (freq/dis)*t_err
    
    U_std = np.ones(len(U))
    #wert der aus den Rauschmessung bei Berechnung des Widerstands als std. 
    #Abweichung angenommen werden kann
    U_std[:]=0.0027
    #Auflösungsfehler
    A_err = 10*2**(-11)/np.sqrt(12)
    print(A_err)
    if (A_err > U_std[0]):
        U_std[:] = A_err
    
    y = analyse.exp_einhuellende(t,U,U_std)
    
    subplot(2,1,1)
    title("Schwingkreis mit Widerstand "+str(R)+"$\Omega$")
    scatter(t,U,marker = ".")
    scatter(t[ind_max],U[ind_max],marker = "x",color = 'red')
    plot(t,y[0]*np.exp(-y[2]*t),color = "red")
    print("delta: ",y[2],"+-",y[3])
    print("frequenz: ",freq,"+-",freq_err)
    
    ylabel('$U$ / V')
    xlabel('t/s')
    grid()
    
    subplot(2,1,2)
    scatter(t,I,marker = ".")
    ylabel('$I$ / A')
    xlabel('t/s')
    grid()
    
    subplots_adjust(hspace = 0.5)
    #Gauss-Fehler-Fortpflanzung
    C = 2.2e-6 #dieser Wert muss evtl. angepasst werden/Charakterisierung Kondensator
    w = (freq*2*np.pi)
    w_err = (freq_err*2*np.pi)
    d = y[2]
    d_err = y[3]

    L = 1./((w**2+y[2]**2)*C)
    L_err = (2/((w**2+d**2)**2*C))**2*(w**2*w_err**2+d**2*d_err**2)
    L_err = np.sqrt(L_err)
    
    R_L = 2*d/((w**2+d**2)*C)-R
    R_L_freq_err = (4*w*d)/((w**2+d**2)**2*C)*w_err
    R_L_d_err = (2/((w**2+d**2)*C)+4*d**2/((w**2+d**2)**2*C))*d_err
    R_L_err = np.sqrt(R_L_freq_err**2+R_L_d_err**2)
    
    print("Induktivität: ",L,"+-",L_err,"H")
    print("R_Rest: ",R_L,"+-",R_L_err,"Ohm")
    print("Herstellerwerte: L = 0.009 R_L = 2.5")
    
    #Residuengraf
    fig = figure()
    errorbar(t[ind_max],y[0]*np.exp(-y[2]*t[ind_max])-U[ind_max],yerr = 0.0043,fmt='.')
    axhline(0, color="RED", linestyle='dashed', linewidth=1)

    d_fin.append(y[2])
    d_fin.append(y[3])
    L_fin.append(L)
    L_fin.append(L_err)
    R_L_fin.append(R_L)
    R_L_fin.append(R_L_err)
    C_fin.append(2.2e-6)
    C_fin.append(0.01e-6)

def analyse_Grenzfall():
    R_ap = 2* np.sqrt(L_fin[0]/C_fin[0])
    R_ap_err = R_ap *np.sqrt((L_fin[1]/(2*L_fin[0]))**2+(C_fin[1]/(2*C_fin[0]))**2)
    R_R_ap_err = np.sqrt(R_ap_err**2 + R_L_fin[1]**2)
    
    print("Aperiodischer Grenzfall bei", R_ap-R_L_fin[0], "+-",R_R_ap_err ," Ohm")

analyse_Schwingung("5.1_Messung_1.lab", 5.1)
analyse_Grenzfall()
