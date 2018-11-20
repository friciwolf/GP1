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

close('all')

'''
Dieses Programm soll den Schwingkreis analysieren 
sowie die verwendete Spule charakterisieren. Die Kabelwerte
könnten aus den gekoppelten Schwingungen bestimmt werden und 
müssten von R_Rest noch abgezogen werden um R_L zu erhalten
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
    ind_max = ind_max[:16]
    dis = calc_mean_distance(ind_max,t)
    
    t_err = 2*(t[1]-t[0])
    freq = 1/dis
    freq_err = (freq/dis)*t_err
    
    
    y = analyse.exp_einhuellende(t,U,0.01*U)
    
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

analyse_Schwingung("5.1_Messung_1.lab",5.1)