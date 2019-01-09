# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 19:57:55 2018

@author: Anwender
"""

import Temperatur
import matplotlib.pyplot as plt
import numpy as np
import praktikum.cassy as cassy
import praktikum.analyse as analyse

#Daten einlesen
data = cassy.CassyDaten("Kalibrierung_Wegaufnehmer.lab")
n = np.array(data.messung(1).datenreihe("n").werte) #Nummer des Messpunktes, wobei n=1 der weiteste Abstand war
R = np.array(data.messung(1).datenreihe("R_B1").werte) #Widerstand in kOhm

#n->L
L=(5*len(n)-5*(n)) #in cm #es wurde in 5cm-Schritten gemessen #TODO: np.arange(0.5,30.5,5)
#So ist L=0 bei der Messung, die am nächsten an der Quelle war und L wird bei größerem Abstand von der Quelle größer

L_err=np.ones(len(L))*0.1/2 #cm
R_err=np.ones(len(R))*0.005/np.sqrt(12) #kOhm
L_esys=0.03+0.02*1 #cm #Güteklasse II:a+b*L mit a=0.3mm, b=0.2mm/m, L= auf ganzen Meter gerundete zu messende Länge
                    #=const.
#lineare Regression
k,k_err,L0,L0_err,chiq,corr=analyse.lineare_regression_xy(R,L,R_err,L_err)
chiq,k,k_err=Temperatur.round_good(chiq,k,k_err)
x,L0,L0_err=Temperatur.round_good(0,L0,L0_err)

if __name__=='__main__':
    print('\nRegression L=R*x+L0=({} +- {})cm/kOhm * R+({} +- {})cm'.format(k,k_err,L0,L0_err))
    
    #Plot: Messpaare und Regression
    plt.title('Regression: Kalibrierung')
    plt.errorbar(R,L,xerr=R_err,yerr=L_err,fmt='ko', label=u'Messpaare                       Chi²/f={}/{}'.format(chiq,len(R)-2))
    x=np.arange(min(R),max(R),0.001)
    plt.plot(x,k*x+L0,label=r'$k*R+L_0=({}\pm{})cm/k\Omega * R+({}\pm{})cm$'.format(k,k_err,L0,L0_err))
    plt.xlabel(r'Widerstand R/k$\Omega$')
    plt.ylabel('Länge L/cm')
    plt.legend()
    plt.savefig('Images/Kalibrierung_Regression.pdf')
    plt.figure() 
    
    #Residuum
    plt.title('Residuenplot')
    plt.plot(x,np.zeros(len(x)),'r--')
    plt.errorbar(R,L-(k*R+L0),yerr=np.sqrt(L_err**2+k**2*R_err**2),fmt='ko',capsize=3)
    #plt.text(1.2,0.08,'Chi²/f={}/{}'.format(chiq,len(R)-2))
    plt.xlabel(r'R/$\Omega$')
    plt.ylabel(r'$L-(k*R+L_0)$ /cm')
    plt.savefig('Images/Kalibrierung_Residuum.pdf')

def Länge(R,Rerr):
    '''
    Rechnet die Länge aus dem Widerstand des Wegaufnehmers aus
    R,Rerr: Widerstand-Array, statistisches Fehler-Array
    returns: L,L_estat,L_esys
    '''
    return (k*R+L0,k*Rerr,np.sqrt(k_err**2*R**2)) #L,L_estat,L_esys    #+(L0_err**2+L_esys**2)))  #TODO: Systematischer mit L0? Was ist bei Differenzen?
