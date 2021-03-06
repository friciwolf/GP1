#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 27 17:24:44 2018

@author: Mate, Patrick
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy
import praktikum.cassy as cassy

def gauss(x, m, s):
    return(np.power(scipy.pi*2*s**2, -0.5)*np.exp(-(x-m)**2/(2*s**2)))

def round_up(x, digits): #für obere Abschätzung der Fehler
    if digits!=0: return np.ceil(x*10**digits)/10**digits
    else: return np.ceil(x)
    
def round_good(m,s,err):
    '''
    Passt Mittelwert- und Standardabweichungsnachkommastellen an den Fehler an (2 signifikante Stellen des Fehlers)
    m=mean
    s=std
    err=std/sqrt(len)
    '''
    i=2
    test=list('{:.20f}'.format(err))
    while test[i]=='0': #TODO: wie groß sind die Fehler auf L? Kommastelle finden oder nicht?
        i+=1
    return np.round([m,s,err],i)

plt.close('all')

t, p, T = [], [], [] #Temperatur in °C, Druck in hPa!
data = cassy.CassyDaten("Rauschmessung.lab")
t = data.messung(1).datenreihe("t").werte
T = data.messung(1).datenreihe("&J_A11").werte
p = data.messung(1).datenreihe("p_B1").werte

print("Temperaturverteilung")
print("m=" + str(np.round(np.average(T),1))+" °C , s="+str(np.round(np.std(T),2))+" °C , err="+str(round_up(np.std(T)/np.sqrt(len(T)),1))+" °C")
plt.hist(T,7, normed=True)
x=np.arange(np.average(T)-0.25, np.average(T)+0.25, 0.001)
plt.plot(x, gauss(x, np.average(T), np.std(T)))
plt.title(u"Rauschmessung bei Raumtemperatur\n $\sigma$=" + str(np.round(np.std(T),2))+"$^{\circ}C$, $\mu$="+str(np.round(np.average(T),1))+ "$^{\circ}C$, $\sigma_{\mu} = $"+str(round_up(np.std(T)/np.sqrt(len(T)),1))+"$^{\circ}C$")
plt.xlabel(u"T/°C")
plt.ylabel('Relatives Vorkommen')
plt.savefig("Images/RauschmessungRT_T_histo.jpg")
plt.figure()

plt.plot(t, T)
plt.title(u"Rauschmessung bei Raumtemperatur\n $\sigma$=" + str(np.round(np.std(T),1))+"$^{\circ}C$, $\mu$="+str(np.round(np.average(T),1))+ "$^{\circ}C$, $\sigma_{\mu} = $"+str(round_up(np.std(T)/np.sqrt(len(T)),1))+"$^{\circ}C$")
plt.ylabel(u"T/°C")
plt.xlabel('Zeit/s')
plt.plot(t,np.array([np.mean(T)]*len(t)),color='green')
plt.savefig("Images/RauschmessungRT_T.jpg")
plt.figure()

print("==="*20)
print("Druckverteilung")
print("m=" + str(np.round(np.average(p),2))+" hPa, s="+str(np.round(np.std(p),2))+" hPa, err="+str(round_up(np.std(p)/np.sqrt(len(p)),2))+ " hPa")

x=np.arange(min(p), max(p), 0.01)
plt.hist(p, 3, normed=True)#, bins=np.arange(984, 986, 0.25))
plt.plot(x, gauss(x, np.average(p), np.std(p)))
plt.title(u"Rauschmessung Druck \n $\sigma$=" + str(np.round(np.std(p),2))+" hPa, $\mu$="+str(np.round(np.average(p),2))+" hPa, $\sigma_{\mu} = $" +str(round_up(np.std(p)/np.sqrt(len(p)),2))+" hPa")
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
plt.figure()

#----------------------------#
#Auswertung der Rauschmessung#
#----------------------------#

#Mittelwerte, Standardabweichungen und Fehler
mT,sT,errT=round_good(np.average(T),np.std(T),np.std(T)/np.sqrt(len(T)))
mp,sp,errp=round_good(np.average(p),np.std(p),np.std(p)/np.sqrt(len(p)))

if __name__ == "__main__": 
    
#------------------------------------------------------------
#Umgebungstemperatur
    print("Temperaturverteilung")
    print('m='+str(mT)+u' °C, std='+str(sT)+u' °C, err='+str(errT)+u' °C')
    
    #Histogramm
    plt.hist(T,7, normed=True)
    x=np.arange(np.average(T)-0.25, np.average(T)+0.25, 0.001)
    plt.plot(x, gauss(x, np.average(T), np.std(T)))
    plt.title(u"Histogramm der Raumtemperatur\n $\sigma$=" + str(sT)+"$^{\circ}C$, $\mu$="+str(mT)+ "$^{\circ}C$, $\sigma_{\mu} = $"+str(errT)+"$^{\circ}C$")
    plt.xlabel(u"T/°C")
    plt.ylabel('Relatives Vorkommen')
    plt.savefig("Images/RauschmessungRT_T_histo.jpg")
    plt.figure()
    
    #Datenplot
    plt.plot(t, T)
    plt.title(u"Rauschmessung bei Raumtemperatur")
    plt.ylabel(u"T/°C")
    plt.xlabel('Zeit/s')
    plt.plot(t,np.array([np.mean(T)]*len(t)),color='green')
    plt.savefig("Images/RauschmessungRT_T.jpg")
    plt.figure()

#------------------------------------------------------------
#Umgebungsdruck
    print("==="*20)
    print("Druckverteilung")
    print('m='+str(mp)+' hPa, std='+str(sp)+' hPa, err='+str(errp)+' hPa')
    
    #Histogramm
    x=np.arange(min(p), max(p), 0.01)
    plt.hist(p, 3, normed=True)#, bins=np.arange(984, 986, 0.25))
    plt.plot(x, gauss(x, np.average(p), np.std(p)))
    plt.title(u"Histogramm Umgebungsdruck \n $\sigma$=" + str(sp)+" hPa, $\mu$="+str(mp)+" hPa, $\sigma_{\mu} = $" +str(errp)+" hPa")
    plt.xlabel(u"Druck/hPa")
    plt.ylabel('Relatives Vorkommen')
    plt.savefig("Images/RauschmessungRT_p_histo.jpg")
    plt.figure()
    
    #Datenplot
    plt.plot(t, p)
    plt.ylabel(u"Druck/hPa")
    plt.xlabel('Zeit/s')
    plt.title("Rauschmessung Druck")
    plt.plot(t,np.array([np.mean(p)]*len(t)),color='green')
    plt.savefig("Images/RauschmessungRT_p.jpg")
