# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 23:18:32 2018

@author: Christine Falter
"""

import matplotlib.pyplot as plt

import numpy as np

import scipy

import praktikum.cassy as cassy

import praktikum.analyse as analyse



def gauss(x, m, s):

    return 1/np.sqrt(np.pi*2)/s*np.exp(-(x-m)**2/(2*s**2))



def round_good(m,s,err): #für obere Abschätzung der Fehler

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



plt.close('all')





data = cassy.CassyDaten("gekoppelt/Rauschmessung.lab")

t = data.messung(1).datenreihe("t").werte

M = data.messung(1).datenreihe("U_A1").werte #Pendelkörper (Masse)

St = data.messung(1).datenreihe("U_B1").werte #Stange allein



#----------------------------#

#Auswertung der Rauschmessung#

#----------------------------#



#Mittelwerte, Standardabweichungen und Fehler

mM,sM,errM=round_good(np.average(M),np.std(M),np.std(M)/np.sqrt(len(M)))

mSt,sSt,errSt=round_good(np.average(St),np.std(St),np.std(St)/np.sqrt(len(St)))



if __name__ == "__main__": 

    

#------------------------------------------------------------


    print("Spannungsverteilung 1")

    print('m='+str(mM)+u' V, std='+str(sM)+u' V, err='+str(errM)+u' V')

    

    #Histogramm

    weightsM=np.ones_like(M)/float(len(M))

    plt.hist(M,7, normed=0, weights=weightsM)

    x=np.arange(np.average(M)-0.03, np.average(M)+0.03, 0.001)

    #plt.plot(x, gauss(x, np.average(M), np.std(M)))

    plt.title(u"Histogramm der Ruhelage der Pendelstange 1 \n $\sigma$=" + str(sM)+"$V$, $\mu$="+str(mM)+ "$V$, $\sigma_{\mu} = $"+str(errM)+"$V$")

    plt.xlabel(u"U/V")

    plt.ylabel('Relatives Vorkommen')

    plt.savefig("Auswertung/Rauschmessung_stange_1_histo.jpg")

    plt.axvline(np.mean(M),color='red',linestyle='--')

    plt.figure()

    

    #Datenplot

    plt.plot(t, M)

    plt.title(u"Rauschmessung Ruhelage Pendel 1")

    plt.ylabel(u"T/°C")

    plt.xlabel('Zeit/s')

    plt.plot(t,np.array([np.mean(M)]*len(t)),color='darkorange')

    plt.savefig("Auswertung/Rauschmessung_stange_1.jpg")

    plt.figure()



#------------------------------------------------------------

#Umgebungsdruck

    print("==="*20)

    print("Spannungsverteilung 2")

    print('m='+str(mSt)+' V, std='+str(sSt)+' V, err='+str(errSt)+' V')

    

    #Histogramm

    weightsSt=np.ones_like(M)/float(len(M))

    plt.hist(St, 7, normed=0,weights=weightsSt)

    x=np.arange(min(St), max(St), 0.0001)

    #plt.plot(x, gauss(x, np.average(St), np.std(St)))

    plt.title(u"Histogramm der Ruhelage der Pendelstange 2 \n $\sigma$=" + str(sSt)+" V, $\mu$="+str(mSt)+" V, $\sigma_{\mu} = $" +str(errSt)+" V")

    plt.xlabel(u"U/V")

    plt.ylabel('Relatives Vorkommen')

    plt.savefig("Auswertung/Rauschmessung_stange_2_hist.jpg")

    plt.axvline(np.mean(St),color='red',linestyle='--')

    plt.figure()

    

    #Datenplot

    plt.plot(t, St)

    plt.ylabel(u"U/V")

    plt.xlabel('Zeit/s')

    plt.title("Rauschmessung Ruhelage Pendel 2")

    plt.plot(t,np.array([np.mean(St)]*len(t)),color='darkorange')

    plt.savefig("Auswertung/Rauschmessung_stange_2.jpg")