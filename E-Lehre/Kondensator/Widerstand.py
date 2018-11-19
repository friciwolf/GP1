# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 16:37:37 2018

@author: Patrick, Denise
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy
import praktikum.cassy as cassy
from praktikum import analyse as ana

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

def sigma_U_sys(U,U_B):
    return (0.01*U+0.005*U_B)/np.sqrt(3)
    
def sigma_I_sys(I,I_B):
    return (0.02*I+0.005*I_B)/np.sqrt(3)

plt.close('all')

U_array=[]
I_array=[]
U_err_array=[]
I_err_array=[]
U_sigma=[]
I_sigma=[]
U_B=20
I_B=0.6 #TODO: Richtig?


sigmaU = U_B / 4096. / np.sqrt(12.)
sigmaI = I_B / 4096. / np.sqrt(12.)

for i in range(1,6):
    #Daten Einlesen
    t, U, I = [], [], [] #Zeit, Spannung, Strom
    data = cassy.CassyDaten("Widerstand_2.lab")
    t = data.messung(i).datenreihe("t").werte
    U = data.messung(i).datenreihe("U_B1").werte
    I = data.messung(i).datenreihe("I_A1").werte
    
    #Mittelwert, Standardabweichung, Fehler auf den Mittelwert
    mU,sU,errU=round_good(np.mean(U),np.std(U),np.std(U)/np.sqrt(len(U)))
    mI,sI,errI=round_good(np.mean(I),np.std(I),np.std(I)/np.sqrt(len(I)))
    
    if __name__=='__main__':
        print('U-werte')
        print('m='+str(mU)+u' V, std='+str(sU)+u' V, err='+str(errU)+u' V')
        print('I-werte')
        print('m='+str(mI)+u' A, std='+str(sI)+u' A, err='+str(errI)+u' A')
        
        #Histogramme
        plt.hist(U, normed=True)
        plt.title('Histogramm der Spannungen um {}'.format(mU))
        plt.xlabel('Spannung in Volt')
        plt.ylabel('rel. Häufigkeiten')
        #plt.savefig('Images/Widerstand_U_hist.jpg')
        plt.figure()
        
        plt.hist(I, normed=True)
        plt.title('Histogramm der Ströme um {}'.format(mI))
        plt.xlabel('Strom in Ampére')
        plt.ylabel('rel. Häufigkeiten')
        #plt.savefig('Images/Widerstand_I_hist.jpg')
        plt.figure()
    
    U_array.append(mU)
    I_array.append(mI)
    if errU>sigmaU:
        U_err_array.append(errU)
    else:
        U_err_array.append(sigmaU)
    if errI>sigmaI:
        I_err_array.append(errI)
    else:
        I_err_array.append(sigmaI)
    U_sigma.append(sU)
    I_sigma.append(sI)
    #Schleifenende

I_array=np.array(I_array)
U_array=np.array(U_array)

#Widerstand+statistischer Fehler
R, eR, b, eb, chiq, corr=ana.lineare_regression_xy(I_array,U_array,np.array(I_err_array),np.array(U_err_array))
x,R,errR_stat=round_good(0.0,R,eR)

#systematischer Fehler
RUp, eRx, bUp, ebx, chiqUp, corrx=ana.lineare_regression_xy(I_array,U_array+sigma_U_sys(U_array,U_B),np.array(I_err_array),np.array(U_err_array))
RUm, eRx, bUm, ebx, chiqUm, corrx=ana.lineare_regression_xy(I_array,U_array-sigma_U_sys(U_array,U_B),np.array(I_err_array),np.array(U_err_array))
RIp, eRx, bIp, ebx, chiqIp, corrx=ana.lineare_regression_xy(I_array+sigma_I_sys(I_array,I_B),U_array,np.array(I_err_array),np.array(U_err_array))
RIm, eRx, bIm, ebx, chiqIm, corrx=ana.lineare_regression_xy(I_array-sigma_I_sys(I_array,I_B),U_array,np.array(I_err_array),np.array(U_err_array))

eR_U_sys=(abs(RUp-R)+abs(RUm-R))/2
eR_I_sys=(abs(RIp-R)+abs(RIm-R))/2

errR_sys=np.sqrt(eR_U_sys**2+eR_I_sys**2)
R,errR_sys,errR_stat=round_good(R,errR_sys,errR_stat)

#Statistische Fehler auf Einzelmessungen
U_err_stat=np.mean(U_sigma)
I_err_stat=np.mean(I_sigma)
x,x,U_err_stat=round_good(0,0,U_err_stat)
x,x,I_err_stat=round_good(0,0,I_err_stat)

if __name__=='__main__':
    #Plotten
    #plt.plot(I_array,U_array,'ko',label='Mittelwerte U,I, unverschoben')
    plt.errorbar(I_array,U_array,yerr=U_err_array,xerr=I_err_array,fmt='ko',label='Mittelwerte U,I, unverschoben',capsize=0.2)
    x=np.linspace(min(I_array),max(I_array),100)
    plt.plot(x,R*x+b, color='red',label=u'Regression Rx+b mit R={}, b={},Chi²/f={}'.format(R,b,chiq/len(U_array)))
    plt.plot(x,RUp*x+bUp, color='lightblue',label=u'U+ mit R+={}, b+={},Chi²/f={}'.format(RUp,bUp,chiqUp/len(U_array)))
    plt.plot(x,RUm*x+bUm, color='lightgreen',label=u'U- mit R-={}, b-={},Chi²/f={}'.format(RUm,bUm,chiqUm/len(U_array)))
    '''
    #alternative Darstellung
    plt.plot(x,RUp*x+bUp, color='lightblue',label = "Regression U+")
    plt.plot(x,RUm*x+bUm,'--', color='lightblue', label = "Regression U-")
    plt.plot(x,RIp*x+bIp, color='lightgreen',label = "Regression I+")
    plt.plot(x,RIm*x+bIm,'--', color='lightgreen', label = "Regression I-")
    '''
    plt.xlabel('I in A')
    plt.ylabel('U in V')
    plt.legend()
    #plt.savefig('Images/Widerstand_Regression.jpg')
    plt.figure()
    
    #Residuen
    plt.title('Residuenplot - Abweichung von der Geraden U-f(x)')
    plt.errorbar(I_array,U_array-(R*I_array+b),yerr=np.sqrt(np.array(U_err_array)**2+np.array(I_err_array)**2),fmt='ko',capsize=3)
    plt.xlabel('I in A')
    plt.ylabel('U-(R*I+b) in V')
    #plt.savefig('Widerstand_Res')
    
    #Resümee
    print(u'Eine lineare Regression ergibt R={} +- {} mit einem Chi²/f={}/{}'.format(R,eR,chiq,len(U_array)))
    print('Insgesamt erhält man R = {} +- {} +- {} Ohm'.format(R,errR_stat,errR_sys))
