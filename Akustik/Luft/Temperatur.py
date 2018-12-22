# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 20:10:54 2018

@author: Anwender
"""


import matplotlib.pyplot as plt
import numpy as np
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

dateien=['messung_vorher','_mitte','_ende']
names=['vor den Messungen','in der Mitte','nach den Messungen']
ts=[]
Ts=[]
T_errs=[]

dig=0.1/np.sqrt(12)

for i,n in enumerate(dateien):
    data = cassy.CassyDaten("Temperatur"+n+".lab")
    t = data.messung(1).datenreihe("t").werte
    T = data.messung(1).datenreihe("&J_A11").werte

    #----------------------------#
    #Auswertung der Rauschmessung#
    #----------------------------#
    
    #Mittelwerte, Standardabweichungen und Fehler
    mT,sT,errT=round_good(np.average(T),np.std(T),np.std(T)/np.sqrt(len(T)))
    
    if __name__ == "__main__": 
        
    #------------------------------------------------------------
    
        print("Temperatur",i+1)
        print('m='+str(mT)+u' °C, std='+str(sT)+u' °C, err='+str(errT)+u' °C')
        
        #Histogramm
        weightsT=np.ones_like(T)/float(len(T))
        plt.hist(T,7, normed=0, weights=weightsT)
        x=np.arange(np.average(T)-0.03, np.average(T)+0.03, 0.001)
        #plt.plot(x, gauss(x, np.average(T), np.std(T))
        plt.title(u"Histogramm der Temperatur {}\n $\sigma$=".format(names[i]) + str(sT)+"$^{\circ}C$, $\mu$="+str(mT)+ "$^{\circ}C$, $\sigma_{\mu} = $"+str(errT)+"$^{\circ}C$")
        plt.xlabel(u"T/°C")
        plt.ylabel('Relatives Vorkommen')
        plt.savefig("Images/Temperatur{}.jpg".format(n))
        plt.axvline(np.mean(T),color='red',linestyle='--')
        plt.figure()
        
        #Datenplot
        plt.plot(t, T)
        plt.title(u"Rauschmessung der Temperatur {}".format(names[i]))
        plt.ylabel(u"T/°C")
        plt.xlabel('Zeit/s')
        plt.plot(t,np.array([np.mean(T)]*len(t)),color='darkorange')
        plt.savefig("Images/Rauschmessung{}.jpg".format(n))
        plt.figure()
    
    ts.append(t)
    Ts.append(mT)
    T_errs.append(max(errT,dig))

#Versuch einer linearen Regression - sollen wir soetwas machen oder irgendwie die Temperaturen bei den einzelnen Versuchen feststellen?
#t1=17:05,t2=18:31, t3=19:45
zeiten=np.array([0,86,160]) #min
Ts=np.log(np.array(Ts))
T_errs=np.array(T_errs)/np.array(Ts)
c,c_err,b,b_err,chiq,corr=analyse.lineare_regression(zeiten,Ts,T_errs)
print('\nErwärmung: {2}°C/min'.format(*round_good(0,0,c)))

x=np.arange(min(zeiten),max(zeiten),0.1)
plt.errorbar(zeiten, Ts,yerr=T_errs,fmt='bo') #eher exponentiell?
plt.plot(x,c*x+b,color='red')
plt.xlabel('t/min')
plt.ylabel(u'T/°C')
plt.title('Temperaturverlauf während der Versuchsreihe')
print(u'Chi² =',chiq)
