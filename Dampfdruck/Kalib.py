#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 27 19:27:20 2018

@author: Mate, Patrick
"""

"""
Daten von Rauschmessung.py:
        Temperaturverteilung
        mT=25.15797 °C, sT=0.04443 °C, errT=0.00077 °C
        ============================================================
        Druckverteilung
        mp=984.7955 hPa, sp=0.1838 hPa, errp=0.0032 hPa
    
    Druck:      statistisch: 1500hPa*0.05% = 0.75 hPa (Auflösungsvermögen) 
    Temperatur: statistisch: (120--20)°C*? = ?? °C    (Auflösungsvermögen nicht gegeben)
                sytematisch: sigma_gerät = +-0.2 °C oder +- 0.4°C 
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy
import praktikum.cassy as cassy
import praktikum.analyse as ana
import Rauschmessung

sigmaT_Gerät_kl_70 = 0.2 #kleiner oder
sigmaT_Gerät_gr_70 = 0.4 #größer als 70 °C
sigma_p_dig=1500*0.05/100 #Auflösung Druck


def gauss(x, m, s):
    return(np.power(scipy.pi*2*s**2, -0.5)*np.exp(-(x-m)**2/(2*s**2)))

def T_th(T_gem, a, b, sT_gem, sa, sb,sT_gemsys,sasys,sbsys):
    T = (T_gem-b-50)/a + 50
    #Gauss auf T_th = (T_gem-b-50)/a + 50) bzgl. a, b, T_gem:
    sT = np.sqrt(((T_gem-b-50)/a**2)**2 * sa**2 + a**-2 * sb**2 + a**-2 * sT_gem**2)
    sTsys = np.sqrt(((T_gem-b-50)/a**2)**2 * sasys**2 + a**-2 * sbsys**2 + a**-2 * sT_gemsys**2)
    return np.array([T, sT, sTsys])

#Erste Messung
t, p, T = [], [], [] #Temperatur in °C!
data = cassy.CassyDaten("Kalibrierung 1 + Dichtigkeit.lab")
t = data.messung(1).datenreihe("t").werte
T = data.messung(1).datenreihe("&J_A11").werte
p = data.messung(1).datenreihe("p_B1").werte


#-----------------------------------------#
#Temperaturkalibrierung und Gasdichtigkeit#
#-----------------------------------------#

#------------------------------------------------------------
#Gasdichtigkeit

#Erstauswertung
mp,sp,errp=Rauschmessung.round_good(np.mean(p),np.std(p),np.std(p)/np.sqrt(len(p)))
errstatp=np.sqrt((Rauschmessung.sp/np.sqrt(len(p)))**2+sigma_p_dig**2)
sigmap=np.sqrt((Rauschmessung.sp)**2+sigma_p_dig**2)
#errT=sT0/np.sqrt(len(T))
#errp=sigmap

#Plotten
print("\nDichtigkeitsmessung")
print("Fehler auf jeden gemessenen Druckwert")
print("s**2 = s_Gerät**2+s_Rausch**2=" + str(Rauschmessung.round_good(0.0,0.0,sigmap**2)[2]) + " hPa")
plt.plot(t, p)

#Regression
m,em,b,eb,chiq,corr = ana.lineare_regression(t, p, sigmap*np.ones(len(p)))
print("\nErgebnis der Regression auf die Dichtigkeitskurve:")
print('m = (%g +- %g) hPa/s,   b = (%g +- %g) hPa,  chi2/dof = %g / %g  corr = %g' % (m, em, b, eb, chiq, len(t)-2, corr))
m,x,em=Rauschmessung.round_good(m,0.0,em)
print('Verlust: ({} +-{}) mbar/s\n'.format(m,em)) #so genehmigt

plt.plot(t, m*t+b, color='green')
plt.title("Dichtigkeitsmessung\nVerlust: (0.010501 +-2.6e-05) mbar/s")
plt.xlabel('Zeit/s')
plt.ylabel(u'Druck/hPa')
plt.savefig("Images/Kalib_Dichtigkeitsmessung.pdf")
plt.figure()

#------------------------------------------------------------
#T=0°C
plt.plot(t, T)

N_Tslice=2000 #gleichmäßiger Bereich
t_sliced=t[N_Tslice:]
T=T[N_Tslice:]
print("\nTemperaturkalibrierung bei 0°C")

#Erstauswertung
mT0,sT0,errT0=Rauschmessung.round_good(np.mean(T),np.std(T),np.std(T)/np.sqrt(len(T)))
errstatT0=Rauschmessung.round_good(0,0,Rauschmessung.sT/np.sqrt(len(T)))[2]
errsysT0=sigmaT_Gerät_kl_70
print("T0=" + str(mT0)+" °C, s="+str(sT0)+" °C,err="+str(errT0)+" °C, errstat="+str(errstatT0)+' °C\n')

#Auftragen der Daten
plt.title("Temperatur bei T=0°C\n $\sigma_{\mu}$=" + str(errT0)+"°C, $\mu$="+str(mT0)+"°C (ab t=200 s berechnet)")
plt.plot(np.arange(min(t), max(t), 0.1),np.ones(len(t))*mT0, color="green")
plt.axvline(x=200, color="red", linestyle = "--")
plt.xlabel('Zeit/s')
plt.ylabel(u'Temperatur/°C')
plt.savefig("Images/Kalib_T_0.pdf")
plt.figure()

#Histogramm
x=np.arange(np.average(T)-0.4, np.average(T)+0.4, 0.01)
plt.hist(T,7, normed=True)
plt.plot(x, gauss(x, np.average(T), np.std(T)))
plt.title("Streuung der gemessenen Temperaturwerte bei T=0°C\n $\sigma_{\mu}$=" + str(errT0)+"°C, $\mu$="+str(mT0)+"°C (ab t=200 s berechnet)")
plt.ylabel('Relatives Vorkommen')
plt.xlabel(u'Temperatur/°C')
plt.savefig("Images/Kalib_T_0_histo.pdf")
plt.figure()

#------------------------------------------------------------
#T=100°C
#Zweite Messung
data = cassy.CassyDaten("Temperatur, Druck bei 100.lab")
t = data.messung(1).datenreihe("t").werte
p = data.messung(1).datenreihe("p_B1").werte
T = data.messung(1).datenreihe("&J_A11").werte

N_slice = [240, 1500] #Betrachte dieses Werteintervall, weil außerhalb von diesem komische Sachen passieren...
t_sliced = t[240:1500]
T_sliced = T[240:1500]
p_sliced = p[240:1500]
#sigmaT100 = np.std(T_sliced)/len(T_sliced)
#sigmap_sliced = np.std(p_sliced)/len(p_sliced)

#Erstauswertung
mT100,sT100,errT100=Rauschmessung.round_good(np.mean(T_sliced),np.std(T_sliced),np.std(T_sliced)/np.sqrt(len(T_sliced)))
mp100,sp100,errp100=Rauschmessung.round_good(np.mean(p_sliced),np.std(p_sliced),np.std(p_sliced)/np.sqrt(len(p_sliced)))
errstatT100=Rauschmessung.round_good(0,0,Rauschmessung.sT/np.sqrt(len(T_sliced)))[2]
errsysT100=sigmaT_Gerät_gr_70
print("\nTemperaturkalibrierung bei 100°C")
print("T100=" + str(mT100)+" °C, s="+str(sT100)+" °C,err="+str(errT100)+" °C, errstat="+str(errstatT100)+' °C\n')

#Plotten
#Druck
x = np.arange(min(t),max(t), 0.01)
plt.plot(t, p)
plt.plot(x, np.ones(len(x))*np.average(p), c="green")
plt.axvline(x=24, color="red", linestyle = "--")
plt.axvline(x=150, color="red", linestyle = "--")
plt.gcf().subplots_adjust(top=0.8)
plt.title("Dampfdruck bei T=100°C\n $\mu$ = " + str(mp100) + " hPa, $\sigma_{\mu} = $"+str(errp100) + " hPa \n (im Intervall t $\in$ [24 s; 150 s])") #Dampfdruck? Ist doch im Prinzip Umgebungsdruck...
plt.xlabel('Zeit/s')
plt.ylabel(u'Druck/hPa')
plt.savefig("Images/Kalib_Dampfdruck_T_100.pdf")
plt.figure()

#Temperatur
plt.plot(t, T)
x = np.arange(min(t),max(t), 0.01)
plt.plot(x, np.ones(len(x))*mT100, color="green")
plt.title("Temperatur bei T=100°C\n $\mu$ ="+ str(mT100)+" °C, " + "$\sigma_{\mu}$ = " + str(errT100) + " °C \n (im Intervall t $\in$ [24 s; 150 s])")
plt.axvline(x=24, color="red", linestyle = "--")
plt.axvline(x=150, color="red", linestyle = "--")
plt.gcf().subplots_adjust(top=0.8)
plt.xlabel('Zeit/s')
plt.ylabel(u'T/°C')
plt.savefig("Images/Kalib_T_100.pdf")
plt.figure()

#Temperaturhistogramm
x = np.arange(mT100-0.75,mT100+0.75, 0.01)
plt.hist(T,9, normed=True)
plt.plot(x, gauss(x, mT100, np.std(T_sliced)))
plt.gcf().subplots_adjust(top=0.8)
plt.title("Streuung der gemessenen Temperaturwerte bei T=100°C\n $\mu$ ="+ str(mT100)+" °C, " + "$\sigma_{\mu}$ = " + str(errT100) + " °C \n (im Intervall t $\in$ [24 s; 150 s])")
plt.ylabel('rel')
plt.xlabel(u'T/°C')
plt.savefig("Images/Kalib_T_histo_100.pdf")

#------------------------------------------------------------
#Kalibrierung
print("\n\nErgebnis auf die Korrekturgleichung:")
print("T bei T=0°C : T_gem = " + str(mT0)+" °C +-"+str(errT0) + " °C +-"+str(errsysT0)+' °C')
print("T bei T=100°C : T_gem = " + str(mT100)+" °C +-"+str(errT100) + " °C +-"+str(errsysT100)+' °C')
#Der statistische Fehler aus der Rauschmessung und der systematische des Geräts werden einzeln fortgepflanzt
    
#Führe lineare Reg. bzgl. T_gem durch, dann stelle auf T_th um, weil wir bei der lin_reg nur Fehler auf y angeben können...
a,ea,b,eb,chiq,corr = ana.lineare_regression(np.array([0, 100])-50,np.array([mT0, mT100])-50, np.array([Rauschmessung.sT,Rauschmessung.sT])) 
asys,easys,bsys,ebsys,chiq,corr = ana.lineare_regression(np.array([0, 100])-50,np.array([mT0, mT100])-50, np.array([errsysT0,errsysT100]))
print("\nErgebnis des linearen Modells:") #asys=a; bsys=b
a,ea,easys=Rauschmessung.round_good(a,ea,easys)
b,eb,ebsys=Rauschmessung.round_good(b,eb,ebsys)
print('a = (%g +- %g +- %g) hPa/s,   b = (%g +- %g +- %g) hPa' % (a, ea,easys, b, eb,ebsys))
print("\nDamit ergibt sich die Korrekturgleichung:")
print("T_gem = "+ str(a)+"*(T_th-50 °C) + " + str(b)+"°C + 50°C")
print("oder nach T_th umgestellt:")
print("T_th = (T_gem-b-50°C)/a + 50°C")
print("mit a = " + str(a) + " " + str(ea))
print("und mit b = " + str(b) + " °C +- " + str(eb) + " °C")

#Prüfung
T_theo, sT, sTsys= Rauschmessung.round_good(*T_th(mT100,a,b,Rauschmessung.sT,ea,eb,sigmaT_Gerät_gr_70,easys,ebsys))
print("\nSomit haben wir T_th(T100) = " + str (T_theo) + " °C +-" + str(sT) + " °C +-" + str(sTsys) + " °C")
T_theo, sT, sTsys= Rauschmessung.round_good(*T_th(mT0,a,b,Rauschmessung.sT,ea,eb,sigmaT_Gerät_kl_70,easys,ebsys))
print("und T_th(T0) = " + str (T_theo) + " °C +-" + str(sT) + " °C +-" + str(sTsys) + " °C")
