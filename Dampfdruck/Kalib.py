#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 27 19:27:20 2018

@author: Mate
"""

"""
Output von Rauschmessung.py:
    Temperaturverteilung bei RT:
    m=25.2 °C , s=0.0 °C , err=0.1 °C
    ============================================================
    Druckverteilung bei  RT:
    m=984.8 hPa, s=0.18 hPa, err=0.01 hPa
    
    Also Druck wird wegen 15000hPa*0.05 = 0.75 hPa (Auflösungsvermögen) auf 2 Nachkommastellen und
    T wird wegen sigma_gerät = +-0.2 °C oder +- 0.4°C auf 1 Nachkommastelle berechnet/aufgerundet
"""
sigmaT_Rausch = 0.1 #siehe obiges Kommentar
sigmaT_Gerät_kl_70 = 0.2 #kleiner oder
sigmaT_Gerät_gr_70 = 0.4 #größer als 70 °C

import matplotlib.pyplot as plt
import numpy as np
import scipy
import praktikum.cassy as cassy
import praktikum.analyse as ana

def gauss(x, m, s):
    return(np.power(scipy.pi*2*s**2, -0.5)*np.exp(-(x-m)**2/(2*s**2)))

def round_up(x, digits): #für obere Abschätzung der Fehler
    if digits!=0: return np.ceil(x*10**digits)/10**digits
    else: return np.ceil(x)
    
def T_th(T_gem, a, b, sT_gem, sa, sb):
    T = (T_gem-b-50)/a + 50
    #Gauss auf T_th = (T_gem-b-50)/a + 50) bzgl. a, b, T_gem:
    sT = np.sqrt((T_gem/a**2)**2 * sa**2 + a**-2 * sb**2 + a**-2 * sT_gem**2)
    return np.array([T, sT])

t, p, T = [], [], [] #Temperatur in °C!
data = cassy.CassyDaten("Kalibrierung 1 + Dichtigkeit.lab")
t = data.messung(1).datenreihe("t").werte
T = data.messung(1).datenreihe("&J_A11").werte
p = data.messung(1).datenreihe("p_B1").werte
sigmaT=np.std(T)
#0.01: wegen Output von Rauschmessung.py (siehe oben), 0.75: Auflösungsvermögen
sigmap=np.sqrt(0.01**2+0.75**2)
errT=sigmaT/np.sqrt(len(T))
errp=sigmap

#Temperaturkalibrierung:

#T=0°C
N_Tslice=2000
t_sliced=t[N_Tslice:]
T=T[N_Tslice:]
print("Temperaturkalibrierung")
T0=np.average(T)
sigmaT0 = np.std(T)/len(T)
print("s=" + str(round_up(sigmaT0,1))+" °C, m="+str(np.round(T0,1))+" °C")
plt.title("Temperatur bei T=0°C\n $\sigma_{\mu}$=" + str(round_up(sigmaT0,1))+"°C, $\mu$="+str(np.round(T0,1))+"°C (ab t=200 s berechnet)")
plt.plot(np.arange(min(t), max(t), 0.1),np.ones(len(t))*T0, color="green")
#Rückkehr an die ursprüngliche Daten wegen des plots...
T = data.messung(1).datenreihe("&J_A11").werte
plt.plot(t, T)
T=T[N_Tslice:]
#plot
plt.axvline(x=200, color="red", linestyle = "--")
plt.xlabel('Zeit/s')
plt.ylabel(u'Temperatur/°C')
plt.savefig("Images/Kalib_T_0.pdf")
plt.figure()

x=np.arange(np.average(T)-0.4, np.average(T)+0.4, 0.01)
plt.hist(T, normed=True)
plt.plot(x, gauss(x, np.average(T), np.std(T)))
plt.title("Streuung der gemessenen Temperaturwerte bei T=0°C\n $\sigma_{\mu}$=" + str(round_up(sigmaT0,1))+"°C, $\mu$="+str(np.round(T0,1))+"°C (ab t=200 s berechnet)")
plt.ylabel('Relatives Vorkommen')
plt.xlabel(u'Temperatur/°C')
plt.savefig("Images/Kalib_T_0_histo.pdf")
plt.figure()

print("\nDichtigkeitsmessung")
print("Fehler auf den gemessenen Druck")
print("s**2 = s_Ger**2+s_Rausch**2=" + str(np.round(sigmap**2, 2)) + " hPa")
plt.plot(t, p)
m,em,b,eb,chiq,corr = ana.lineare_regression(t, p, sigmap*np.ones(len(p)))
print("\nErgebnis der Regression auf die Dichtigkeitskurve:")
print('m = (%g +- %g) hPa/s,   b = (%g +- %g) hPa,  chi2/dof = %g / %g  corr = %g' % (m, em, b, eb, chiq, len(t)-2, corr))
print('Verlust: ({} +-{}) mbar/min\n'.format(np.round(60*m, 2),round_up(60*em,2))) #so genehmigt

plt.plot(t, m*t+b, color='green')
plt.title("Dichtigkeitsmessung")
plt.xlabel('Zeit/s')
plt.ylabel(u'Druck/hPa')
plt.savefig("Images/Kalib_Dichtigkeitsmessung.pdf")
plt.figure()

#T=100°C:
data = cassy.CassyDaten("Temperatur, Druck bei 100.lab")
t = data.messung(1).datenreihe("t").werte
p = data.messung(1).datenreihe("p_B1").werte
T = data.messung(1).datenreihe("&J_A11").werte

N_slice = [240, 1500] #Betrachte dieses Werteintervall, weil außerhalb von diesem komische Sachen passieren...
t_sliced = t[240:1500]
T_sliced = T[240:1500]
p_sliced = p[240:1500]
T100 = np.average(T_sliced)
sigmaT100 = np.std(T_sliced)/len(T_sliced)
sigmap_sliced = np.std(p_sliced)/len(p_sliced)

x = np.arange(min(t),max(t), 0.01)
plt.plot(t, p)
plt.plot(x, np.ones(len(x))*np.average(p), c="green")
plt.axvline(x=24, color="red", linestyle = "--")
plt.axvline(x=150, color="red", linestyle = "--")
plt.gcf().subplots_adjust(top=0.8)
plt.title("Dampfdruck bei T=100°C\n $\mu$ = " + str(np.round(np.average(p),2)) + " hPa, $\sigma_{\mu} = $"+str(round_up(sigmap_sliced,2)) + " hPa \n (im Intervall t $\in$ [24 s; 150 s])")
plt.xlabel('Zeit/s')
plt.ylabel(u'Druck/hPa')
plt.savefig("Images/Kalib_Dampfdruck_T_100.pdf")
plt.figure()

plt.plot(t, T)
x = np.arange(min(t),max(t), 0.01)
plt.plot(x, np.ones(len(x))*T100, color="green")
plt.title("Temperatur bei T=100°C\n $\mu$ ="+ str(np.round(T100,1))+" °C, " + "$\sigma_{\mu}$ = " + str(round_up(sigmaT100,2)) + " °C \n (im Intervall t $\in$ [24 s; 150 s])")
plt.axvline(x=24, color="red", linestyle = "--")
plt.axvline(x=150, color="red", linestyle = "--")
plt.gcf().subplots_adjust(top=0.8)
plt.xlabel('time/s')
plt.ylabel(u'T/°C')
plt.savefig("Images/Kalib_T_100.pdf")
plt.figure()

x = np.arange(T100-0.75,T100+0.75, 0.01)
plt.hist(T, normed=True)
plt.plot(x, gauss(x, T100, sigmaT100*len(T_sliced)))
plt.gcf().subplots_adjust(top=0.8)
plt.title("Streuung der gemessenen Temperaturwerte bei T=100°C\n $\mu$ ="+ str(np.round(T100,1))+" °C, " + "$\sigma_{\mu}$ = " + str(round_up(sigmaT100,2)) + " °C \n (im Intervall t $\in$ [24 s; 150 s])")
plt.ylabel('rel')
plt.xlabel(u'T/°C')
plt.savefig("Images/Kalib_T_histo_100.pdf")
plt.figure()

print("Ergebnis auf die Korrekturgleichung:")
#Die Fehler kommen aus der Rauschmessung, aus der Fluktuation der Werte bei Temperatur T und aus dem Fehler des Geräts:
#sigmaT_tot**2 = sigmaT_Rausch**2 + sigmaT_Gerät**2 + sigmaT0**2
sigmaT0_tot = np.sqrt(sigmaT_Rausch**2+sigmaT_Gerät_kl_70**2+sigmaT0**2)
sigmaT100_tot = np.sqrt(sigmaT_Rausch**2+sigmaT_Gerät_gr_70**2+sigmaT100**2)
print("T bei T=0°C : T_gem = " + str(np.round(T0,1))+" °C +-"+str(round_up(sigmaT0_tot,1)) + " °C")
print("T bei T=100°C : T_gem = " + str(np.round(T100,1))+" °C +-"+str(round_up(sigmaT100_tot,1)) + " °C")
      
#Führe lineare Reg. bzgl. T_gem durch, dann stelle auf T_th um, weil wir bei der lin_reg nur Fehler auf y angeben können...
a,ea,b,eb,chiq,corr = ana.lineare_regression(np.array([0, 100])-50,np.array([T0, T100])-50, np.array([sigmaT0_tot,sigmaT100_tot]))
print("\nDamit ergibt sich die Korrekturgleichung:")
print("T_gem = "+ str(a)+"*(T_th-50) + " + str(b)+" + 50")
print("oder auf T_th umgestellt:")
print("T_th = (T_gem-b-50)/a + 50")
print("mit a = " + str(np.round(a, 1)) + " C^-1 +- " + str(round_up(ea,1)) + " °C^-1")
print("und mit b = " + str(np.round(b, 1)) + " °C +- " + str(round_up(eb,1)) + " °C")

T_gem=99.8
T_theo, sT = T_th(T_gem,a,b,sigmaT_Gerät_gr_70,ea,eb)
print("damit haben wir T_th(99.8°C) = " + str (np.round(T_theo,1)) + " °C +-" + str(round_up(sT,1)) + " °C")
T_gem=2.5
T_theo, sT = T_th(T_gem,a,b,sigmaT_Gerät_kl_70,ea,eb)
print("und T_th(2.5°C) = " + str (np.round(T_theo,1)) + " °C +-" + str(round_up(sT,1)) + " °C")
