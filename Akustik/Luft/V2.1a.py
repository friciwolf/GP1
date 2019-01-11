# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 23:29:44 2018

@author: Anwender
"""

import Kalibrierung

R=[]
dt=[]
dt_err=[]
R_err=float(Kalibrierung.R_err[0])
dt_errdig=0.00000025 #s
v_0T_0=np.sqrt(8.3145*7/5/28.984*10**3)

ns=[]
#Daten einlesen
for i in range(1,7):
    data = cassy.CassyDaten('V2.1A.lab')
    R_i = np.array(data.messung(i).datenreihe('R_B1').werte)
    dt_i = list(data.messung(i).datenreihe('&Dt_A1').werte)
    
    mR,sR,errR=Temperatur.round_good(R_i[2],0,R_err) #erster R_i-Wert manchmal verschoben
    
    mdt=np.mean(dt_i)
    
    #Abweichungen aussortieren #TODO: letzte Messreihe komisch; vlt. lieber Std benutzen, statt zu sortieren? 
    while np.argmax(abs(np.array(dt_i)-mdt)>0.001*mdt)!=0:
        dt_i.pop(np.argmax(abs(np.array(dt_i)-mdt)>0.001*mdt))
    
    mdt,sdt,errdt=Temperatur.round_good(np.mean(dt_i),np.std(dt_i),dt_errdig)
    
    R.append(mR)
    dt.append(mdt)
    dt_err.append(max(sdt/(len(dt_i)-1),dt_errdig))
    ns.append(len(dt_i))

#Lineare Regression    
dt=np.array(dt)
#dt[-1]=0.0006200#6205
#dt_err[-1]=dt_errdig
dt_err=np.array(dt_err)
L,L_estat,L_esys=Kalibrierung.Länge(np.array(R),R_err*np.ones(len(R)))
v,v_err,s0,s0_err,chiq,corr=analyse.lineare_regression_xy(dt,L,dt_err,L_estat)

#systematischer Fehler durch Verschiebemethode
vplus,v_errplus,s0plus,s0_errplus,chiqplus,corr=analyse.lineare_regression_xy(dt,L+L_esys,dt_err,L_estat)
vmin,v_errmin,s0min,s0_errmin,chiqmin,corr=analyse.lineare_regression_xy(dt,L-L_esys,dt_err,L_estat)
v_esys=(abs(vplus-v)+abs(vmin-v))/2

v_esys,v,v_err=Temperatur.round_good(v_esys,v,v_err)
chiq,s0,s0_err=Temperatur.round_good(chiq,s0,s0_err)
print('\nv_Schall=({} +- {} +- {})cm/s\nAbstand von der Quelle:({} +- {})cm\nChi²/f={}/{}'.format(v,v_err,v_esys,s0,s0_err,chiq,len(dt)-2)) #dL/dt #cm/s

#Plot
plt.title('Laufzeit gegen Laufstrecke')
plt.plot(dt,L,'ko')
x=np.arange(min(dt),max(dt),0.0000001)
plt.plot(x,v*x+s0,'r-', label='$v\cdot x+s_0=({}\pm{})cm/s*t+({}\pm{})cm $ \n Chi²/f={}/{} ≈ {}'.format(v,v_err,s0,s0_err, chiq,len(R)-2,np.round(chiq/(len(R)-2),2)))
plt.xlabel(r'Laufzeit $\Delta$t/s')
plt.ylabel('Strecke L/cm')
plt.legend(loc=2)
plt.savefig('Images/VA_Regression.pdf')
plt.figure()

#Residuum
plt.title('Residuum V2.1A')
plt.plot(x,np.zeros(len(x)),'r--')
plt.errorbar(dt,L-(v*dt+s0),yerr=np.sqrt(v**2*dt_err**2+L_estat**2),fmt='ko')
plt.xlabel(r'$\Delta$t/s')
plt.ylabel('L-(v*dt+s0)/cm')
plt.savefig('Images/VA_Residuum.pdf')

#Theoretische Geschwindigkeit
T,T_estat,T_esys=Temperatur.Temperatur_mit_Fehlern(82)
#v_theo=331.6+0.6*T
v_theo=v_0T_0*np.sqrt(T+273.15) #T in K
v_theo_estat=v_theo*0.5/(T+273.15)*T_estat
v_theo_esys=v_theo*0.5/(T+273.15)*T_esys

T,T_estat,T_esys=Temperatur.round_good(T,T_estat,T_esys)
v_theo,v_theo_estat,v_theo_esys=Temperatur.round_good(v_theo,v_theo_estat,v_theo_esys)
print('Theoretische Geschwindigkeit bei T=({} +- {} +- {})°C: ({} +- {} +- {})m/s'.format(T,T_estat,T_esys,v_theo,v_theo_estat,v_theo_esys))
