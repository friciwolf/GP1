# -*- coding: utf-8 -*-

import numpy as np
import praktikum.analyse as analyse
import praktikum.cassy as cassy
import Rauschmessung

offsetM=Rauschmessung.mM
offsetSt=Rauschmessung.mSt

l=0.68
r=0.08035/2
errl_stat=0.001
errr_stat=0.0001
errl_sys=0.0003+0.0002*1 # Güteklasse II:a+b*L mit a=0.3mm, b=0.2mm/m, L= auf ganzen Meter gerundete zu messende Länge
T=[]
errT=[]

for i in range(1,4):
    data = cassy.CassyDaten("Erdbeschleunigung/perfekte_anpassung_messung_"+str(i)+".lab")
    t = data.messung(1).datenreihe("t").werte
    M = data.messung(1).datenreihe("U_A1").werte-offsetM #Pendelkörper (Masse)
    St = data.messung(1).datenreihe("U_B1").werte-offsetSt #Stange allein
    
    if __name__=='__main__':
        #Rohdatenplot bis Stange zu gedämpft
        bis=np.argmax(t>100)
        plt.plot(t[:bis],M[:bis], color='black')
        plt.plot(t[:bis],St[:bis],color='red')
        plt.title('Schwingung mit körperloser Frequenz {}'.format(i))
        plt.grid()
        plt.savefig('Images/Erdbeschleunigung_Roh_'+str(i)+'.jpg')
        if i!=3:
            plt.figure()
    
    #Nullstellen zählen
    count=0
    for i,volt in enumerate(M):
        if volt!=M[-1] and ((volt<=0 and M[i+1]>0) or (volt>=0 and M[i+1]<0)):
            count+=1
            ende=i
            if count==1:
                start=i
    if count%2==1:
        count-=1
        ende-=1
        
    #Periodendauer berechnen
    errt=(t[-1]-t[0])/len(t)
    T.append(2*(t[ende]-t[start])/count)
    errT.append(2*np.sqrt(2)*errt/count) #TODO: Richtig?

T,errT=analyse.gewichtetes_mittel(np.array(T),np.array(errT))

#Frequenz und Erdbeschleunigung + Fehler
f=1/T
errf=f*errT/T
print('\nFrequenz: f=({1} +- {2}) Hz'.format(*Rauschmessung.round_good(0,f,errf)))

g=4*np.pi**2*f**2*l*(1+r**2/l**2/2)
errg_stat=np.sqrt(((g/l-4*np.pi**2*f**2*r**2/l**2)*errl_stat)**2+(2*g/f*errf)**2) #TODO: vlt.genauer?
errg_sys=(g/l-4*np.pi**2*f**2*r**2/l**2)*errl_sys
if errg_sys>errg_stat: 
   g,errg_stat,errg_sys=Rauschmessung.round_good(g,errg_stat,errg_sys)
else: 
   g,errg_sys,errg_stat=Rauschmessung.round_good(g,errg_sys,errg_stat)
print('Erdbeschleunigung: g=({} +- {} +- {})m/s²'.format(g,errg_stat,errg_sys))
