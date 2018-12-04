# -*- coding: utf-8 -*-
"""
Created on Sat Dec  1 14:31:28 2018

@author: Christine Falter
"""

import Rauschmessung
import Rauschmessung_gekoppelt
import Erdbeschleunigung
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy
import praktikum.cassy as cassy
import praktikum.analyse as analyse
import matplotlib.gridspec as gridspec


plt.close('all')

offset1=Rauschmessung_gekoppelt.mM
offset2=Rauschmessung_gekoppelt.mSt

g = Erdbeschleunigung.g
errg = np.sqrt((Erdbeschleunigung.errg_sys)**2 + (Erdbeschleunigung.errg_stat)**2) #als Systematischer
l_s=0.68
errl_s=0.001 #in m
m=1.0725+0.0872+0.1528
errm=0.0001 #in kg
l_F=np.array([16.5,16.5+12,16.5+12+12.4,16.5+12+12.4+12.5]) #in cm
errl_F=np.ones(4)*errl_s # in m

bis = 8000
name = ["Schwebung","Gleichsinnig","Gegensinnig"]

w_plus=[[0,0,0],[0,0,0],[0,0],[0]];w_minus=[[0,0,0],[0,0,0],[0,0],[0]]
errw=[[0,0],0,0,0]
ks=[];errks=[]

for i in range(1,5):
    plt.close('all')
    
    #Anzahl der Dateien
    path, dirs, files = next(os.walk("gekoppelt/pos_"+str(i)))
    file_count = len(files)
    t = np.empty((file_count,bis))
    U1 = np.empty((file_count,bis))
    U2 = np.empty((file_count,bis))
    #Daten einlesen
    data = cassy.CassyDaten("gekoppelt/pos_"+str(i)+"/schweb_pos_"+str(i)+".lab")
    t1 = data.messung(1).datenreihe("t").werte
    t[0] = t1[:bis]
    U1_1 = data.messung(1).datenreihe("U_A1").werte-(offset1+(i-1)*0.05) #TODO: Systematik im offset?
    U1[0] = U1_1[:bis]
    U1_2 = data.messung(1).datenreihe("U_B1").werte-offset2 #Pendel2
    U2[0] = U1_2[:bis]

    if file_count>1:
        data = cassy.CassyDaten("gekoppelt/pos_"+str(i)+"/gleich_pos_"+str(i)+".lab")
        t2 = data.messung(1).datenreihe("t").werte
        U2_1 = data.messung(1).datenreihe("U_A1").werte #Pendel1
        U2_1 = U2_1-np.mean(U2_1) #TODO: OK?
        U1[1] = U2_1[:bis]
        U2_2 = data.messung(1).datenreihe("U_B1").werte #Pendel2
        U2_2 = U2_2-np.mean(U2_2) #TODO: OK?
        U2[1] = U2_2[:bis]
        t[1] = t2[:bis]
        
    if file_count==3:
        data = cassy.CassyDaten("gekoppelt/pos_"+str(i)+"/gegen_pos_"+str(i)+".lab")
        t3 = data.messung(1).datenreihe("t").werte
        U3_1 = data.messung(1).datenreihe("U_A1").werte #Pendel1
        U3_1 = U3_1-np.mean(U3_1) #TODO: OK?
        U1[2] = U3_1[:bis]
        U3_2 = data.messung(1).datenreihe("U_B1").werte #Pendel2
        U3_2 = U3_1-np.mean(U3_1)
        U2[2] = U3_2[:bis]
        t[2] = t3[:bis]
    
    print(len(t1))
        
    #figure()
    #gs = gridspec.GridSpec(file_count, 1)
    for j in range(0,file_count):
        plt.figure()
        gs = gridspec.GridSpec(2, 1)
        ax1 = plt.subplot(gs[0, 0])
        plt.plot(t[j],U1[j],color='black')
        plt.title("Messung "+str(i)+": "+name[j])
        plt.ylabel(u"U1/V")
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.grid()
        
        ax2 = plt.subplot(gs[1, 0],sharex = ax1)
        plt.plot(t[j],U2[j],color='red')
        plt.ylabel(u"U2/V")
        plt.xlabel(u"t/s")
        ax2.yaxis.tick_right()
        plt.grid()
        plt.subplots_adjust(hspace=.0)
        
        plt.savefig("Images/Messung "+str(i)+"_"+name[j]+".jpg")
        
        #Nullstellen zählen
        count=0
        M=U1[j]
        for k,volt in enumerate(M):
            if k!=len(M)-1 and ((volt<=0 and M[k+1]>0) or (volt>=0 and M[k+1]<0)):
                count+=1
                if count==1:
                    a=(t[j][k+1]-t[j][k])/(M[k+1]-M[k]) #Die Zeiten sind die Nullstellen der Geraden durch die beiden Punkten um die Nullstelle
                    start=-M[k]/a+t[j][k]
                if count%2==1:
                    ende=-M[k]/a+t[j][k]
                #plt.plot(t[k],volt,'bo',color='lightblue')
        if count%2==1:
            count-=1
        else:
            count-=2
            
        print('Nullstellenzahl:',count)
        
        #Periodendauer berechnen
        errt=(t[-1]-t[0])/(len(t)-1)
        T[i-1][j]=(2*(ende-start)/count)
        #errT.append(2*np.sqrt(2)*errt/count)
        
        #Fourieranalyse
        w,A = analyse.fourier_fft(t[j],U1[j])
        ind_w1=analyse.peak(w,A,w[np.argmax(w>0.5)],w[np.argmax(w>0.6)])
        ind_w2=analyse.peak(w,A,w[np.argmax(w>0.7)],w[np.argmax(w>0.8)])
        
        range1=[(0.6,0.63),(0.6,0.63),(0.6,0.63),(0.5,0.64)]
        range2=[(0.6,0.63),(0.63,0.64),(0.64,0.7),(0.69,0.8)]
        ind_w1=list(A).index(max(A[np.argmax(w>range1[i-1][0]):np.argmax(w>range1[i-1][1])]),np.argmax(w>range1[i-1][0]))
        ind_w2=list(A).index(max(A[np.argmax(w>range2[i-1][0]):np.argmax(w>range2[i-1][1])]),np.argmax(w>range2[i-1][0]))
        
        plt.figure()
        plt.title('Fourieranalyse von Messung '+str(i)+": "+name[j])
        plt.plot(w,A/max(A))
        plt.ylabel(u"rel. Häufigkeit")
        plt.xlabel(u"w/Hz")
        plt.axvline(x=w[int(ind_w1)], color="darkred", linestyle = "--") 
        plt.text(w[int(ind_w1)],0.7,'max_1:{} Hz'.format(np.round(w[int(ind_w1)],4)))
        plt.axvline(x=w[int(ind_w2)], color="darkred", linestyle = "--") 
        plt.text(w[int(ind_w2)],0.6,'max_2:{} Hz'.format(np.round(w[int(ind_w2)],4)))
        plt.xlim(0,1.3)
        
        w_plus[i-1][j]=w[ind_w1] #Schweb, gleich, gegen für j=0,1,2
        w_minus[i-1][j]=w[ind_w2] #Kopplung: i=1 schwach, i=4 stark    
        
        if i==1 and j>0:
            errw[i-1][j-1]=(max(w)-min(w))/(len(w)-1)
    #Bestimmung des Kopplungsgrades in pos_1 aus gleich und gegen
    if i==1:
        k=(w_minus[i-1][2]**2-w_plus[i-1][1]**2)/(w_minus[i-1][2]**2+w_plus[i-1][1]**2)
        errk=(2*(w_minus[i-1][2]/(w_minus[i-1][2]**2+w_plus[i-1][1]**2)-w_minus[i-1][2]*(w_minus[i-1][2]**2-w_plus[i-1][1]**2)/(w_minus[i-1][2]**2+w_plus[i-1][1]**2)**2)*errw[i-1][2-1])**2
        errk+=(-2*(w_plus[i-1][1]/(w_minus[i-1][2]**2+w_plus[i-1][1]**2)-w_plus[i-1][1]*(w_minus[i-1][2]**2-w_plus[i-1][1]**2)/(w_minus[i-1][2]**2+w_plus[i-1][1]**2)**2)*errw[i-1][1-1])**2
        errk=np.sqrt(errk)
        x,k,errk=Rauschmessung.round_good(0,k,errk)
        print(str(i),'Kopplungsgrad:',k,'+-',errk)
        ks.append(k)
        errks.append(errk)

    #Bestimmung des Kopplungsgrades in pos_2-4 aus Schwebung
    if i>1: #TODO: Bei i=2 auch gleich und gegen auswerten?
        errw[i-1]=(max(w)-min(w))/(len(w)-1)
        k=(w_minus[i-1][0]**2-w_plus[i-1][0]**2)/(w_minus[i-1][0]**2+w_plus[i-1][0]**2)
        errk=(2*(w_minus[i-1][0]/(w_minus[i-1][0]**2+w_plus[i-1][0]**2)-w_minus[i-1][0]*(w_minus[i-1][0]**2-w_plus[i-1][0]**2)/(w_minus[i-1][0]**2+w_plus[i-1][0]**2)**2)*errw[i-1])**2
        errk+=(-2*(w_plus[i-1][0]/(w_minus[i-1][0]**2+w_plus[i-1][0]**2)-w_plus[i-1][0]*(w_minus[i-1][0]**2-w_plus[i-1][0]**2)/(w_minus[i-1][0]**2+w_plus[i-1][0]**2)**2)*errw[i-1])**2
        errk=np.sqrt(errk)
        x,k,errk=Rauschmessung.round_good(0,k,errk)
        print(str(i),'Kopplungsgrad:',k,'+-',errk)
        ks.append(k)
        errks.append(errk)
        
#Plot k gegen l
plt.figure()
plt.errorbar(l_F,ks,yerr=errks,xerr=errl_F,fmt='ko')
plt.xlabel('Abstand der Feder von der Pendelachse $l_F$ in cm')
plt.ylabel('Kopplungsgrad $k$')
plt.title('Abhängigkeit des Kopplungsgrades vom Drehachsenabstand')
plt.grid()
plt.savefig('Images/Gekoppelt_Kopplung')

ks=np.array(ks)
errx=np.array(errks/ks**2)
#erry=m*g*l_s/(l_F/100)**2*np.sqrt((errm/m)**2+(errg/g)**2+(errl_s/l_s)**2+(2*errl_F/(l_F/100))**2)
#TODO: Sind das wirklich statistische?
erry=m*g*l_s/(l_F/100)**2*2*errl_F/(l_F/100)

D,errD,b,errb,chiq,corr=analyse.lineare_regression_xy(1/ks-1,m*g*l_s/(l_F/100)**2,errx,erry)
print('Federkonstante:(',D,'+-',errD,')N/m, Chi²/f=',chiq/2)
#TODO: Chi² zu klein

#Regressionsplot
plt.figure()
plt.errorbar(1/ks-1,m*g*l_s/(l_F/100)**2,yerr=erry,xerr=errx,fmt='ko')
steps=np.arange(min(1/ks-1),max(1/ks-1),0.001)
plt.plot(steps,D*steps+b,color='red')
plt.ylabel('$m*g*l_s/l_F^2$ in N/m')
plt.xlabel('$1/k-1$')
plt.title('Regression zur Berechnung von D')
plt.grid()
plt.savefig('Images/Gekoppelt_Federkonstante')

#Residuum
plt.figure()
plt.errorbar(1/ks-1,m*g*l_s/(l_F/100)**2-(D*(1/ks-1)+b),yerr=np.sqrt(erry**2+errx**2),fmt='ko')
plt.ylabel('$\\frac{mgl_s}{l_F^2}$ in N/m')
plt.xlabel('$1/k-1$')
plt.title('Residuum zur Berechnung von D')
plt.grid()
plt.savefig('Images/Gekoppelt_Federkonstante_Residuum')