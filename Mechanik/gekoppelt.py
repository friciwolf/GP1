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
err_g = np.sqrt((Erdbeschleunigung.errg_sys)**2 + (Erdbeschleunigung.errg_stat)**2) #TODO: Richtig?

bis = 8000
name = ["Schwebung","Gleichsinnig","Gegensinnig"]
w_plus=[[0,0,0],[0,0,0],[0,0],[0]];w_minus=[[0,0,0],[0,0,0],[0,0],[0]]
ks=[];errks=[]
errw=[[0,0],0,0,0]
l_F=[1/0.165**2,1/(0.165+0.12)**2,1/(0.165+0.12+0.124)**2,1/(0.165+0.12+0.124+0.125)**2]

T=[[0,0,0],[0,0,0],[0,0],[0]]

for i in range(1,5):
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
        plt.savefig("Images/Position "+str(i)+"_"+name[j]+"Fourier.jpg")
        
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
    if i>1:
        errw[i-1]=(max(w)-min(w))/(len(w)-1)
        k=(w_minus[i-1][0]**2-w_plus[i-1][0]**2)/(w_minus[i-1][0]**2+w_plus[i-1][0]**2)
        errk=(2*(w_minus[i-1][0]/(w_minus[i-1][0]**2+w_plus[i-1][0]**2)-w_minus[i-1][0]*(w_minus[i-1][0]**2-w_plus[i-1][0]**2)/(w_minus[i-1][0]**2+w_plus[i-1][0]**2)**2)*errw[i-1])**2
        errk+=(-2*(w_plus[i-1][0]/(w_minus[i-1][0]**2+w_plus[i-1][0]**2)-w_plus[i-1][0]*(w_minus[i-1][0]**2-w_plus[i-1][0]**2)/(w_minus[i-1][0]**2+w_plus[i-1][0]**2)**2)*errw[i-1])**2
        errk=np.sqrt(errk)
        x,k,errk=Rauschmessung.round_good(0,k,errk)
        print(str(i),'Kopplungsgrad:',k,'+-',errk)
        ks.append(k)
        errks.append(errk)
        
    
#plt.close('all')
plt.figure()
ks = np.array(ks)
errks = np.array(errks)
ks = ks[1:]
errks = errks [1:]
ks = 1/ks
errys = errks*ks**2
l_F = np.array(l_F)
l_F = l_F[1:]
plt.errorbar(l_F,ks-1,yerr = errys, fmt = '.')
plt.ylabel('1/k')
plt.xlabel('1/l_F**2 in 1/m**2')
a,ea,b,eb,chiq,corr = analyse.lineare_regression(l_F,ks-1,errks)
plt.plot(l_F,a*l_F+b)
print("a: ",a)
print("ea: ",ea)
print("chiq: ",chiq)
