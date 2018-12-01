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
name = ["schwebung","gleichsinnig","gegensinnig"]

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
        
        #Fourier analyse
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
