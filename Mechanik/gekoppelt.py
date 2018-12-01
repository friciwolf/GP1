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
err_g = Erdbeschleunigung.errg_sys + Erdbeschleunigung.errg_stat #TODO: Richtig?

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
    print(len(t1))
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
    
        
    #figure()
    #gs = gridspec.GridSpec(file_count, 1)
    for j in range(0,file_count):
        figure()
        gs = gridspec.GridSpec(2, 1)
        ax1 = plt.subplot(gs[0, 0])
        plt.plot(t[j],U1[j],color='black')
        plt.title("Messung "+str(i)+": "+name[j])
        plt.ylabel(u"U1/V")
        plt.setp(ax1.get_xticklabels(), visible=False)
        grid()
        
        ax2 = plt.subplot(gs[1, 0],sharex = ax1)
        plt.plot(t[j],U2[j],color='red')
        plt.ylabel(u"U2/V")
        plt.xlabel(u"t/s")
        ax2.yaxis.tick_right()
        grid()
        plt.subplots_adjust(hspace=.0)
        
        plt.savefig("Images/Messung "+str(i)+"_"+name[j]+".jpg")
