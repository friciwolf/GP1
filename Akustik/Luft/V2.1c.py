# -*- coding: utf-8 -*-
"""
Created on Sun Dec 30 12:06:59 2018

@author: Anwender
"""

import Kalibrierung
import scipy.optimize as opt


def sort(eins,zwei='none',drei='none'):
    '''
    Sortiert ein Array oder eine Liste aufsteigend
    Optional: Sortiert ein zweites Array oder eine zweite Liste so, dass zusammengehörige Werte
    noch an gleichen Stellen stehen
    return: array, (array)
    '''
    if len(eins)!=len(zwei) or len(eins)!=len(drei):
        raise ValueError('Arrays must have the same lenghth')
        
    uno=[]
    dos=[]
    tres=[]
    eins=list(eins)
    zwei=list(zwei)
    drei=list(drei)
    
    while eins!=[]:
        index=eins.index(min(eins))
        uno.append(eins.pop(index))
        if zwei!='none':
            dos.append(zwei.pop(index))
        if drei!='none':
            tres.append(drei.pop(index))        
    
    if zwei!='none' and drei!='none':
        return (np.array(uno),np.array(dos),np.array(tres))
    elif zwei!='none':
        return (np.array(uno),np.array(dos))    
    elif drei!='none':
        return (np.array(uno),np.array(tres))
    else:
        return np.array(uno)


L_Rohr=42.1 #cm
L_Rohr_estat=0.05 #cm
L_Rohr_esys=0.0003+0.0002*1 #Güteklasse II:a+b*L mit a=0.3mm, b=0.2mm/m, L= auf ganzen Meter gerundete zu messende Länge
f=1607.4 #????! #Hz
f_err=0.1 #Hz


#Daten einlesen
data = cassy.CassyDaten("V2.1C.lab")
n = np.array(data.messung(1).datenreihe("n").werte) #Nummer des Messpunktes, wobei n=1 der weiteste Abstand war
R = np.array(data.messung(1).datenreihe("R_B1").werte) #Widerstand in kOhm
U = np.array(data.messung(1).datenreihe('U_A1').werte) #Spannung in Volt

L,L_estat,L_esys=Kalibrierung.Länge(R,Kalibrierung.R_err[0]*len(R)) #in cm
L,U,L_esys=sort(L,U,L_esys) #s.o.

#Versuch eines Fits (gescheitert)
def func(x,k):
    return 2.3*abs(np.sin(k*(x-L[0])))

par,par_err=opt.curve_fit(func,L,U,p0=0.4)

#Wellenlänge
mitte=int(len(L)/2)
min1=np.argmax(U==min(U[:mitte]))
min2=np.argmax(U==min(U[mitte:]))

lamb=2*(L[min2]-L[min1]) #cm
lamb_err=0.5 #cm #großzügig abgeschätzt
lamb_esys=(L_esys[min1]+L_esys[min2])/2 #TODO: So?

#Rohdatenplot
plt.title('Stehende Welle im inneren des Rohrs bei {} Hz'.format(f))
plt.plot(L,U,'b-') #'bo'?
plt.xlabel('Länge in cm') #unbekannter Nullpunkt?
plt.ylabel('Amplitude in V')
#plt.plot(L,func(L,*par),'r--')
plt.axvline(L[min1],linestyle='--',color='lightgreen')
plt.text(L[min1+1],1.3,'Knoten bei {:.2f} cm'.format(round(L[min1],2)))
plt.axvline(L[min2],linestyle='--',color='lightgreen')
plt.text(L[min2+1],1.3,'Knoten bei {:.2f} cm'.format(round(L[min2],2)))
plt.savefig('Images/VC_Roh')

#Geschwindigkeit
v=lamb*f
v_estat=v*np.sqrt((lamb_err/lamb)**2+(f_err/f)**2)
v_esys=f*lamb_esys
v,v_esys,v_estat=Temperatur.round_good(v,v_esys,v_estat)
print('\nv_Schall=({} +- {} +- {}) cm/s'.format(v,v_estat,v_esys))

#Theoretische Geschwindigkeit
T=Temperatur.Temperatur(157)
print('Theoretische Geschwindigkeit bei T={}°C: {}m/s'.format(T,331.6+0.6*T))

