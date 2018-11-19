# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 20:16:31 2018

@author: Anwender
"""

import Widerstand


#Aufladung
tA, UA, IA = [], [], [] #Zeit, Spannung, Strom
dataA = cassy.CassyDaten("Kondensator_aufladen.lab")
tA = dataA.messung(1).datenreihe("t").werte[0:545]
UA = dataA.messung(1).datenreihe("U_B1").werte[0:545]
IA = dataA.messung(1).datenreihe("I_A1").werte[0:545]
                
#Entladung
tE, UE, IE = [], [], [] #Zeit, Spannung, Strom
dataE = cassy.CassyDaten("Kondensator_entladen.lab")
tE = dataE.messung(1).datenreihe("t").werte[0:545]
UE = dataE.messung(1).datenreihe("U_B1").werte[0:545]
IE = dataE.messung(1).datenreihe("I_A1").werte[0:545]

#Theo-Werte
Uoffset=min(UE)
Ioffset=min(IA)
wvm=np.argmax(IA==max(IA))+2 #weg vom Maximum bei IA, damit nur monoton fallend (siehe Rohdatenplot)

C_Hersteller=10**(-5) #Farad 
tau_theo=-1/(Widerstand.R*C_Hersteller)

UE=UE-Uoffset+0.0215 #0.0215
IA=IA-Ioffset+0.00032 #0.00025      #Hier rumschrauben


if __name__=='__main__':
    #Rohdaten-Plot
    plt.title('Aufladen des Kondensators')
    plt.xlabel('Zeit in s')
    plt.plot(tA,UA, color='darkblue', label='Spannung')
    plt.ylabel('Spannung in Volt',color='darkblue')
    plt.twinx()
    plt.plot(tA,IA, color='darkorange', label='Strom')
    plt.axvline(tE[wvm], color="red", linestyle = "--")
    plt.ylabel('Strom in Ampere', rotation=270,verticalalignment='bottom',color='darkorange')
    #plt.savefig('Images/Kondensator_Auf.jpg')
    plt.figure()
    
    plt.title('Entladen des Kondensators')
    plt.xlabel('Zeit in s')
    plt.plot(tE,UE, color='darkblue', label='Spannung')
    plt.ylabel('Spannung in Volt',color='darkblue')
    plt.twinx()
    plt.plot(tE,IE, color='darkorange', label='Strom')
    plt.ylabel('Strom in Ampere', rotation=270,verticalalignment='top',color='darkorange')
    #plt.savefig('Images/Kondensator_Ent.jpg')
    plt.figure()

'''
#Zwei Werte aus jeder Messung herausgreifen
N1=50
N2=100
#U_1A=UA[N1]
#U_2A=UA[N2]
U_1E=UE[N1]
U_2E=UE[N2]
I_1A=IA[N1]
I_2A=IA[N2]
I_1E=IE[N1]
I_2E=IE[N2]

#C_UA=(tA[N2]-tA[N1])/Widerstand.R/ana.log(U_1A/U_2A) #TODO: anders umzustellen
C_UE=(tA[N2]-tA[N1])/Widerstand.R/ana.log(U_1E/U_2E)
C_IA=(tA[N2]-tA[N1])/Widerstand.R/ana.log(I_1A/I_2A)
C_IE=(tA[N2]-tA[N1])/Widerstand.R/ana.log(I_1E/I_2E)
'''
#linus
#UA_log=ana.log(UA)
UE_log=ana.log(UE)
IA_log=ana.log(IA[wvm:])
tA=tA[wvm:]
#IE_log=ana.log(IE)

#Parameter Tau bestimmen: #Hier ist tau=-1/tau
tauE, etauE, bE, ebE, chiqE, corr=ana.lineare_regression(np.array(tE),np.array(UE_log),np.array(Widerstand.U_err_stat/UE))
tauA, etauA, bA, ebA, chiqA, corr=ana.lineare_regression(np.array(tA),np.array(IA_log),np.array(Widerstand.I_err_stat/IA[wvm:]))
x,tauE,etauE=Widerstand.round_good(0,tauE,etauE)
x,tauA,etauA=Widerstand.round_good(0,tauA,etauA)

if __name__=='__main__':
    
    #Plot Tau-Bestimmung Entladung
    plt.plot(tE,UE_log,'bo', label='logarrithmierte Messreihe U_Entladung')
    plt.plot(tE,tauE*tE+bE, color='brown', label=u'Geradenanpassung $\\tau$*x+b mit $\\tau$={}, b={} mit Chi²/f={}'.format(tauE,round(bE,2),round(chiqE/len(tE),2)))
    plt.title('Bestimmung des Parameters $\\tau$ aus der Entladungsspannung')
    plt.xlabel('Zeit in s')
    plt.ylabel('log(U_Entladung)')
    plt.legend()
    #plt.savefig('Images/Kondensator_TauE.jpg')
    plt.figure()
    
    #Plot Tau-Bestimmung Aufladung
    plt.plot(tA,IA_log,'bo', label='logarrithmierte Messreihe I_Aufladung')
    plt.plot(tA,tauA*tA+bA,color='brown', label=u'Geradenanpassung $\\tau$*x+b mit $\\tau$={}, b={} mit Chi²/f={}'.format(tauA,round(bA,2),round(chiqA/len(tA),2)))
    plt.title('Bestimmung des Parameters $\\tau$ aus dem Aufladungsstrom')
    plt.xlabel('Zeit in s')
    plt.ylabel('log(I_Aufladung)')
    plt.legend()
    #plt.savefig('Images/Kondensator_TauA.jpg')
    plt.figure()
    
#Gewichteter Mittelwert
tau_mittel,etau_mittel=ana.gewichtetes_mittel(np.array([tauA,tauE]),np.array([etauA,etauE]))

#Kapazität bestimmen
C_s=-1/(tau_mittel*Widerstand.R)
C_err_stat=C_s*(etau_mittel/tau_mittel)
C_err_sys=C_s*np.sqrt((Widerstand.errR_stat/Widerstand.R)**2+(Widerstand.errR_sys/Widerstand.R)**2) #TODO: So?
C_s,C_err_sys,C_err_stat=Widerstand.round_good(C_s,C_err_sys,-C_err_stat)

if __name__=='__main__':
    print('Entladung:  -1/tau = {} +- {} +- {} /s'.format(tauE,etauE,0.0))
    print('Aufladung:  -1/tau = {} +- {} +- {} /s'.format(tauA,etauA,0.0))
    print('Gesamt: C = {} +- {} +- {} Farad'.format(C_s,C_err_stat,C_err_sys))
    
    #Residuum
    s=5 #step, damit nicht ein Gewusel an Residuenpunkten
    
    plt.title('Residuenplot - Abweichung von der Geraden ln(U_Ent)-f(t)')
    plt.errorbar(tE[::s],(UE_log-(tauE*tE+bE))[::s],yerr=np.array(Widerstand.U_err_stat/UE)[::s],fmt='ko',capsize=2, markersize=3)
    plt.xlabel('t in s')
    plt.ylabel('ln(U)-t/$\\tau$+b in V')
    #plt.savefig('Images/Kondensatur_ResE')
    plt.figure()
    
    plt.title('Residuenplot - Abweichung von der Geraden ln(I_Auf)-f(t)')
    plt.errorbar(tA[::s],(IA_log-(tauA*tA+bA))[::s],yerr=np.array(Widerstand.I_err_stat/IA[wvm:])[::s],fmt='ko',capsize=2, markersize=3)
    plt.xlabel('t in s')
    plt.ylabel('ln(I)-t/$\\tau$+b in V')
    #plt.savefig('Images/Kondensatur_ResA')
    
