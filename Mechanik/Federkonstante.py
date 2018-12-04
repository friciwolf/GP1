# -*- coding: utf-8 -*-

import Erdbeschleunigung
plt.close('all')

#Wertetabelle und Fehler
m=np.array([0,59.9,109.9,130.1])*0.001
x=(54-np.array([54,42.5,33,29]))*0.01

errm=np.array(np.ones(len(m))*0.0001) #kg
errx=np.array(np.ones(len(m))*0.001) #m

g=Erdbeschleunigung.g
errg=np.sqrt(Erdbeschleunigung.errg_stat**2+Erdbeschleunigung.errg_sys**2) #nur noch systematisch
errD_sys=0.0 #0, da konstante Verschiebung von g? #TODO: 

#Lineare Regression
D,errD_stat,b,errb,chiq,corr=analyse.lineare_regression_xy(x,m*g,errx,errm*g) #Chi² zu klein
print('Federkonstante: D=({} +- {} +- {})N/m mit Chi²/f={}'.format(D,errD_stat,errD_sys,chiq/(len(m-2))))

#Plot
plt.plot(x,m*g,'ko',label='Messpunkte')
t=np.arange(0,0.3,0.001)
plt.plot(t,D*t+b,color='red',label='Regression D*x+b='+str(round(D,2))+'*x+'+str(round(b,2)))
plt.title("Messpunkte zur Ermittlung der Federkonstanten D aus dem Hooke'schen Gesetz")
plt.xlabel('Verlängerung der Feder in m')
plt.ylabel('Gewichtskraft in N')
plt.legend()
plt.grid()
plt.savefig('Images/Federkonstante_Regression.jpg')
plt.figure()

#Residuum
plt.errorbar(x,m*g-(D*x+b),yerr=np.sqrt((errm*g)**2+errx**2),capsize=0.5,fmt='ko')
plt.title('Residuenplot')
plt.xlabel('x in m')
plt.ylabel('m*g-(Dx+b)') #Gast mag große Bilder
plt.plot(x,0*np.ones(len(x)),'-.r')
plt.savefig('Images/Federkonstante_Residuum.jpg')
