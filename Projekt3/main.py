import numpy as np
import matplotlib.pyplot as plt
import math as mat

a = 6378137
e2 = 0.00669438002290

def Np(B):
    N = a/(1-e2*(np.sin(B)**2))**0.5
    return N

def hirvonen(X,Y,Z):
    r = (X**2 + Y**2)**0.5
    B = mat.atan(Z/(r*(1-e2)))
    
    while 1:
        N = Np(B)
        H = r/np.cos(B) - N
        Bst = B
        B = mat.atan(Z/(r*(1-(e2*(N/(N+H))))))    
        if abs(Bst-B)<(0.00001/206265):
            break
    L = mat.atan2(Y,X)
    N = Np(B)
    H = r/np.cos(B) - N
    return B, L, H 

#Wczytanie danych
wyniki = np.genfromtxt('pos_MIMA1.pos', comments='%')



xyzref=[3655333.847, 1403901.067,  5018038.047]
dxyz=wyniki[:,2:2+3]-xyzref

phi, lam, h= hirvonen(xyzref[0],xyzref[1],xyzref[2])
Rneu= np.array([[-np.sin(phi) * np.cos(lam), -np.sin(lam), np.cos(phi) * np.cos(lam)],
                        [-np.sin(phi) * np.sin(lam), np.cos(lam), np.cos(phi) * np.sin(lam)],
                        [np.cos(phi), 0, np.sin(phi)]])

dneu=[]
for dx in dxyz:
    dn=Rneu.T@dx
    dneu.append(dn)

dneu=np.array(dneu)


# Różnica między współrzędnymi wynikowymi i referencyjnymi w NEU
num_points = dneu.shape[0]
time_in_hours = np.arange(0, num_points * 0.5 / 60, 0.5 / 60)
fig, ax = plt.subplots(3, 1, figsize=(10, 8))
for i in range(3):
    ax[i].plot(time_in_hours, dneu[:, i])
    ax[i].set_xlabel('Czas (godziny)')
fig.suptitle('Różnice między współrzędnymi wynikowymi i referencyjnymi w NEU')
plt.tight_layout()
plt.show()


#Wykres testu ratio
ratio=wyniki[:,-1]
plt.figure(figsize=(10, 5))
plt.plot(time_in_hours, ratio,color='b')
plt.xlabel('Czas (godziny)')
plt.title('Wartości testu ratio')
plt.legend()
plt.grid(True)
plt.show()

fix=wyniki[:,5]==1
floats=wyniki[:,5]==2
fig, ax=plt.subplots(3,1)
t=np.arange(0,len(wyniki))
tfix=t[fix]
tfloat=t[floats]


for i in range(3):
    ax[i].plot(tfix,dneu[fix,i], color='b', label='fix')
    ax[i].plot(tfloat,dneu[floats,i], color='r', label='float')
fig.suptitle('Liczba poszczególnych typów rozwiązań wektora(fix, float)')
plt.legend(loc='upper left')
plt.show()

fig, ax= plt.subplots()
ax.scatter(dneu[fix,1],dneu[fix,0], color='r', s=5, label='Fixed')
ax.scatter(dneu[floats,1],dneu[floats,0], color='b', s=5, label='Float')
ax.axhline(0)
ax.axvline(0)
ax.axis('equal')
ax.set_title('Rozkład punktów')

plt.legend()
plt.show()


std=np.std(dneu[fix,:], axis=0)
rms=np.sqrt(np.sum(dneu[fix,:]**2, axis=0))
x=np.array([1,2,3])
neu = np.array(['N', 'E', 'U'])

# Wykres dla odchylenia standardowego
plt.figure(figsize=(8, 4))
plt.bar(x, std, width=0.5, color='b', align='center')
plt.title('Odchylenie standardowe')
plt.xticks(x, neu)
plt.show()

# Wykres dla RMS
plt.figure(figsize=(8, 4))
plt.bar(x, rms, width=0.5, color='r', align='center')
plt.title('RMS')
plt.xticks(x, neu)
plt.show()

numfix=len(dneu[fix,:])
numfloat=len(dneu[floats,:])

fig, ax= plt.subplots()
colors = ['blue', 'lightblue'] 
ax.pie([numfix, numfloat], labels=['Fix', 'Float'], autopct='%1.1f%%', colors=colors)
ax.set_title('Procentowy udział punktów Fix i Float')
plt.show()


max_bledy=np.max(np.abs(dneu), axis=0)
fig, ax= plt.subplots()
ax.bar(['N', 'E', 'U'], max_bledy)
plt.title('Wartości maksymalnych błędów')

plt.show()
