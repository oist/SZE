#SZE_ENTROPYCALC

import math
import numpy as np
import matplotlib.pyplot as plt
print 'yes'
from scipy.integrate import nquad
from scipy.integrate import nquad


#wave function solution, infinite square well
def psi(x,n, l):
	return math.sin((x+L/2)*n*math.pi / l)
	
def psiL(x1, x2, l):
	return math.fabs(psi(x1, 1,  l)*psi(x2, 2, l)-psi(x2, 1, l)*psi(x1, 2, l))
	
def psiL2(x1, x2, l):
	return math.fabs(psi(x1, 1,  l)*psi(x2, 1, l)-psi(x2, 1, l)*psi(x1, 1, l))

#calculate entropy from probability array, includes indistinuishability

def calcEntr(pArr):
	entropyL = []
	newNorm = 0
	for i, row in enumerate(pArr):
			for j in range(i+1):
				pn =pArr[i][j]
				newNorm = newNorm + pn
				#print "i, j: " + str(i) + ' ' + str(j)
				#print ("pn " + str(pn))
				s = -pn*math.log(pn)
				#print 'entropy single ' + str(s)
				entropyL.append(s)	
	entropytot= sum(entropyL)
	#print newNorm
	#print entropytot
	entropytot = entropytot / newNorm			
	print entropytot
	return entropytot

def inInterval(num, intMin, intMax):
	return intMin <= num <= intMax

def getProbArrayAnalytic(l):
	pArray = np.zeros((3,3))
	norm = 0
	for i in range(nb+1):
		for j in range(nb+1):
			pv = nquad(psiL, [[barrierPos[i], barrierPos[i+1]], [barrierPos[j], barrierPos[j+1]]], args = [l])
			pArray[i][j]= pv[0]
			norm = norm + pv[0]
	pArr = pArray/norm
	return pArr
def getProbArrayAnalyticL2(l):
	pArray = np.zeros((3,3))
	norm = 0
	for i in range(nb+1):
		for j in range(nb+1):
			pv = nquad(psiL2, [[barrierPos[i], barrierPos[i+1]], [barrierPos[j], barrierPos[j+1]]], args = [l])
			pArray[i][j]= pv[0]
			norm = norm + pv[0]
	pArr = pArray/norm
	return pArr

#case: specify L
L = 10
#case: barrier positions for even spaced gamma and nbar = 2
nb = 2
barrierPos = [-1.666666666665, 1.666666666665, -5, 5]
barrierPos.sort()

#Infinite square well - L
xv = np.linspace(-L/2, L/2, 100)
yv1 = np.zeros(len(xv))
yv2 = np.zeros(len(xv))
for i in range(len(xv)):
	yv1[i] = psi(xv[i], 1, L)
	yv2[i] = psi(xv[i], 2, L)
#plot first two states
plt.plot(xv, yv1, 'r-', xv, yv2, 'g-')
for entry in barrierPos:
	plt.axvline(x=entry, alpha = 0.3, color = 'k', linestyle = '--')

#calculate prob of position states from wf (normalize to 1)
pArrL = getProbArrayAnalytic(L)
#pArrL2 = getProbArrayAnalyticL2(L/2)

#half infinite square well - L/2
xv = np.linspace(-L/2, L/2, 100)
yv1 = np.zeros(len(xv))
yv2 = np.zeros(len(xv))
for i in range(len(xv)):
	yv1[i] = psi(xv[i], 1, L/2)
	yv2[i] = psi(xv[i], 2, L/2)
#plot first two states
plt.plot(xv, yv1, 'r-', xv, yv2, 'g-')
for entry in barrierPos:
	plt.axvline(x=entry, alpha = 0.3, color = 'k', linestyle = '--')


#probability array for nb =2, gamma = 1/3
b2arr = [[6.562148641429094e-7, 0.16666574720514407, 0.16666692379332065],
[0.16666574720514407, 1.8449779245197027e-6, 0.16666574764966313],
[0.16666692379332065, 0.16666574764966313, 6.56214867280625e-7]]

bsplitArr = [[1.30974e-10,0.0786478,0.323589],[0.0786478,0.0382305,0.0786478],
				[0.323589,0.0786478,1.30974e-10]]
#entropy for L = 10
print 'L'
print pArrL
SL = calcEntr(pArrL)
print 'L/2'
print bsplitArr
SL2 = calcEntr(bsplitArr)
#entropy for 2 bar case
print 'nb2'
print b2arr
S2bar = calcEntr(b2arr)
#calculate entropy fpr particles in different gaps scenario
#extract prob array from nbar 2, gamma 0 

#new section for energies

def getE(k):
	#E = (hb*k)**2/(2m)
	hb =1
	m =1
	return (hb*k)**2/(2*m)
	
	
def Estate(k1,k2):
	return getE(k1) + getE(k2)
	
#energy state L
kL_1 = 1*math.pi/L
kL_2 = 2*math.pi/L
EL = Estate(kL_1, kL_2)
#energy state L2
kL2_1 = math.pi/(L/2)
kL2_2 = kL2_1
EL2 = Estate(kL2_1, kL2_2) 
#energy state nb2
knb2_1 = 0.940361999
knb2_2 = 0.941771471
Enb2 = Estate(knb2_1, knb2_2)

EArr = [EL, Enb2, Enb2, EL2, EL]
SArr = [SL, S2bar, 0, SL2, SL]
plt.clf()
plt.plot(EArr, SArr, 'g-o')
plt.xlabel('Energy')
plt.ylabel('Entropy')

labels = ['1 - L', '2 - nb2', ' 3 - nb2_meas', '4 - L2']


for label, x, y in zip(labels, EArr, SArr):
    plt.annotate(
        label,
        xy=(x, y), xytext=(-20, 20),
        textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

print '-'*20
print 'cycle so far'
print '-'*20
print 'insertion process'
print 'energy changes from ' + str(EL) + ' to ' + str(Enb2)
print 'entropy changes from ' + str(SL) + ' to ' + str(S2bar)
print 'Entropy change : ' + str(S2bar - SL)
print '-'*20
print 'measurement process'
print 'energy does not change ??'
print 'entropy changes from ' + str(S2bar) +' to ' + str(0)
print 'Entropy change : ' + str(0 -S2bar)
print '-'*20
print 'expansion process'
print 'energy changes from ' + str(Enb2) + ' to ' + str(EL2)
print 'entropy changes from ' + str(0) +' to ' + str(SL2)
print 'Entropy change : ' + str(SL2 - 0)
print '-'*20
print 'removal process'
print 'energy changes from ' + str(EL2) + ' to ' + str(EL)
print 'entropy changes from ' + str(SL2) +' to ' + str(SL)
print 'Entropy change : ' + str(SL - SL2)


plt.show()



