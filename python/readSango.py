#readSango
import glob
import numpy as np
import math
import matplotlib.pyplot as plt

fileDir = 'sango/thermalSZE_Int/'
run = 'run2/'
print fileDir + run + '*.out'
filelist = glob.glob(fileDir + run + '*.out*')
print filelist

dataTable = []

for c,f in enumerate(filelist):
	#open
	#read
	thisFile =open(f, "r")
	text = thisFile.read()
	entries = text.split('\n')
	values = entries[10:-3]
	thisV = []
	for v in range(0,4):
		val = values[v].split('_')[1]
		thisV.append(val)
	dataTable.append(thisV)

print dataTable

#sort

Tvals = [0.1, 1, 10 ]
gvals = [-1, 0, 1]
nbvals = [1,2,3,4]

gpos = 0
Tpos = 1
nbpos = 2
wpos = 3

constnb = 4
constT = 0.1

cTList = []
cnbList= []

Escale = 0.0493

for ar in dataTable:
	if float(ar[Tpos]) == float(constT):
		cTList.append(ar)
	if float(ar[nbpos]) == float(constnb):
		cnbList.append(ar)
			
		
print cTList
print '--------'
print cnbList
print '--------'
gTList = np.zeros([len(gvals), len(Tvals)])
gnbList = np.zeros([len(gvals), len(nbvals)])


for cg, g in enumerate(gvals):
	for ar in cnbList:
		for ct,t in enumerate(Tvals):
			if float(ar[gpos]) == float(g) and float(ar[Tpos]) == float(t):
				gTList[cg][ct] = float(ar[wpos])/t/math.log(constnb+1)/Escale
	for ar in cTList:
		for cnb,nb in enumerate(nbvals):
			if float(ar[gpos]) == float(g) and float(ar[nbpos]) == float(nb):
				gnbList[cg][cnb] = float(ar[wpos])/constT/Escale
print gTList
print gnbList

plt.figure(1)
for cv,vals in enumerate(gTList):
	plt.plot(Tvals, vals, label=str(gvals[cv]))
plt.xlabel('temperature')
plt.ylabel('W/T/log(nb+1)/Escale')
plt.title('WT, g, constant nb ' + str(constnb))
plt.xscale('log')
plt.legend()

plt.figure(2)
for cv,vals in enumerate(gnbList):
	plt.plot(nbvals, vals, label=str(gvals[cv]))
plt.xlabel('number of barriers')
plt.ylabel('W/T/Escale')
plt.title('Wnb, g, constant T ' + str(constT))
plt.legend()
plt.show()
	

	
