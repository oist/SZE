#SZE_2P
#determine work output of szilard engine using probability calculated
#from WF

#Note: L not stored in filename!

#to do
#plot w(nbar)
#output cleanup

import math
import numpy as np
import os, sys
import matplotlib.pyplot as plt

#Settings
scriptDir = "Documents/OIST/"
nbarList = [9,10,11,12,13,14,15,16]
gammaList = [] 					#if not empty there must be at least one entry for ALL nb values
if not gammaList:
	for i in range(len(nbarList)):
		gammaList.append([])	
for cb,nb in enumerate(nbarList):
	gammaList[cb].append((float(nb)-1)/(nb+1))
L = 10
print gammaList

#calcWFequal = True #shall the WF and probs for the last state in case particles are in different wells be calculated?

def outName(nb, gam, txtname):
	return scriptDir + "sze_out/" +txtname+"_nbar_" +str(nb)+ "_gamma_"+str(gam) + ".txt"

def generateWF(nb, gam):
	print("generating WF for nb: " + str(nb) + " gamma: " + str(gam))
	command = "wolframscript -file " + scriptDir + "mathematica/altered_code/RandomLatticeTG_ShortMoveOut_mod2.wl " + str(nb) + " " + str(gam) + " " + str(L) + " >" + outName(nb,gam,"WF2p")
	print(command)
	os.system(command)
	
 
def generateProb(nb, gam):
	print("generating probability arrays for nb: " + str(nb) + " gamma: " + str(gam))
	command = "wolframscript -file " + scriptDir + "mathematica/my_code/probList.wl " + outName(nb, gam, "WF2p") + " >" + outName(nb, gam, "probArray")
	print(command)
	os.system(command)

def generate(nbarL, gammaL):
	for cn, nbar in enumerate(nbarL):
		print("nbar " + str(nbar))
		for gamma in gammaL[cn]:
			print("gamma " + str(gamma))
			generateWF(nbar, gamma)
			generateProb(nbar, gamma)

#calculates QSE work output for one process with set nbar and gamma
def calcWork(pnArr, const):
	work = []
	for i, row in enumerate(pnArr):
		for j, pn in enumerate(row):
			if i ==j:
				#pout = 1
				red = 0
			else:
				red = const 
			#if i < j:
			#	pout = poutArr[0][1]
				
			#if i >j: 
			#	pout = poutArr[1][0]
			w = -pn*math.log(pn) + pn*red
			print 'i' + str(i) + ' j' +str(j)
			#print 'pout ' + str(pout)
			print 'pn '+ str(pn)
			print 'consRed '+ str(red)
			print 'w '+ str(w)
			work.append(w)
	wtot = sum(work)
	wrel = wtot / math.log(2) 
	print work
	print wtot
	print wrel
	return [wtot, wrel]


def convertMathToPythonArray(mathematicaArray):
	print mathematicaArray
	splitArr = mathematicaArray.split('},')
	print(splitArr)
	while '' in splitArr:
		splitArr.remove('')
	print(splitArr)
	for i in range(len(splitArr)):
		print splitArr[i]			
		splitArr[i] = splitArr[i].replace('{', '').replace('}', '').replace('\n', '').split(',')
	print splitArr
	for i in range(len(splitArr)):
		for j in range(len(splitArr[i])):
			splitArr[i][j] = float(splitArr[i][j].replace('*^', 'e'))
	print splitArr
	return splitArr

#convert matlab output to python arrays
def getValTab(nbarL, gammaL):
	#valTable = np.zeros((len(nbarL)*len(gammaL), 3))
	valTable = []
	for cb, nbar in enumerate(nbarL):
		for cg, gamma in enumerate(gammaL[cb]):
			matpnArr = extractArrayFromFile(outName(nbar, gamma, "probArray"))
			#pos = cb +cg
			valTable.append([nbar, gamma, matpnArr])
			#valTable[pos][0] = nbar
			#valTable[pos][1] = gamma
			#valTable[pos][2] = matpnArr
			
	return valTable
	
def extractArrayFromFile(filename):
	probFile =open(filename, "r")
	text = probFile.read()
	return convertMathToPythonArray(text.split("probList")[1])
	 
def partFunc(kList):
	beta = 1
	Esum = 0
	for k in kList:
		Esum = Esum + getE(k)
	return math.exp(-beta*Esum)
	
def getE(k):
	#E = (hb*k)**2/(2m)
	hb =1
	m =1
	return (hb*k)**2/(2*m)
#some hard coding first, generalize later!


def workingSteps():
	nbar = nbarList[0]
	gamma = gammaList[0]
	#wIns
	#Z(L), simple infininte well, k = npi/L
	ZL = partFunc([math.pi/L, 2*math.pi/L])
	sumZsigmaLins = a
	
	wIns = math.log(partFunc(L))
	#wExp
	wExp
	#wRem
	wRem

def main():
	#if calcWFequal:
	#	generate(nbarList, [0])
	#check if file exists...
	#poutArr = extractArrayFromFile(outName(nbarList[0],0, "probArray"))
	generate(nbarList, gammaList)
	kListL2 = [0.62800452871160091378254069275765869245	,  0.6283185307179586476925286766760615023, 
	 1.25600905897040747407039528562767678525, 1.2566370614359172953850573533521230046, 
	 1.88401359232357960112130950218395361334,   1.88495559215387594307758602996770173052, 
	2.51201813031818576609764169980682356048, 2.51327412287183459077011470670424600921, 
	 3.14002267450115727495157720559353193228, 3.14159265358979323846265123557554033663, 
	3.76802722641924256245126564242825732829, 3.76991118430775188615517205993540346104]
	ZL = partFunc([math.pi/L, 2*math.pi/L])
	ZL2 =partFunc([kListL2[0], kListL2[1]])
	#print ZL
	#print ZL2
	#print(math.log(ZL2/ZL))
	#print '------'
	EL = getE(math.pi/L) + getE(2*math.pi/L)
	EL2 = (getE(kListL2[0]) + getE(kListL2[1]))
	#print EL
	#print EL2 
	#print(- EL2+ EL)
	#SAVE VALUE TABLE
	constRed = -EL2+ EL
	print constRed
	valueTable = getValTab(nbarList, gammaList)
	#print poutArr
	print valueTable[0][-1]
	for cs, scenario in enumerate(valueTable):
		work = calcWork(valueTable[cs][-1], constRed)
		scenario.append(work)
	print valueTable
	wtotList = []
	wrelList = []
	for i in range(len(valueTable)):
		wtotList.append(valueTable[i][-1][0])
		wrelList.append(valueTable[i][-1][1])
	print wtotList
	print wrelList
	plt.plot(nbarList, wtotList, 'ro')
	plt.show()
	#calcWork(valueTable[0][-1], constRed)
	
def figureGamma():
	a = gam*L/(nb-1) #spacing between barriers
	y0 = -L/2+ (L-(nb-1)*a)/2 #first barrier position
	y0 = -L/2+ (L-(nb-1)*gam*L/(nb-1))/2 #first barrier position
	return 42
	
	

	      
main()
