#SZE_2P_MOD
#determine work output of szilard engine using probability calculated
#from WF

#Note: L not stored in filename!
#to do
#why crash at nb = 12
#check formulas
#what about working precision?
#output cleanup
print 'starting sze_2p_mod.py'

import math
import numpy as np
import os, sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from decimal import *
getcontext().prec = 32

print 'done importing'

def gammaEven(nb):
	return (Decimal(nb-1)/Decimal(nb+1))

#Settings
scriptDir = "Documents/OIST/"
nbarList = [7]
dList = [5.523] 	

#gammaList = [[Decimal(1)/Decimal(3)]]
gammaList = [[]]
allGammaEntry = False
evenGamma = True
if allGammaEntry:
	for entry in nbarList:
		gammaList.append([0.1,0.2,0.3,0.4,0.5, 0.6,0.7,0.8,0.9])		#enter elements that shall be available for all nbar
		#gammaList.append([0.01, 0.02, 0.03, 0.04, 0.045, 0.05, 0.08])							#enter elements that shall be available for all nbar
#gammaList = [[0.1,0.2,0.3,0.4,0.5, 0.6,0.7,0.8,0.9]]*len(nbarList) 		#if not empty there must be at least one entry for ALL nb values
#nbarList = [2]


if not gammaList:
	for i in range(len(nbarList)):
		gammaList.append([])	
if evenGamma:
	for cb,nb in enumerate(nbarList):
		thisg = gammaEven(nb)
		if not thisg in gammaList[cb]:
			gammaList[cb].append(thisg)
			gammaList[cb].sort()
			
pltstyle = ['r-o', 'g-s','b-*', 'c-^', 'm-v', 'k-8']		#min length: length nbarList
pairpltstyle = [['r-o','r-^'],['b-o','b-^'],['g-o', 'g-s'],['c-o', 'c-^'],['m-o', 'm-^'],['k-o', 'k-^']]

L = 10
kListL2 = [0.62800452871160091378254069275765869245	,  0.6283185307179586476925286766760615023, 
	 1.25600905897040747407039528562767678525, 1.2566370614359172953850573533521230046, 
	 1.88401359232357960112130950218395361334,   1.88495559215387594307758602996770173052, 
	2.51201813031818576609764169980682356048, 2.51327412287183459077011470670424600921, 
	 3.14002267450115727495157720559353193228, 3.14159265358979323846265123557554033663, 
	3.76802722641924256245126564242825732829, 3.76991118430775188615517205993540346104]
	


def outName(nb, gam, d, txtname):
	return scriptDir + "sze_out/" +txtname+"_nbar_" +str(nb)+ "_gamma_"+str(gam) + "_height_"+str(d) + ".txt"
def outFigName(nb, gam, figname):
	return scriptDir + "sze_out/" +figname+"_nbar_" +str(nb)+ "_gamma_"+str(gam) + "_height_"+str(d) + ".png"

#beware: script not updated to deal with varyig height, only for TG limit!
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


def generateWFTG(nb, gam, d):
	print("generating WF TG for nb: " + str(nb) + " gamma: " + str(gam) + " bar height: " + str(d))
	command = "wolframscript -file " + scriptDir + "mathematica/altered_code/RandomLatticeTG_ShortMoveOut_mod2_realTG_d.wl "\
				+ str(nb) + " " + str(gam) + " " + str(L) +  " " + str(d) + " >" + outName(nb,round(gam,2),d,"WF2p_TG_d")
	print(command)
	os.system(command)
	
def generateProbTG(nb, gam, d):
	print("generating TG probability arrays for nb: " + str(nb) + " gamma: " + str(gam) + " bar height: " + str(d))
	command = "wolframscript -file " + scriptDir + "mathematica/my_code/probList_d.wl " + outName(nb, round(gam,2),d, "WF2p_TG_d") \
			+ " >" + outName(nb, round(gam,2),d, "probArray_TG_d")
	print(command)
	os.system(command)

def generateTG(nbarL, gammaL, dList):
	for d in dList:
		for cn, nbar in enumerate(nbarL):
			print("nbar " + str(nbar))
			for gamma in gammaL[cn]:
				print("gamma " + str(gamma))
				generateWFTG(nbar, gamma,d)
				#generateProbTG(nbar, gamma, d)

#calculates QSE work output for one process with set nbar and gamma
def calcWork(pnArr, const):
	work = []
	newNorm = float(0)
	print pnArr
	#calculate the new probability norm taking indistinguishability into account
	for i, row in enumerate(pnArr):
		for j in range(i+1):
			pn =pnArr[i][j]
			newNorm = newNorm + pn
			if i ==j:
				red = 1
			else:
				red = const 
			print "i, j: " + str(i) + ' ' + str(j)
			print ("pn " + str(pn))
			w = -pn*math.log(pn) + pn*math.log(red)
			print 'work single ' + str(w)
			work.append(w)	
		
	print newNorm
	print work
	print sum(work)
	wtot = sum(work)/newNorm
	print wtot
	return wtot


def calcWorkDisting(pnArr, const):
	work = []
	for i, row in enumerate(pnArr):
		for j in range(len(row)):
			pn =pnArr[i][j]
			if i ==j:
				red = 1
			else:
				red = const 
			print "i, j: " + str(i) + ' ' + str(j)
			print ("pn " + str(pn))
			w = -pn*math.log(pn) + pn*math.log(red)
			print 'work single ' + str(w)
			work.append(w)	
	print work
	print sum(work)
	wtot = sum(work)
	return wtot



def convertMathToPythonArray(mathematicaArray):
	print mathematicaArray
	splitArr = mathematicaArray.split('},')
	while '' in splitArr:
		splitArr.remove('')
	for i in range(len(splitArr)):			
		splitArr[i] = splitArr[i].replace('{', '').replace('}', '').replace('\n', '').split(',')
	for i in range(len(splitArr)):
		for j in range(len(splitArr[i])):
			splitArr[i][j] = float(splitArr[i][j].replace('*^', 'e'))
	return splitArr

#convert matlab output to python arrays
def getValTab(nbarL, gammaL):
	valTable = []
	for cb, nbar in enumerate(nbarL):
		for cg, gamma in enumerate(gammaL[cb]):
			print 'nb, gamma '+ str(nbar) + ' ' + str(gamma)
			matpnArr = extractArrayFromFile(outName(nbar, gamma, "probArray"))
			valTable.append([nbar, gamma, matpnArr])			
	return valTable
	
def getValTabTG(nbarL, gammaL):
	valTable = []
	for cb, nbar in enumerate(nbarL):
		for cg, gamma in enumerate(gammaL[cb]):
			print 'nb, gamma '+ str(nbar) + ' ' + str(gamma)
			matpnArr = extractArrayFromFile(outName(nbar, gamma, "probArray_TG"))
			valTable.append([nbar, gamma, matpnArr])			
	return valTable
	
#extract mathematica array from txt file and return in python format
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


def plotGamma(valTab):
	a = 0.7
	#horizontal line for each even gamma
	for cb, nb in enumerate(nbarList):
		#my quick and dirty solution
		firstIndex = 0 
		for i in range(cb):
			firstIndex = firstIndex + len(gammaList[i])
		print firstIndex
		ng = gammaList[cb]
		w = []
		for i in range(len(ng)):
			w.append(valTab[firstIndex + i][-1])
		plt.plot(ng, w, pltstyle[cb],alpha = a, label = str(nb))
		plt.axvline(x=gammaEven(nb), alpha = 0.3, color = pltstyle[cb][0], linestyle = '--')
	plt.legend(title = 'nbar')
	plt.xlabel('gamma')
	plt.ylabel('Work output')
	plt.savefig(outFigName([gammaList[0][0], gammaList[0][-1]], 'gammaRange', 'new_workPlot_g_TG'))
	plt.show()
	

def plotNbar(valueTable, extraPlot):
	wtotList = []
	wrelList = []
	gamList = [] 
	nbL = []
	#plot for gamma - even, diff n
	#plot for gamma changing, same n
	#3d plot gamma and nbar
	wclassList = []
	for i in range(len(valueTable)):
		wtotList.append(valueTable[i][-1])
		wrel = math.log(valueTable[i][0]+1)
		wrelList.append(valueTable[i][-1]/wrel)
		wclassList.append(wrel)
		#gamList.append(valueTable[i][1])
		#nbL.append(valueTable[i][0])
		
	#for cb, nb in enumerate(nbarList):
	#	nbL.append([nb])
	#	for cv, thisList in enumerate(valueTable):
	print wtotList
	print wrelList
	
	xs = np.linspace(nbarList[0],nbarList[-1],20)
	#####horizontal line
	horiz_line_data = np.array([1 for i in xrange(len(xs))])
	if extraPlot:
		plt.plot(nbarList, wtotList, 'ro',label = 'W_q_indist')
		#plt.plot(nbarList, extraPlot, 'r^', label = 'W_q_disting')
		plt.plot(nbarList, wclassList, 'bo', label = 'W_class')
		plt.plot(nbarList, wrelList, 'go',label = 'W_q_indist/ W_class')
		plt.plot(xs, horiz_line_data, 'b--')
		plt.xlabel('nbar')
		plt.ylabel('Work output TG')
		
		
	else:	
		workPlot = plt.plot(nbarList, wtotList, 'ro', nbarList, wrelList, 'go',xs, horiz_line_data, 'b--',
				nbarList, wclassList, 'bo')
	#plt.show()
	plt.legend()
	plt.savefig(outFigName([nbarList[0], nbarList[-1]], 'even', 'newWorkPlot_TG'))
	plt.show()


	

def plotKgamma(nbList, gammaL):
	for cn, nbar in enumerate(nbList):
		band1 = []
		band2 =[]	
		for gamma in gammaL[cn]:
	#read WFplot file for gamma and nbar
			kArray = extractKarray(nbar, gamma)
			band1.append(kArray[0])
			band2.append(kArray[1])
		gammaL[cn].sort()
		a = 0.5
		plt.axvline(x=gammaEven(nbar), color = pairpltstyle[cn][0][0],alpha = 0.2,  linestyle = '--')
		plt.plot(gammaL[cn], band1, pairpltstyle[nbar-2][0], alpha = a, label= str(nbar) + '_k1')
		plt.plot(gammaL[cn], band2,pairpltstyle[nbar-2][1], alpha= a, label = str(nbar) + '_k2')
		plt.xlabel('gamma')
		plt.ylabel('k')
		plt.legend()
	plt.savefig(outFigName([nbList[0], nbList[-1]], [gammaList[0][0], gammaList[0][-1]], 'kGamma'))
	plt.show()
	#obtain first two k values an save to array
	#plot array vs nbar


def plotKnbar(nbList, gammaL):
	band1 = []
	band2 =[]
	for cn, nbar in enumerate(nbList):
		for gamma in gammaL[cn]:
			print 'nb, gamma' + str(nbar) + ' ' + str(gamma)
	#read WFplot file for gamma and nbar
			kArray = extractKarray(nbar, gamma)
			band1.append(kArray[0])
			band2.append(kArray[1])
	plt.plot(nbarList, band1, 'ro', label= 'band1')
	plt.plot( nbarList, band2,'g^', label = 'band2')
	plt.show()
	#obtain first two k values an save to array
	#plot array vs nbar
	
	
def extractKarray(nbar, gamma):
	probFile =open(outName(nbar, gamma, "probArray"), "r")
	text = probFile.read()
	snippet1 = text.split('{{')
	snippet2 = snippet1[1].split('}}')[0]
	mathematicaArray = snippet2
	splitArr = mathematicaArray.split('},')
	while '' in splitArr:
		splitArr.remove('')
	for i in range(len(splitArr)):			
		splitArr[i] = splitArr[i].replace('{', '').replace('}', '').replace('\n', '').split(',')
		splitArr[i] = splitArr[i][0].replace('k ->', '').replace("`20.",'')
		splitArr[i] = float(splitArr[i])
	return splitArr
	
def main():
	print 'starting main'
	print 'nbarList:'
	print nbarList
	print 'gammaList:'
	print gammaList
	print 'dList:'
	print dList
	#actual code
	
	#plotKgamma(nbarList, gammaList)
	#plotKnbar(nbarList, gammaList)
	#exit()
	generateTG(nbarList, gammaList, dList)
	exit()
	#print gammaList
	valueTable = getValTabTG(nbarList, gammaList)
	#print valueTable
	wDistingList = []
	for cs, scenario in enumerate(valueTable):
		work = calcWork(valueTable[cs][-1], float(1)/3) #change to 1/4 for distinguishable limit
		wDistingList.append(calcWorkDisting(valueTable[cs][-1], float(1)/4))
		scenario.append(work)
	print valueTable[0][-1]
	print valueTable[1][-1]
	print math.log(2)
	
	#plot: x- gamma, one plot for each n
	plotGamma(valueTable)
	#plotNbar(valueTable, wDistingList) #out name set to TG!
	exit()
	plotKnbar(nbarList, gammaList)
	plotDiagNbar(nbarList, gammaList)
	plotNentries(nbarList, gammaList)
	      
main()
