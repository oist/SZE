% obtain le work code, stick to python script

%{
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
%}


%prob Array needs to be in form of 1 value per square, not psi!
function[wtot] = calcWork(posArray)
w = 0;
for i=1:length(posArray)
    for j=1:i
        pn = posArray(i,j);
        w = w -pn*log(pn);
    end
end
wtot = w;
end

