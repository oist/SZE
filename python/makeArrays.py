

gA = [-0.1, 0, 0.1]
nbA = [1,2,3,4]
TA = [0.1]

gL = []
nbL = []
TL = []
for g in gA:
	for nb in nbA:
		for T in TA:
				gL.append(g)
				nbL.append(nb)
				TL.append(T)
gs = str(gL)
nbs = str(nbL)			
Ts = str(TL)

gs = gs.replace(',', '')
nbs = nbs.replace(',', '')
Ts =Ts.replace(',', '')

print gs
print nbs
print Ts
