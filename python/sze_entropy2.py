import numpy.linalg as lin
#import numpy.matrix as mx
import numpy as np
import math

b2arr = [[6.562148641429094e-7, 0.16666574720514407, 0.16666692379332065],
[0.16666574720514407, 1.8449779245197027e-6, 0.16666574764966313],
[0.16666692379332065, 0.16666574764966313, 6.56214867280625e-7]]

#some manual stuff due to lack of brainpower
psame = 0
p12 = b2arr[0][1]
p13 = b2arr[0][2]
p23 = b2arr[1][2]
thisTr = 2*(p12+p13+p23)
print thisTr
s_v1 = 2*(p12*math.log(2*p12)+p13*math.log(2*p13) + p23*math.log(2*p23))
print s_v1
print 'other'
 
thisT = p12+p13+p23
p12 = p12/thisT
p13 = p13/thisT
p23 = p23/thisT
print p12+p13+p23
s_v2 = p12*math.log(p12)+p13*math.log(p13)+p23*math.log(p23)
print s_v2
#sameresult!!! yayyy
exit()


zeroCutOff = 1e-5
for row in b2arr:
	for i in range(0, len(row)):
		if row[i] < zeroCutOff:
			row[i] =0
print b2arr

#normalize 
thisNorm = np.trace(b2arr)
ro = b2arr/thisNorm
print lin.norm(ro)

[ero, vro] = lin.eig(ro)

print ero
