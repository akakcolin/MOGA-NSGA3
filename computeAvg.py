import numpy as np

ffile = open("ZerothFront.txt", "r")

lines = ffile.readlines()

MeanValues = []
StdValues = []

counter = 0
var = []
sum1 = np.zeros(5)

file1 = open("avgErr.txt", "w")
k = 0
for line in lines:
	line = line.strip().split()
	if len(line) == 0:
		counter+=1
		#print(counter)
		var = np.matrix(var)
		#print(var.shape)
		#print(np.mean(var,axis=0))
		MeanValues.append(np.array(np.mean(var,axis=0)[0]))
		StdValues.append(np.array(np.std(var,axis=0)[0]))
		var = []
	elif len(line) > 15:
		l = []
		c = 0
		for elt in line[1:]:
			c+=1
			#if c == 15: l.append(float(elt)/100.0)	
			l.append(float(elt))
		#print(l)
		var.append(l)

	if counter != 1 and counter%50 == 0 and len(line) == 0:
		k += 1
		sum1 = sum1/50.0
		file1.write("%d \t %12.6f \t %12.6f \t %12.6f \t %12.6f \t %12.6f\n" %(k,sum1[0], sum1[1], sum1[2], sum1[3], sum1[4]))
		sum1 = np.zeros(5)	

	elif len(line) == 0:

		sum1[0] += MeanValues[-1][0][14]
		sum1[1] +=  MeanValues[-1][0][15]
		sum1[2] +=  MeanValues[-1][0][16]
		sum1[3] +=  MeanValues[-1][0][17]
		sum1[4] +=  MeanValues[-1][0][18]		

#print(np.matrix(MeanValues).shape)
for i in range(counter):
	print("****** Iteration %04d" %i)
	print("%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f       %10.3f%10.3f%10.3f%10.3f%10.3f" %(MeanValues[i][0][0], MeanValues[i][0][1], MeanValues[i][0][2], MeanValues[i][0][3], MeanValues[i][0][4], MeanValues[i][0][5], MeanValues[i][0][6], MeanValues[i][0][7], MeanValues[i][0][8], MeanValues[i][0][9], MeanValues[i][0][10], MeanValues[i][0][11], MeanValues[i][0][12], MeanValues[i][0][13], MeanValues[i][0][14], MeanValues[i][0][15], MeanValues[i][0][16], MeanValues[i][0][17], MeanValues[i][0][18]))
	print("%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f       %10.3f%10.3f%10.3f%10.3f%10.3f" %(StdValues[i][0][0],StdValues[i][0][1], StdValues[i][0][2],StdValues[i][0][3],StdValues[i][0][4],StdValues[i][0][5],StdValues[i][0][6],StdValues[i][0][7],StdValues[i][0][8],StdValues[i][0][9],StdValues[i][0][10],StdValues[i][0][11],StdValues[i][0][12],StdValues[i][0][13],StdValues[i][0][14],StdValues[i][0][15],StdValues[i][0][16],StdValues[i][0][17],StdValues[i][0][18]))
