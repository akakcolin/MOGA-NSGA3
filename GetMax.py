import numpy as np

ffile = open("last400", "r")
lines = ffile.readlines()
ffile1 = open("dist.txt", "w")

def computedist(obj):
	dist = 0.0
	for i in range(len(obj)):
		dist += obj[i]*obj[i]
	return np.sqrt(dist)

obj = [0.0, 0.0, 0.0, 0.0, 0.0]

for line in lines:
	line = line.strip().split()
	pop = int(line[0])
	obj = [float(line[-1]), float(line[-2]), float(line[-3]), float(line[-4]), float(line[-5])]
	print(obj)
	dist = computedist(obj)
	print(dist)
	if pop < 500:
		ffile1.write("%d \t " %(pop))
		for i in range(1,15):
			ffile1.write("%12.6f " %(float(line[i])))
		ffile1.write("%12.6f \n" %(float(dist)))


