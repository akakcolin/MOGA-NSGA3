

all : ga moga
	mv ga moga ../.

ga : nsga3.c
	mpicc nsga3.c -lm -o ga

moga : moga.c
	mpicc moga.c -lm -o moga

clean:
	mv ../ga ../moga .
	rm ga moga
