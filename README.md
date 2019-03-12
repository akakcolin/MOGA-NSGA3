# MOGA-NSGA3

To run the workflow follow the steps as mentioned below:

## Step 1: Provide data for running the force field fitting
1. Put reference band structure data in UTIL folder as done in this repo.
2. Put cell data corresponding to different states as cell1, cell2, cell3 
3. Put the template forcefield which you wish to train


## Step 2: Compile the nsga3 and worflow code:
1. Go into src/ directory: `cd src/`
2. Modify the following lines in moga.c:
  `int iteration_num = 500;   // number of training iterations`
  
  `int population_num = 300; // population size for ga.in training`
