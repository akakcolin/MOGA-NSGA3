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
  
  *to the total number of epochs and total population size required for the training*
  
  For example: If the total number of epochs is 300 and population size is 200 then modify as follows:
  
  `int iteration_num = 300;   // number of training iterations`
  
  `int population_num = 200; // population size for ga.in training`
  
  *Once these modifications are done compile the code by running Makefile*
  
  `make`
  
## Step 3: Modify var.in to specify all the input parameters
The structure of var.in is as follows:
300
14
5
0.8
0.2
1.0
10
2.2     5.5     2.8     1.3     0.4     0.24    30.0    10.7    27.0    50.70   15.00   1.60    4.50    74.0
2.4     5.9     3.0     2.0     0.7     0.38    45.2    16.1    40.7    76.20   23.00   2.41    6.80    85.0
  
  
  
