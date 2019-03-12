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

`300        # Size of the population (same as population_num in moga.c) `<br/>
`14         # Total number of variables in the force field file that need modification`<br/>
`5          # Total number of objectives` <br/>
`0.8        # Cross-Over Probability`<br/>
`0.2        # Mutation Probability`<br/>
`1.0        # Degree of Mutation between [0-1]`<br/>
`10         # Total number of divisions (will be read from hyperplane.in file, any value here will be ignored)`<br/>
`2.2     5.5     2.8     1.3     0.4     0.24    30.0    10.7    27.0    50.70   15.00   1.60    4.50    74.0 # Lower bound of variables`<br/>
`2.4     5.9     3.0     2.0     0.7     0.38    45.2    16.1    40.7    76.20   23.00   2.41    6.80    85.0 # Upper bound of variables`<br/>
  
 Please modify each line with suitable values consistent with the comments corresponding each line
  
  
