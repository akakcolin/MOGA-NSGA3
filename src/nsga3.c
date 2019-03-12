#include <stdio.h>
#include<stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <string.h>

int N = 500; // size of  population
int nvar=14; // number of variables to be manipulated
int nobj = 5; // number of objectives 
double crossover_P=0.5; // crossover percentage of population
double mutation_P=0.5; // mutation percentage of population
double mutation_rate = 1; // rate of mutation
int nDivision=10; // number of divisions for generating hyperplane
int size_hyperplane = 84;
int it_num=-1;
int suppress_output = 1;
struct bound {
    double x, y; // for reading in lower and upper bounds
};

/* functions to generate */
void create_hyperplane();
int partition(int f_p[][nobj],int index);
void permute(FILE *output,double sing_part[],int i);
void remove_duplicates(FILE *output, double hyperplane[][nobj],int size_hyperplane);

/* function to read in hyperplane */
void read_hyperplane(FILE *input,double hyperplane[][nobj],int size);

/* function to read in program's input file */
int read_var_input(char filename[],struct bound bounds[]);

/* function for reading population input files */
int read_pop_input(char fileName[],double x[][N],double y[][N]); 

/* nondominated sort for child + parent population and selection of fronts >= N */
int non_dom_sort(int pop,double y[][pop],int rank[pop],int dom_by[pop][pop]); 

/* selection of N population from 2N combined population */
int select_pop(int rank[2*N],int new_pop[N],int dom_by[2*N][2*N], int front[2*N]); // fills population as close to N by rank
int normalize(double obj[nobj][2*N], double norm_obj[nobj][2*N]); // normalization of rank for niche preservation
int associate(double norm_obj[nobj][2*N],double hyperplane[nobj][size_hyperplane],int nearest_pointer[2*N],double nearest_distance[2*N]); // association of normalized points to hyperplane
int niching(int pop_filled,int front[2*N],int nearest_pointer[2*N], int new_pop_point[N], double nearest_distance[2*N]);

/* creates output files */
int write_pop_file(char fileName[],char allName[],double x[][N], double y[][N], int type);

/* roulette selection for reproduction */
int select_crossover(int r[],double var_parent[nvar][N],double var_child[nvar][N], struct bound lower_upper[]);

/* mutant generation for remaining population */ 
int select_mutants(double var_parent[nvar][N], double var_child[nvar][N],struct bound lower_upper[]);

/* Matrix math functions */
double determinant(double [][nobj], double);
void cofactor(double [][nobj], double, double inv[nobj][nobj]);
void transpose(double [][nobj], double [][nobj], double, double inverse[nobj][nobj]);

int main(int argc, char *argv[]) {
    int iterationNum;
    sscanf(argv[1], "%d", iterationNum);
    int i,j;
    char old_input[] = "ga_old.in"; // stores previous generation
    char all_inputs[] = "ga_all.in"; // stores all generations
    char new_input[] = "ga.in"; // stores current generation
    char var_input[] = "var.in"; // stores variables
    char hyper[] = "hyperplane.in"; // stores hyperplane
    struct bound lower_upper[50]; 

    srand(time(0));

    printf("IN GA: Reading var.in file \n");
    read_var_input(var_input, lower_upper);
    printf("IN GA: Reading var.in file finished\n");

    
    /*printf("IN GA: Creating hyperplane file \n");
    if(access(hyper, F_OK) == -1) {
        create_hyperplane();
    }
    printf("IN GA: Creating hyperplane finished \n");
    */
    FILE *hyper_input, *create;
    char point[9];
    hyper_input = fopen("hyperplane.in","r");
    fgets(point,9,hyper_input);
    sscanf(point,"%d",&size_hyperplane);
    double hyperplane[size_hyperplane][nobj];
    fgets(point,9,hyper_input);
    for(i=0;i<size_hyperplane;i++) {
        if(i>0) fgets(point,9,hyper_input);
        for(j=0;j<nobj;j++) {
            fgets(point,9,hyper_input);
            sscanf(point,"%lf",&hyperplane[i][j]);
        }
    }    


    double var_old[nvar][N], var_new[nvar][N], var_pop[nvar][2*N];
    double obj_old[nobj][N], obj_new[nobj][N], obj_pop[nobj][2*N];
    double norm_obj_pop[nobj][2*N];
    int rank2[2*N], rank[N],dom_by_combined[2*N][2*N],dom_by_single[N][N];
    int front[2*N];
    int hyper_pointer[2*N];
    double hyper_distance[2*N];
    int point_parent[N];
    double var_parent[nvar][N],var_child[nvar][N];
    double obj_parent[nobj][N];
    double empty_obj[nobj][N];
    int current_pop_size;
   
  
    for(i=0;i<N;i++) {
        for(j=0;j<nobj;j++) {
            empty_obj[j][i] = 0.0;
        }
    }

    if(access(old_input, F_OK) != -1) { // returns TRUE if file exists
  	printf("Reading old population file \n");
        read_pop_input(old_input,var_old,obj_old);
	printf("Reading new population file \n");
        read_pop_input(new_input,var_new,obj_new); 

        /* combines old and new into single '2N' population */
        for(i=0;i<N;i++) {
            for(j=0;j<nvar;j++) var_pop[j][i] = var_old[j][i]; 
            for(j=0;j<nobj;j++) obj_pop[j][i] = obj_old[j][i];
        }
        for(i=0;i<N;i++) {
            for(j=0;j<nvar;j++) var_pop[j][N+i] = var_new[j][i];
            for(j=0;j<nobj;j++) obj_pop[j][N+i] = obj_new[j][i];
        }
        if(suppress_output==0) {
            printf("\nVariables and Objectives: \n");
            for(i=0;i<2*N;i++) {
                for(j=0;j<nvar;j++) {
                    printf("%lf ",var_pop[j][i]);
                }
                printf(" | ");
                for(j=0;j<nobj;j++) {
                    printf("%lf ",obj_pop[j][i]);
                }
                printf("\n");
            }
            printf("\n");
        }

        printf("Non-dominated sorting of population \n");

        non_dom_sort(2*N,obj_pop,rank2,dom_by_combined); // sorts population into ranks
	//FILE *ZerothFront = fopen("ZerothFront.txt", "a");
        //fprintf(ZerothFront, "%d\n", iterationNum);
	FILE *ZerothFront = fopen("ZerothFront.txt", "a");   
 
	printf("\nRank:\n");
        for(i=0;i<2*N;i++)
        {
                printf("%d ",rank[i]);
        }
        printf("\n\nDominated By:\n");
        int count = 0;
        for(i=0;i<2*N;i++)
        {
		count = 0;
                int domination = 0;
                for(j=0;j < 2*N;j++)
                {
                    printf("%d ",dom_by_single[i][j]);
                    if (dom_by_single[j][i] == 1) domination = 1;
                }
                if (domination == 0 && i >= N)
                {
                        count++;
                        fprintf(ZerothFront, "%d \t ", i+1 - N);
                        int iter = 0;
                        for (iter = 0; iter < nvar; iter++)
                                fprintf(ZerothFront, "%12.6f \t ", var_pop[iter][i]);
                        for (iter = 0; iter < nobj; iter++)
                                fprintf(ZerothFront, "%12.6f \t" ,obj_pop[iter][i]);
                        fprintf(ZerothFront, "\n");
                }
                printf("\n");
        }
        //fseek(ZerothFront, 0, SEEK_SET);
        fprintf(ZerothFront, "\n"); 
	fclose(ZerothFront);

        printf("Select population size\n");        
        /* selects population according to NSGA-III */
        current_pop_size = select_pop(rank2,point_parent,dom_by_combined,front); // populations with ranks until > N
        if(suppress_output==0) {
            printf("\nSelected Population:\n");
            printf("Current Population Size: %d\n",current_pop_size);
            printf("\nNew Population Pointers:\n");
            for(i=0;i<N;i++) {
                printf("%d ",point_parent[i]);
            }
            printf("\n\nFront:\n");
            for(i=0;i<2*N;i++) {
                printf("%d ",front[i]);
            }
            printf("\n");
        }
        printf("Normalize objectives \n");
        if (current_pop_size != N) {
            normalize(obj_pop,norm_obj_pop);
            if(suppress_output==0) {
                printf("\nNormalized Objectives\n");
                for(i=0;i<nobj;i++) {
                    for(j=0;j<2*N;j++) {
                        printf("%lf ",norm_obj_pop[i][j]);
                    }
                    printf("\n");
                }
                printf("\n");
                printf("\nNormalized Objectives\n");
                for(i=0;i<nobj;i++) {
                    for(j=0;j<2*N;j++) {
                        printf("%lf ",norm_obj_pop[i][j]);
                    }
                    printf("\n");
                }
                printf("\n");
            }
	    printf("Performing association operation\n");
            associate(norm_obj_pop,hyperplane,hyper_pointer,hyper_distance);
            if(suppress_output==0) {
                printf("ID#, Front, Nearest Pointers, and Distance:\n");
                for(i=0;i<2*N;i++) {
                    printf("%d %d %d %lf\n",i,front[i],hyper_pointer[i],hyper_distance[i]);
                }
                printf("\n");
            }
            printf("Performing niching operation\n");
            niching(current_pop_size,front,hyper_pointer,point_parent,hyper_distance);
            if(suppress_output==0) {
                printf("Filled Population Pointers:\n");
                for(i=0;i<N;i++) {
                    printf("%d ",point_parent[i]);
                }
                printf("\n");
            }     
        }

        /* translates pointer to selected population into arrays for variable and rank */
        if(suppress_output==0) printf("\nParent Variables and Objectives:\n");
        for(i=0;i<N;i++) {
            for(j=0;j<nvar;j++) {
                var_parent[j][i] = var_pop[j][point_parent[i]];
                if(suppress_output==0) printf("%lf ",var_parent[j][i]);
            }
            if(suppress_output==0) printf(" | ");
            for(j=0;j<nobj;j++) {
                obj_parent[j][i] = obj_pop[j][point_parent[i]];
                if(suppress_output==0) printf("%lf ",obj_parent[j][i]);
            }
            rank[i] = rank2[point_parent[i]];
            if(suppress_output==0) printf("\n");
        }
        if(suppress_output==0) printf("\n");

    } 
    else 
    {
        read_pop_input(new_input,var_new,obj_new); 
        if(suppress_output==0) {
            printf("\nVariables and Objectives: \n");
            for(i=0;i<N;i++) {
                for(j=0;j<nvar;j++) {
                    printf("%lf ",var_new[j][i]);
                }
                printf(" | ");
                for(j=0;j<nobj;j++) {
                    printf("%lf ",obj_new[j][i]);
                }
                printf("\n");
            }
            printf("\n");
    	}
       
	FILE *ZerothFront = fopen("ZerothFront.txt", "a");
	fprintf(ZerothFront, "%d\n", iterationNum); 
        non_dom_sort(N, obj_new,rank,dom_by_single); // sorts population into ranks
        
        printf("\nRank:\n");
        for(i=0;i<N;i++) 
	{
        	printf("%d ",rank[i]);
        }
        printf("\n\nDominated By:\n");
	int count = 0;
        for(i=0;i<N;i++) 
	{
		int domination = 0;
                for(j=0;j<N;j++) 
		{
                    printf("%d ",dom_by_single[i][j]);
		    if (dom_by_single[j][i] == 1) domination = 1;
                }
		if (domination == 0) 
		{
			count++;
			fprintf(ZerothFront, "%d \t ", i+1);
			int iter = 0;
			for (iter = 0; iter < nvar; iter++)
				fprintf(ZerothFront, "%12.6f \t ", var_new[iter][i]);
			for (iter = 0; iter < nobj; iter++) 
				fprintf(ZerothFront, "%12.6f \t" ,obj_new[iter][i]);
			fprintf(ZerothFront, "\n");
		}
                printf("\n");
        }
	//rewind(ZerothFront);
	fprintf(ZerothFront, "\n");
	fclose(ZerothFront);	
   
        //fprintf(ZerothFront, "\n\n");
        if(suppress_output==0) printf("\nParent Variables and Objectives:\n");
        for(i=0;i<N;i++) {
            for(j=0;j<nvar;j++) {
                var_parent[j][i] = var_new[j][i]; //copies var_new to var_parent for consistency
                if(suppress_output==0) printf("%lf ",var_parent[j][i]); 
            }
            if(suppress_output==0) printf(" | ");
            for(j=0;j<nobj;j++) {
                obj_parent[j][i] = obj_new[j][i];
                if(suppress_output==0) printf("%lf ",obj_parent[j][i]);
            }
            if(suppress_output==0) printf("\n");
        }
        if(suppress_output==0) printf("\n");
        
    }
    /*  SAVES PARENT POPULATION  */
    if(access(all_inputs, F_OK) == -1) {
        create = fopen(all_inputs,"w");
        fclose(create);
    }
    printf("Write ga_all.in file \n");
    write_pop_file(old_input,all_inputs,var_parent,obj_parent,1);
    /* generates child population via crossover and mutation */
    printf("Crossover operation \n");
    select_crossover(rank,var_parent,var_child, lower_upper);
    if(suppress_output==0) {
           printf("\nChild Post-Crossover:\n");
        for(i=0;i<N;i++) {
            for(j=0;j<nvar;j++) {
                printf("%lf ",var_child[j][i]);
            }
            printf("\n");
        }
    }

    printf("Mutation operation \n");
    select_mutants(var_parent,var_child,lower_upper);
    if(suppress_output==0) {

        printf("\nFinal Child Population:\n");
        for(i=0;i<N;i++) {
            for(j=0;j<nvar;j++) {
                printf("%lf ",var_child[j][i]);
            }
            printf("\n");
        } 
    }
    
    printf("Writes new population\n"); 
    /* Saves child population */
    write_pop_file(new_input,all_inputs,var_child,empty_obj,0);
    
    return(0);
}


int non_dom_sort(int pop,double y[][pop],int rank[pop],int dom_by[pop][pop]) { 
    int i,j,k;
    int obj_dom = 0; // number of objective over which one member dominates another
   
    //FILE *zeroth = fopen("ZerothFront.txt", "a"); 
    /* vector initialization */
    for(i=0;i<pop;i++) rank[i]=0;
    for(i=0;i<pop;i++) {
        for(j=0;j<pop;j++) dom_by[i][j]=0;
    }
    /* first pass pareto sort */
    for(i=0;i<pop;i++) { // loop through all members of the child + parent population 
	int domination = 0;
        for(j=i+1;j<pop;j++) { // loop through all other members for comparison purposes
            for(k=0;k<nobj;k++) { // loop through all objectives
                if(y[k][i]<y[k][j]) {
                    obj_dom += 1; // increments if objective is better for i than j
                } else if (y[k][i]>y[k][j]) {
                    obj_dom -= 1; // decrements if objective is worse for i than j
                }
            }
            if(obj_dom==nobj) { // i dominates j
                    dom_by[i][j] = 1; 
                    rank[j] += 1; 
            } else if (obj_dom==-nobj) { // j dominates i
                    dom_by[j][i] = 1;
                    rank[i] += 1;
            }
            obj_dom = 0; // reset counter
        }
    }
    return(0);
}

int select_pop(int rank[2*N],int new_pop[N],int dom_by[2*N][2*N],int front[2*N]) {
    int i,j;
    int current_rank_size = 0; // counter for current rank size
    int pop_size = 0; // counter for size of accepted population
    int rank_temp[2*N];
    for(i=0;i<2*N;i++) front[i]=0;
    for(i=0;i<N;i++) new_pop[i]=-1;
    for(i=0;i<2*N;i++) rank_temp[i]=rank[i]; // copies rank over to temp array


    //FILE *zeroth = fopen("ZerothFront.txt", "a");
    for(i=0;i<2*N;i++) {
        if(rank_temp[i] == 0) {
             front[i] = 1; // adds to first pareto front
             current_rank_size += 1;
	     //fprintf(zeroth, "%d \n", i);
        }
    }
    //rewind(zeroth);
    //fprintf(zeroth, "%d\n", current_rank_size);
    //fprintf(zeroth, "\n\n");

    while ((pop_size+current_rank_size) < N) {
        for(i=0;i<2*N;i++) {
            if(front[i]==1) { // if member i is part of the front
                front[i] = 0; // removes from current front
                rank_temp[i] = -1; // marks as removed from population
                new_pop[pop_size] = i; // add "pointer" to member i to final population
                pop_size += 1; // increment size counter
                for(j=0;j<2*N;j++) {
                    if(dom_by[i][j]==1) {
                        rank_temp[j] -= 1; // reduces rank of member dominated by removed pareto front member
                    }
                }       
            }
        }
        current_rank_size = 0;
        for(i=0;i<2*N;i++) {
            if(rank_temp[i] == 0) {
                front[i] = 1; // adds to first pareto front
                current_rank_size += 1;
            }
        }
    }
    return(pop_size);
}

int normalize(double obj[nobj][2*N],double norm_obj[nobj][2*N]) {
    int i,j,k,index[nobj],repeat;
    double min[nobj], stemp=0.0, s[2*N], smin=0.0, max[nobj][nobj], obj_dummy[nobj][2*N];
    double max_inv[nobj][nobj], hyper_intercepts[nobj];
    double det;
    double epsilon = 0.00000000001;
    double epsilon_prime;
    for(i=0;i<nobj;i++) { // initialization
        for(j=0;j<2*N;j++) {
            obj_dummy[i][j] = obj[i][j]; // holds objs for manipulation
            max[i][j] = 0.0;
        }
        min[i] = 0.0;
        hyper_intercepts[i] = 0.0;
        index[i] = 0;
    }
    for(i=0;i<nobj;i++) {
        for(j=0;j<2*N;j++) {
            if(j==0) {
                min[i] = obj_dummy[i][j]; // sets baseline for minimum
            } else if (obj[i][j]<min[i]) min[i] = obj_dummy[i][j]; // finds minimum value across combined population
        }
        for(j=0;j<2*N;j++) {
            obj_dummy[i][j] -= min[i]; // normalize objective values
        }
    }
    /* computes extreme point of each objective axis */

    for(i=0;i<nobj;i++) 
    { // loop through objective axes
        for(j=0;j<2*N;j++) 
        { // loop through population
            for(k=0;k<nobj;k++) 
            { // loop through objectives

                if ( i == k ) epsilon_prime = 1;
                else epsilon_prime = epsilon;
                
                if((obj_dummy[k][j]/epsilon_prime)>stemp) stemp = obj_dummy[k][j]/epsilon_prime; // stores maximum value among objectives

            }

            s[j] = stemp; // one maximum non-axis objective per member of population
            stemp = 0.0;
        }

        for(j=0;j<2*N;j++) 
        { // finds index of the minimum of the set of maximums
            if(j==0) 
            {
                for(k=0;k<i;k++) 
                    if(j==index[k]) repeat = 1;

                if(repeat==1) smin = (1/epsilon)*(1/epsilon);
                else 
                {
                    smin = s[j];
                    index[i] = j;
                }

                repeat = 0;
            }
            if(s[j]<smin) 
            {
                for(k=0;k<i;k++) 
                    if(j==index[k]) repeat = 1;
         
                if(repeat!=1) 
                {
                    smin = s[j];
                    index[i] = j;
                }
 
                repeat = 0;
            }
        }
        for(k=0;k<nobj;k++) max[k][i] = obj_dummy[k][index[i]];
    }
    
    det = determinant(max,nobj);

    if (det==0) {
           printf("\nInverse of maximum matrix is not possible\n");
    }

    cofactor(max,nobj,max_inv);

    for(i=0;i<nobj;i++) {
        for(j=0;j<nobj;j++) {
            hyper_intercepts[i] += max_inv[j][i];
        }
        hyper_intercepts[i] = 1/hyper_intercepts[i];
        for(j=0;j<2*N;j++) {
            norm_obj[i][j] = obj_dummy[i][j]/hyper_intercepts[i];
        }
    }
    /* INSERT ASSOCIATION TO USER PROVIDED HYPERPLANE HERE (?) */

    return(0);
}

int associate(double norm_obj[nobj][2*N],double hyperplane[nobj][size_hyperplane],int nearest_pointer[2*N],double nearest_distance[2*N]) {
    // input: normalized objective points (norm_obj), hyperplane
    int i,j,k;
    double dist[2*N][size_hyperplane]; // perpendicular distance to reference lines
    double dist_hold = 0.0; // perpendicular distance placeholder
    double vector_mult_acc = 0.0; // accumulates for vector multiplication
    double norm[nobj]; // 2-norm normalized hyperplane vector
    double two_norm = 0.0; // 2-norm of current vector (distance to vector from origin)
    double sing_norm_obj[nobj]; // single population member's normalized objectives

    FILE *fp = fopen("ClosestPoint.txt", "w");
    int nearestReference[size_hyperplane];
    for (i=0;i<size_hyperplane;i++) nearestReference[i] = 0;

    /* Initialization */
    for(i=0;i<size_hyperplane;i++) {
        for(j=0;j<2*N;j++) dist[j][i] = 0.0; 
    }
    for(i=0;i<nobj;i++) {
        norm[i] = 0.0;
        sing_norm_obj[i] = 0.0;
    }

    /* Association */
    for(i=0;i<2*N;i++) { // loop over population
        if(suppress_output==0) printf("In loop for population index %d\n", i);
        for(j=0;j<size_hyperplane;j++) { // loop over hyperplane 
            for(k=0;k<nobj;k++) { // loop over objectives
                two_norm += hyperplane[k][j]*hyperplane[k][j]; // accumulate for calculation of two_norm
            }  
            two_norm = sqrt(two_norm); // calculate two_norm
            for(k=0;k<nobj;k++) { // loop over objectives
                norm[k] = hyperplane[k][j]/two_norm; // populate 2-norm'd vector
                sing_norm_obj[k] = norm_obj[k][i]; // extract current normalized objectives for ease of reading
                vector_mult_acc += norm[k]*sing_norm_obj[k]; // accumulate for vector muliplication of above
            }
            for(k=0;k<nobj;k++) {
                dist_hold += (sing_norm_obj[k]-vector_mult_acc*norm[k])*(sing_norm_obj[k]-vector_mult_acc*norm[k]);
            }
            dist[i][j] = sqrt(dist_hold); // distance between hyperplane point and population member 
            if(j==0) {
                nearest_distance[i] = dist[i][j]; // initializes pointer to nearest hyperplane point and the distance
                nearest_pointer[i] = j;
            } else {
                if(dist[i][j]<nearest_distance[i]) {
                    nearest_distance[i] = dist[i][j]; // sets lowest distance and associated pointer
                    nearest_pointer[i] = j;
                }
            }
            two_norm = 0.0; // reset placeholder variables
            dist_hold = 0.0;
            vector_mult_acc = 0.0;   
        }
    }

    double minDist = 999999999.0;
    int minCounter=-1;
    for (i=0; i < size_hyperplane; i++)
    {
	minDist = 999999999.9;
	minCounter = -1;
	for (j=N; j<2*N;j++)
	{
	    if ( minDist >= dist[i][j] ) 
	    {
		minDist = dist[i][j];
		minCounter=j;
	    }
	}
	fprintf(fp, "%d  %d   %10.6f\n", i, minCounter-N, minDist);
    }

    return(0);
}

int niching(int pop_filled,int front[2*N],int nearest_pointer[2*N], int new_pop_point[N], double nearest_distance[2*N]) {
    int i,j,count = 0;
    int rho[size_hyperplane]; // amount each niche is filled
    int niche_min, num_niche_min=1;
    double unit_rand;
    int nth_niche;
    int niche_point = 0;
    int initialized = 0;
    double min_dist = 0.0;
    int new_point = -1;
    int loop = 0;
    for(i=0;i<size_hyperplane;i++) {
        rho[i] = 0; // initialization
    }
    for(i=0;i<pop_filled;i++) {
        rho[nearest_pointer[new_pop_point[i]]] += 1; // fills niches according to current population
    }
    
    
    while(pop_filled<N) 
    { // loops until parent population is... populated
        niche_min = rho[0]; // initializes minimum niche 
        for(j=1;j<size_hyperplane;j++) 
        {
            if(rho[j]==niche_min) num_niche_min += 1; // increments for multiple niches of smallest size
            if(rho[j]<niche_min) 
            {
                niche_min = rho[j]; // resets if it encounters a smaller niche
                num_niche_min = 1;
            }
        }
        unit_rand = (double)rand()/(double)(1.0 + RAND_MAX); // random number [0,1)
        nth_niche = (int)(unit_rand*(double)num_niche_min); // randomly selects which niche to fill 
     
        for(j=0;j<size_hyperplane;j++) 
	{ 
            if(rho[j]==niche_min)
             { // if niche is a minimum niche
                if(nth_niche==count) 
                {
                    niche_point = j; // if niche is selected niche, set pointer and break
                    break;
                }
                if(nth_niche>count) count += 1; // else, increment niche counter
            }
        }

        if(niche_min == 0) 
        { // if minimum niche is empty
            for(j=0;j<2*N;j++) 
            { // loop over population
                if(front[j]==1)
                { // if point is member of the front
                    /* if the niche of the member of the front matches the selected niche */
                    if(nearest_pointer[j]==niche_point) 
                    { 
                        if(initialized == 0) 
                        { // if this is first member of niche
                            min_dist = nearest_distance[front[j]]; // sets first distance to minimum
                            new_point = j; // sets first pointer to shortest distance member of niche
                            initialized = 1; // sets boolean to true
                        } 
                        else if (nearest_distance[j]<min_dist) 
                        {
                            min_dist = nearest_distance[j]; // sets new distance to minimum
                            new_point = j; // sets new pointer to shortest distance member of niche 
                        }
                    }
                }
            }
            if(new_point<0) 
            { // if no point in niche exists within front
                rho[niche_point] = 2*N+1; // niche  removed from consideration
            } 
            else 
            {
                new_pop_point[pop_filled] = new_point; // point added to parent population
                pop_filled += 1; // parent population size incremented
                rho[niche_point] += 1; // population of niche incremented
                front[new_point] = 0; // point removed from front
            }
        } 
        else 
        { // if miniumum niche is not empty
            count = 0; // reset counter
            for(j=0;j<2*N;j++) 
            {
                if(front[j]==1) 
                { // if point is member of front
                    /* if the niche of the member of the front matches the selected niche */
                    if(nearest_pointer[j]==niche_point) count += 1; // increment counter of number of members in niche
                }
            }
            if(count==0) 
            { // no point in niche exists within front
                rho[niche_point] = 2*N+1; // removes niche from consideration
            } 
            else 
            {
                unit_rand = (double)rand()/(double)(1.0 + RAND_MAX); // random number [0,1)
                new_point = (int)(unit_rand*(double)count); // randomly picks member to add to population
                count = 0; // reset counter
                for(j=0;j<2*N;j++) 
                {
                    if(front[j]==1) 
                    {   // if point is member of front
                        /* if the niche of the member of the front matches the selected niche */
                        if(nearest_pointer[j]==niche_point) 
                        {
                            if(count==new_point) 
                            {
                                new_pop_point[pop_filled] = j; // adds to population
                                pop_filled += 1;
                                rho[niche_point] += 1;
                                front[j] = 0;
                                break;
                            }
                            count += 1; // increment counter of number of members in niche
                        }
                    }
                }   
            }
        }
        num_niche_min = 1;
        count = 0;
        initialized = 0;
        min_dist = 0.0;
        new_point = -1;
    }
    return(0);
}

int select_crossover(int r[],double var_parent[nvar][N],double var_child[nvar][N], struct bound lower_upper[]) {

	int  i,j,n;
	double alpha = 0.5, pi;
	int crossover_n = crossover_P*N;
    	int mating_pair[2];
	int xlow,xhigh;
	double rsum=0,rmax=0,rand_num,rprob[N],raprob[N];

	for (i=0;i<N;i++) 
	{
		if(rmax<(r[i]+1)) rmax = r[i]+1; // finds latest rank included in population
	}
	for (i=0;i<N;i++) 
	{
        	rsum+=(rmax-r[i]-1); // sums up ranks (+1 accounts for first rank = 0)
        	for(j=0;j<nvar;j++) var_child[j][i] = 0.0; // initializes var_child
        }
    
	for(i=0;i<N;i++)  rprob[i]=(rmax-r[i]-1)/rsum; // lower ranks have higher probability of acceptance
	raprob[0]=rprob[0];
	for(i=1;i<N;i++) raprob[i]=raprob[i-1]+rprob[i];

	for (i=0;i< crossover_n;i++) 
	{
	
        		/* selects two parents for mating via roulette-wheel selection */
        		j = 0;
		
			do 
			{ 
            			rand_num=(double)rand()/(double)RAND_MAX; // random number [0,1]
            			n=(int)(rand_num*(N+1));
				if (n >= N) n = N-1;
				if (j == 0) { mating_pair[j]=n; j++;}
				else if (j == 1)
				{
					if (n != mating_pair[j]) { mating_pair[j+1] = n; j++; }
				}
        		} while (j < 2);
  			for(j=0;j<nvar;j++) 
			{
				do
				{
          				//rand_num=(double)rand()/(double)RAND_MAX;
            				/* crossover between mating pair via weighted addition */
            				//var_child[j][i] = rand_num*var_parent[j][mating_pair[0]]+(1-rand_num)*var_parent[j][mating_pair[1]];
					double ui = (double)rand()/(double)RAND_MAX;
					double yi = (1+2*alpha)*ui - alpha;
					var_child[j][i] = yi*var_parent[j][mating_pair[0]]+(1-yi)*var_parent[j][mating_pair[1]];
       				} while (var_child[j][i] < lower_upper[j].x || var_child[j][i] > lower_upper[j].y);
			}
	}	
	return(0);
}

int select_mutants(double var_parent[nvar][N], double var_child[nvar][N],struct bound lower_upper[]) {
    int i,j;
    int n_mutants = (int)ceil(mutation_P*N); // calculates portion of population that are mutants
    int n_var_mutate = (int)ceil(mutation_rate*nvar); // number of variables to mutate per mutation
    double unit_rand1, unit_rand2; // holds generated value of rands scaled 0-1
    int parent_mutate, var_mutate; // holds scaled values of rand
    double BM_uniform_1, BM_uniform_2, BM_normal; // rands for box-muller method - normal distributed rand
    int n_start = N-n_mutants;
    double sigma;
    double mutant;

    printf("Population to mutate %d \n", n_mutants);

    for(i=0;i<n_mutants;i++) {
        unit_rand1 = (double)rand()/(double)(1.0 + RAND_MAX);
        parent_mutate = (int)(unit_rand1*(double)(N)); // picks parent to mutate
        printf("Population to mutate %d \n", parent_mutate);
        for(j=0;j<nvar;j++) var_child[j][n_start+i] = var_parent[j][parent_mutate]; // copies variables of parent
        printf("Total variables in a population to be mutated %d \n", n_var_mutate);

        for(j=0;j<nvar;j++) {
            
            var_mutate = j;//(int)(unit_rand2*(double)(nvar)); // picks variable to be mutated
            printf("Selected variable is %d \n", var_mutate);
            printf("Value of selected populations %d variable %d is %lf\n", n_start+i, var_mutate, var_child[var_mutate][n_start+i]);
            printf("Lower and upper bound of variable %d is %lf and %lf\n", var_mutate, lower_upper[var_mutate].x, lower_upper[var_mutate].y);
            do {
                sigma = 0.5*(lower_upper[var_mutate].y-lower_upper[var_mutate].x);
                BM_uniform_1 = (double)rand()/(double)RAND_MAX;
                BM_uniform_2 = (double)rand()/(double)RAND_MAX;
                BM_normal = sqrt(-2*log(BM_uniform_1))*cos(2*M_PI*BM_uniform_2);
                mutant = var_child[var_mutate][n_start+i] + sigma*BM_normal;
            } while (mutant < lower_upper[var_mutate].x || mutant > lower_upper[var_mutate].y);
            printf("Lower and Upper bound of selected variable %d is %lf and %lf \n", var_mutate, lower_upper[var_mutate].x, lower_upper[var_mutate].y);
            var_child[var_mutate][n_start+i] += sigma*BM_normal;
            printf("Value of the selected population %d and corresponding variable %d after mutation is %lf\n", n_start+i, var_mutate, var_child[var_mutate][n_start+i]);
        }
    }
    return(0);
}

int read_var_input(char filename[],struct bound bounds[]) {
    FILE *input;
    int i;
    input = fopen(filename,"r");

    if (input == NULL) {
        printf("%s file does not exist\n",filename); // null pointer handler
        exit(0);
    }    

    fscanf(input,"%d",&N);
    fscanf(input,"%d",&nvar);
    fscanf(input,"%d",&nobj);
    fscanf(input,"%lf",&crossover_P);
    fscanf(input,"%lf",&mutation_P);
    fscanf(input,"%lf",&mutation_rate);
    fscanf(input,"%d",&nDivision);
    //fscanf(input,"%d",&size_hyperplane);
    for(i=0;i<nvar;i++) {
        fscanf(input,"%lf",&bounds[i].x);
    }
    for(i=0;i<nvar;i++) {
        fscanf(input,"%lf",&bounds[i].y);
    }
    fscanf(input,"%d",&suppress_output);

    fclose(input);
    return(0);
}


int read_pop_input(char fileName[],double x[][N],double y[][N]){

	FILE *fin; 
	int i, j;
    int count=0;
    char point[10];
	fin = fopen(fileName,"r"); // open input file stream
    if (fin == NULL) {
        printf("No population input file exists\n"); // null pointer handler
        exit(0);
    }    
	fscanf(fin, "%d",&it_num);
	for(i=0;i<N;i++){ // loop over population
		for (j = 0; j < nvar; j++) { // loop over line in population
            fscanf(fin,"%lf",&x[j][i]);
            //printf("%0.4lf ",x[j][i]);
        }
        //printf("\n");
        for(j=0;j<nobj;j++) {
            fscanf(fin,"%lf",&y[j][i]);
            //printf("%0.4lf ",y[j][i]);
        }
        //printf("\n");
        //fgets(point,10,fin);
  	}
    
	fclose(fin);

	return(0);
}   

int write_pop_file(char fileName[],char allName[],double x[][N], double y[][N], int type) {
    FILE *output,*all_out; 
    int i, j;
    output = fopen(fileName,"w"); // creates output file, overriding previous file if applicable
    all_out = fopen(allName,"a");

    fprintf(output,"%d\n",N);
    if(type==1) fprintf(all_out,"%d\n",N);
    for(i=0;i<N;i++) {
        for(j=0;j<nvar;j++) {
            fprintf(output,"%f ",x[j][i]); // prints variables in a line
            if(type==1) fprintf(all_out,"%f ",x[j][i]);
        }
        if(type==1) {
            for(j=0;j<nobj;j++) {
                fprintf(output,"%f ",y[j][i]); // prints objectives in a line for parent only
                fprintf(all_out,"%f ",y[j][i]);
            }
        }
        fprintf(output,"\n"); // line break
        if(type==1) fprintf(all_out,"\n");
    }
    fclose(output);
    fclose(all_out);
    return(0);
}


void create_hyperplane() {
    int i,j;
    int n_part = 0; // number of partitions (minus one [loop])
    int test;
    double sing_part[nobj]; // temp for single partition
    /* Variables to hold all partitions of size no more than nobj */
    int part_divis[nDivision][nobj]; // if num of divisions is larger than num of objectives
    int part_obj[nobj][nobj]; // if num of objectivies is larger than num of divisions
    FILE *output,*create; 
    create = fopen("hyper_wd.in","w"); // creates hyperplane w/ duplicates file
    fclose(create);
    output = fopen("hyper_wd.in","r+"); // read/write to hyperplane file
    int size_hyperplane; // size of hyperplane
    int obj_factorial = 1; 

    for(i=1;i<nobj+1;i++) {
        obj_factorial *= i; // calculates number of permutations per partition
    }

    if(nDivision>nobj) { // if the number of divisions is larger than the number of objectives
        for(i=0;i<nobj;i++){
            for(j=0;j<nDivision;j++) part_divis[j][i] = 0; // initialization
        }
        n_part = partition(part_divis,n_part);
    } else {
        for(i=0;i<nobj;i++){
            for(j=0;j<nobj;j++) part_obj[j][i] = 0; // initialization
        }
        n_part = partition(part_obj,n_part);
    }

    size_hyperplane = obj_factorial*(n_part+1); // size of hyperplane is number of unique partitions * permutaions per partition
    for(i=0;i<n_part;i++) {
        if(nDivision>nobj) {
            for(j=0;j<nobj;j++) sing_part[j] = part_divis[i][j]; //nDivision; // extract single partition
        } else {
            for(j=0;j<nobj;j++) sing_part[j] = part_obj[i][j]; //nDivision; // extract single partition
        }
        permute(output,sing_part,0);  
    }
    double hyperplane[size_hyperplane][nobj]; // hold variable for hyperplane
    for(i=0;i<size_hyperplane;i++) {
        for(j=0;j<nobj;j++) hyperplane[i][j] = 0.0; // initialization
    }
    remove_duplicates(output,hyperplane,size_hyperplane); // removes duplicates from hyperplane document
    fclose(output);

}


int partition(int f_p[][nobj],int index) {
    int i,k = 0;  // Index of last element in a partition
    int p[nDivision]; // An array to store a partition
    for(i=0;i<nDivision;i++) p[i] = 0;
    
    p[k] = nDivision;  // Initialize first partition as nDivisionber itself
    int rem_val;
    // This loop first prints current partition, then generates next
    // partition. The loop stops when the current partition has all 1s
    while (1) {
        if(nDivision>nobj) {  // save current partition
            if(p[nobj]!=0) break; // partition has size > nobj
        }
        for(i=0;i<nDivision;i++) f_p[index][i] = p[i]; // copy single partition to set
        index++; // increment set counter
        rem_val = 0;
        while (k >= 0 && p[k] == 1) {
            rem_val += p[k];
            k--;
        }
 
        if (k < 0)  return(index); // if k < 0, all the values are 1 so there are no more partitions
        
        p[k]--; // Decrease the p[k] found above and adjust the rem_val
        rem_val++;
 
        /*  If rem_val is more, then the sorted order is violated.  Divide
            rem_val in different values of size p[k] and copy these values at
            different positions after p[k] */
        while (rem_val > p[k]) {
            p[k+1] = p[k];
            rem_val = rem_val - p[k];
            k++;
        }
        p[k+1] = rem_val; // Copy rem_val to next position and increment position
        k++;
    }
    return(index); // return number of partitions
}

void swap(double arr[],int a,int b) { //function to swap the variables
    int temp;
    temp = arr[a];
    arr[a] = arr[b];
    arr[b] = temp;
    return;
}

void permute(FILE *output,double sing_part[],int start) {
    int i,j;
    double read_array[nobj],unique_check;
    

    if(start==nobj) { // if a full permutation is filled      
        for(i=0;i<nobj;i++) {
            fprintf(output,"%.5lf ",sing_part[i]/nDivision); // print permutation to hyperplane doc
        }
        fprintf(output,"\n");
        return;
    }
    for(i=start;i<nobj;i++) {
        swap(sing_part,i,start); //swaps numbers
        permute(output,sing_part,start+1); /* fixes one first digit and calls permutation on the rest */
        swap(sing_part,i,start); // swaps back
    }

}

void remove_duplicates(FILE *output, double hyperplane[][nobj],int size_hyperplane) {
    int i,j,k=0,current_line=0,current_hyper=0,match=0;
    int duplicates[size_hyperplane+1];
    for(i=0;i<size_hyperplane+1;i++) duplicates[i] = -1;
    char line[200],point[9];
    double temp[nobj];
    int file_length;
    FILE *true_hyper;

    rewind(output); // rewinds document to beginning for reading
    while(1) {
        SKIP: if(feof(output)) break; // checks for EoF, breaks if found. goto position
        if(current_line>0) { // for all but first loop
            fgets(point,9,output); // throws out "\n" between partitions
        }
        for(i=0;i<nobj;i++) { // extracts a partition and saves to temporary storage
            fgets(point,9,output);
            //printf("%s\n",point);
            sscanf(point,"%lf",&temp[i]);
        }

        if(current_hyper>0) { // all but first loop
            for(i=0;i<current_hyper;i++) { // loop through current extracted partitions
                for(j=0;j<nobj;j++) {
                    if(temp[j]==hyperplane[i][j]) match++; // check is partition is duplicate
                }
                if(match==nobj) { // if partition is a duplicate
                    match = 0; // reset match counter
                    duplicates[k] = current_line; // flag current position in doc
                    current_line++;
                    k++; // increment duplicate counter
                    goto SKIP; // use goto to skip step in which partition is added to hyperplane 
                }
                match = 0; // reset match counter
            }
        }
        current_line++; // increment counter for current line in document
        for(i=0;i<nobj;i++) {
            hyperplane[current_hyper][i] = temp[i]; // update hyperplane
        }
        current_hyper++;

    }
    true_hyper = fopen("hyperplane.in","w"); // open doc for hyperplane w/o duplicates 
    file_length = current_line - k; 
    fprintf(true_hyper,"%d\n",file_length); // print ultimate size of hyperplane for reference
    
    current_line = 0; // reset variables
    k = 0;
    rewind(output); // rewind document

    while(fgets(line,200,output) != NULL) { // while not EoF

        if(current_line==duplicates[k]) {
            k++; // move to next duplicate, do not print
        } else {
            fputs(line,true_hyper); // add to document
        }
        current_line++;

    }
    
    fclose(true_hyper);

}


/* CALCULATION OF MATRIX INVERSE */

/*For calculating Determinant of the Matrix */
double determinant(double a[nobj][nobj], double k) {
  double s = 1, det = 0, b[nobj][nobj];
  int i, j, m, n, c;
  if (k == 1)
    {
     return (a[0][0]);
    }
  else
    {
     det = 0;
     for (c = 0; c < k; c++)
       {
        m = 0;
        n = 0;
        for (i = 0;i < k; i++)
          {
            for (j = 0 ;j < k; j++)
              {
                b[i][j] = 0;
                if (i != 0 && j != c)
                 {
                   b[m][n] = a[i][j];
                   if (n < (k - 2))
                    n++;
                   else
                    {
                     n = 0;
                     m++;
                     }
                   }
               }
             }
          det = det + s * (a[0][c] * determinant(b, k - 1));
          s = -1 * s;
          }
    }
 
    return (det);
}
 
void cofactor(double num[nobj][nobj], double f, double inv[nobj][nobj]) {
 double b[nobj][nobj], fac[nobj][nobj];
 int p, q, m, n, i, j;
 for (q = 0;q < f; q++)
 {
   for (p = 0;p < f; p++)
    {
     m = 0;
     n = 0;
     for (i = 0;i < f; i++)
     {
       for (j = 0;j < f; j++)
        {
          if (i != q && j != p)
          {
            b[m][n] = num[i][j];
            if (n < (f - 2))
             n++;
            else
             {
               n = 0;
               m++;
               }
            }
        }
      }
      fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
    }
  }
  transpose(num, fac, f, inv);
}
/*Finding transpose of matrix*/ 
void transpose(double num[nobj][nobj], double fac[nobj][nobj], double r, double inverse[nobj][nobj]) {
  int i, j;
  double b[nobj][nobj], d;
 
  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
         b[i][j] = fac[j][i];
        }
    }
  d = determinant(num, r);
  for (i = 0;i < r; i++)
    {
     for (j = 0;j < r; j++)
       {
        inverse[i][j] = b[i][j] / d;
        }
    }

}


