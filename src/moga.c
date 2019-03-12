#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include "mpi.h"

struct Range
{
        double min_x, min_f, max_x, max_f;

};

// function in which the majority of running a single thread is completed
void processor_run(int id, int iter_num, char folder[], char inputs[]);

// initializes population based on hardcoded (in function) initial conditions
void initialize_population(int population_num);
// modifies gulp forcefield input file
void modify_input(double A[3], double rho[3], double B[3], double lambda[3], double gamma_3b[2], double theta);
// read output from gulp 
void read_output(double* err_a,double* elast);
// calculate error in phonon dispersion curve
double ErrorPhononDispersion(int state);
// spile interpolation for use in dispersion curve error
double CubicSplineInterp(int N, double *x, double *f, struct Range *limit, int ref_N, double *ref_x, double *ref_f, int id);
// read in band data for dispersion curve, contains hardcoded values for band data structures
void ReadData(int numBands, double x_vasp_data[][20000], double f_vasp_data[][20000], double x_gulp_data[][20000], double f_gulp_data[][20000], int state);
// creates a copy of a file, and populates it accordingly in case of gulp inputs
void file_copy(char file[], char copyTo[], int bands, int state);

int main(int argc, char *argv[]) {
    
    
    int myid; /* My rank */
    int nprocs; /* Number of processors */
    int iteration_num = 500;   // number of training iterations
    int population_num = 300; // population size for ga.in training
    int restart = 0; // set to 1 if restarting an optimization
    int initialized; // tracks MPI instantiation 

    if(restart==0) {
        initialize_population(population_num); // initializes population
    } 

    
    int i,j,line; // iterators
    char c,variables[500],folder[200];
    char test[200];
    
    for(j=0;j<iteration_num;j++) { // recursively optimizes for number of iterations specified


        MPI_Initialized(&initialized);
        if (!initialized) {
            MPI_Init(&argc, &argv);
        }
        MPI_Comm_rank(MPI_COMM_WORLD, &myid); // starts MPI
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

      
        FILE *input;
        input = fopen("ga.in","r"); // opens ga.in file for reading
        
        if (input == NULL)
        {
            printf("Cannot open file \n"); // null handler
            exit(0);
        }

	//printf("Collecting ga_lines \n");
        for (i=myid;i<population_num;i+=nprocs) { // loops through population list, assigning members to processors

            sprintf(folder,"%d",i); // creates a folder to hold all analysis on current line

            // gets line specified by loop and reads it into memory
            
            line = 0;
            rewind(input);

            while (line!=(i+1)) {
                if((c = fgetc(input)) == '\n') line++;
            }
            fgets(variables,500,input);
            
            // end
            printf("%s\n", variables);
            processor_run(myid,j, folder,variables); // runs analysis on line input

        }
        fclose(input);

	    MPI_Barrier(MPI_COMM_WORLD);

        if (myid == 0) {
            FILE *ga_input,*ga_line;
            //sprintf(test,"ga_%d.in", j);
            ga_input = fopen("ga.in","w");

            if (ga_input == NULL)
            {
                printf("Cannot open file \n");
                exit(0);
            }

            fprintf(ga_input,"%d\n",population_num);
	    printf("Concatenating ga input files \n");
            for (i=0;i<population_num;i++) {
                sprintf(folder,"%d",i);
                chdir(folder);

                ga_line = fopen("ga_line","r");
                fgets(variables,500, ga_line);
                fputs(variables, ga_input);

                fclose(ga_line);

                chdir("..");
                // code that removes directories files generated above
                //printf("Deleting Directories\n");
                char rm_dir[200];
                if(j!=iteration_num-1) {
                    sprintf(rm_dir,"rm -rf %s", folder);
                    system(rm_dir);
                }
            }
            fclose(ga_input);
        
	    printf("%d Running ga executable \n", j+1);
            char cmdexec1[200];
	        sprintf(cmdexec1,"./ga ga.in %d", j+1);
	        system(cmdexec1);
            // char cmdexec[200];
	        // sprintf(cmdexec, "cp value.d ga.in");
	        // system(cmdexec);
            printf("Iteration = %d done \n\n", j);

        }
    MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Finalize();

    return 0;
}


void file_copy(char file[], char copyTo[], int bands, int state) 
{ // creates a copy of a file

    FILE *in_stream, *out_stream;
    char c;
    
    int regions[3][3];
    regions[0][0] = 92;		regions[0][1] = 53;		regions[0][2] = 106;
    regions[1][0] = 87;		regions[1][1] = 51;		regions[1][2] = 100;
    regions[2][0] = 83;		regions[2][1] = 48;		regions[2][2] = 96;    

    in_stream = fopen(file,"r");
    if (in_stream == NULL)
    {
        printf("Cannot open file %s \n", file); // null pointer handler
        exit(0);
    }
    
    out_stream = fopen(copyTo, "a");
    if (out_stream == NULL)
    {
        printf("Cannot open file %s \n", copyTo); // null pointer handler
        exit(0);
    }
    
    c = fgetc(in_stream);
    while (c != EOF)
    {
        fputc(c, out_stream);
        c = fgetc(in_stream);
    }
    
    switch (bands)
    {
	case 1:
		// G_M region
		fprintf(out_stream,"\ndispersion 1 %d \n", regions[state][0]);
		fprintf(out_stream,"0.0 0.0 0.0 to 0.5 0.0 0.0\n");
		break;
	case 2:
		// M_K region
		fprintf(out_stream,"\ndispersion 1 %d \n", regions[state][1]);
		fprintf(out_stream,"0.5 0.0 0.0 to 0.333333 0.333333 0.0\n");
		break;
	case 3:
		// K_G region
		fprintf(out_stream,"\ndispersion 1 %d \n", regions[state][2]);
		fprintf(out_stream,"0.333333 0.333333 0.0 to 0.0 0.0 0.0\n");
		break;
	//default:
		//printf("Copied %s to %s \n", file, copyTo);

    }

    if (bands > 0) 
    {
        fprintf(out_stream,"output phon phonon\n");
        fprintf(out_stream,"shrink 5 5 5\n\n");
    }
    
    fclose(in_stream);
    fclose(out_stream);
    
}

void modify_input(double A[3], double rho[3], double B[3], double lambda[2], double gamma_3b[2], double theta) {
    
    double theta0 = 80.8150; // geometry-dependent constant
    
    FILE *ff;
    char filename[] = "forcefield"; // Input file name
    ff = fopen(filename, "r+");
    
    char c, str[200]; //current character, string to hold number
    int line=0; // current line, current position in file
    
    // moves to 5th line for first modification

    while (line !=4) 
    {
        if((c = fgetc(ff)) == '\n') line++;
    }

    int i, target;
    
    for (i=0;i<3;i++) { //modifies 2-body terms 3-pairs in total Se-Se, Mo-Se, Mo-Mo
        
        fseek(ff,11,SEEK_CUR); //goes to first modified value
        sprintf(str,"%6.3lf   %6.3lf   %7.4lf",A[i],rho[i],B[i]); // Modifies A, rho and B to every pair corresponding to the tuple in ga.in file

        fputs(str,ff); // updates first line w/ correct values, current cursor position at end of B
    
        target = line+1;

        while (line != target) 
        {
        	if((c = fgetc(ff)) == '\n') 
			line++; // Loops moves the counter to next line
        } 
    }
    
    while (line !=11) 
    {
        if((c = fgetc(ff)) == '\n') line++; // moves to 11th line for 3-body terms
    }
  
    for (i=0;i<2;i++) 
    { 
        fseek(ff,10,SEEK_CUR); //goes to first modified value

	sprintf(str,"%7.4lf  %7.4lf %6.3lf %6.3lf", lambda[i], theta, gamma_3b[i], gamma_3b[i]);
	
        fputs(str,ff); // updates first line w/ correct values, current cursor position at end of rmax
        
        target = line + 1;
        while (line != target) {
            if((c = fgetc(ff)) == '\n') line++; // moves to next line
        }
    }

    fclose(ff); //close input file
}

void read_output(double* err_a,double* elast) {

    FILE *af;
   
    char output[] = "afterfit.out";
    af = fopen(output,"r"); // read output file
    
    if (af == NULL)
    {
        printf("Cannot open file %s \n", output); // null pointer handler
        exit(0);
    }
    
    
    double c_11, c, el_c;
    double lat_const = 13.97, elast_real = 128.7; // elastic constant in GPa
    double elast_calc;

    char tempStr1[100], tempStr2[100];
    double temp1, temp2, temp3;

    char lines[200], search_string[] = "  Comparison of initial and final structures : ", ch;
    while (fgets(lines,200,af) != NULL) { // search output file for ^ test string
        int i;
        if (strstr(lines,search_string)) {
	    //printf(lines);
            break;
        }
    }
    
    int line_down=0;
    while (line_down !=5) {
        if((ch = fgetc(af)) == '\n') line_down++; // moves to constant "a"
    }
    fgets(lines,200,af); // reads line in
    int len = strlen(lines);
    double err_dum, err_dum1; // dummy error
    //printf("Extracting err_a\n");
    sscanf(lines, "%s %lf %lf %lf %s %lf", tempStr1, &temp1, &temp2, &temp3, tempStr2, &err_dum); // assigns value of percent error in a
    fgets(lines,200,af); // reads line in
    sscanf(lines, "%s %lf %lf %lf %s %lf", tempStr1, &temp1, &temp2, &temp3, tempStr2, &err_dum1);
    *err_a = sqrt(pow(err_dum,2) + pow(err_dum1,2));

    //printf("Extracting err_a done\n");    

    /*line_down = 0;
    while (line_down !=1) { // moves down to elastic constant "c" (for elastic fit)
        if((ch = fgetc(af)) == '\n') line_down++;
    }*/
    
    fgets(lines,200,af);

    //printf("Extracting C\n");
    sscanf(lines,"%s %lf %lf %lf %s %lf", tempStr1, &temp1, &c, &temp2, tempStr2, &temp3); // assigns value of elastic constant "c" (c_final)
    //printf("Extracting C done\n");

    line_down = 0;
    while (line_down !=48) { // moves down to elastic constant C11
        if((ch = fgetc(af)) == '\n') line_down++;
    }
    
    fgets(lines,200,af);

    //printf("Extracting C11\n");
    int temp;
    //printf("%s\n", lines);
    sscanf(lines, "%*d %lf", &c_11); // assigns value of elastic constant C11
    //printf("Extracting C11 done\n");    

    //printf("c = %f, c_11 = %f\n", c, c_11);

    elast_calc = c_11*c*(2/lat_const); // elastic constant for comparison
    *elast = fabs(elast_calc-elast_real)/elast_real*100; // percent error in elastic constant
    
    fclose(af);
    
    return;
}

void ReadData(int numBands, double x_vasp_data[][20000], double f_vasp_data[][20000], double x_gulp_data[][20000], double f_gulp_data[][20000], int state)
{
        int regionsPerBand = 3, pointsPerRegion_vasp = 101;
	int regions[3][3];
	regions[0][0] = 92;         regions[0][1] = 53;             regions[0][2] = 106;
        regions[1][0] = 87;         regions[1][1] = 51;             regions[1][2] = 100;
        regions[2][0] = 83;         regions[2][1] = 48;             regions[2][2] = 96;

        int points_GM_gulp = regions[state][0], points_MK_gulp = regions[state][1], points_KG_gulp = regions[state][2]; 

	double temp_x, temp_y;

        double phonon_shift[4] = {0, regions[state][0], regions[state][0] + regions[state][1], regions[state][0] + regions[state][1] + regions[state][2]};
        double conv_x = 1000, conv_y = 33.35641;
        int band_offset = 3, band_segm_length = 21, gulp_offset = 3;
	char *band_file;

        int i = 0, j = 0, k = 0, lineNum_vasp = 0, lineNum_gulp = 0, bands_vasp = 0, bands_gulp = 0;

        const char *phonon_file[3] = {"G_M/phonon.disp", "M_K/phonon.disp", "K_G/phonon.disp"};
         
        char lines_gulp[200], lines_band[200], lines_band_ph[200];

        FILE *vasp, *gulp[4];

        if (state == 0) vasp = fopen("../../UTIL/band0.dat", "r");
        else if (state == 1) vasp = fopen("../../UTIL/band1.dat", "r");
	else if (state == 2) vasp = fopen("../../UTIL/band2.dat", "r");

        /* Reading Vasp File band.dat */
	//vasp = fopen(band_file, "r");
        if (vasp == NULL) { printf("Unable to open file for reading dispersion data .... exiting \n"); exit(1); }

        while (lineNum_vasp < band_offset) { fgets(lines_band, 200, vasp); lineNum_vasp++; }

        for (i = 0; i < numBands; i++)
        {
                for (j = 0; j < regionsPerBand; j++)
                {
                        fgets(lines_band, 200, vasp);
                        lineNum_vasp++;
                        while(strcmp(lines_band,"\n") != 0)
                        {
                                sscanf(lines_band, "%lf %lf", &temp_x, &temp_y);
				x_vasp_data[i][k] = temp_x;
				f_vasp_data[i][k] = temp_y;
				//printf("i %d j %d = %10.6f %10.6f \n", i + 1, k + 1 , x_vasp_data[i][k], f_vasp_data[i][k]);
				k++;
                                fgets(lines_band, 200, vasp);
                                lineNum_vasp++;
                        }
			//printf("Line Number = %d \n", lineNum_vasp);
		}
		fgets(lines_band, 200, vasp);
                lineNum_vasp++;
                k=0;
        }

	int curPoints = 0, prevPoint = 0;
        for (i = 0; i < regionsPerBand; i++)
        {
                lineNum_gulp = 0;
                gulp[i] = fopen(phonon_file[i], "r");
                if (gulp[i] == NULL) { printf("Unable to open file %s for reading dispersion data .... exiting \n", phonon_file[i]); exit(1); }

                while (lineNum_gulp < gulp_offset) {fgets(lines_gulp, 200, gulp[i]); lineNum_gulp++; }

                if (i == 0) curPoints = points_GM_gulp;
                else if (i == 1) curPoints = points_MK_gulp;
                else if (i == 2) curPoints = points_KG_gulp;
                for (j = 0; j < curPoints; j++)
                {
                        for (k = 0; k < numBands; k++)
                        {
                                fgets(lines_gulp, 200, gulp[i]);
                                sscanf(lines_gulp, "%lf %lf", &temp_x, &temp_y);
                                x_gulp_data[k][prevPoint + j] = temp_x + prevPoint;
                                f_gulp_data[k][prevPoint + j] = temp_y;
			}
		}

		prevPoint += curPoints;
        }

        fclose(vasp); fclose(gulp[0]); fclose(gulp[1]); fclose(gulp[2]);

	return;
}

double CubicSplineInterp(int N, double *x, double *f, struct Range *limit, int ref_N, double *ref_x, double *ref_f, int id)
{

        int i = 0, j = 0, k = 0, rangeFound = 0;
        char fileName[100];
        sprintf(fileName, "interp-%d.txt", id);

        FILE *fw = fopen(fileName, "w");

        if (fw == NULL) { printf("Cannot open file to write.... exiting \n"); exit(1);}

        double e[N], h[N], d[N], r[N], z[N];
        double df, dx, xs, x1, y, sqErr = 0.0;
        char buf[100];

        N = N-1;
        for (i = 0; i <= N-1; i++)
        {
                if (i < N-1)
                {
                        e[i] = 2.0*(x[i+2] - x[i]);
                        h[i] = x[i+2] - x[i+1];
                        d[i] = x[i+1] - x[i];
                        r[i] = 6*(f[i+2] - f[i+1])/h[i] - 6*(f[i+1] - f[i])/d[i];
		}
                else if (i == N-1)
                {
                        d[i] = x[i+1] - x[i];
                }
        }

        for (i = 0; i <= N-3; i++)
        {
                df = d[i]/e[i];
                e[i+1] = e[i+1] - df*h[i];
                r[i+1] = r[i+1] - df*r[i];
        }

        for (i = N-3; i >= 0; i--)
        {
                df = h[i]/e[i+1];
                r[i+1] = r[i] - r[i+1]*df;
        }

        for (i = 0; i <= N-2 ; i++)
                z[i+1] = r[i]/e[i];
        z[0] = 0.0; z[N] = 0.0;

        for (j = 0; j < ref_N; j++)
        {
	
                x1 = ref_x[j];
                if (limit->min_x > x1 || limit->max_x < x1) continue;

		k = 0;
		rangeFound = 0;
		while (k < N)
                {
                        if (x[k] < x1 && x[k+1] > x1)
                        {
                                i = k;
                                rangeFound = 1;
                                break;
                        }

                        else if (x[k] == x1)
                        {
                                y = f[k]; break;
                        }

                        else if (x[k+1] == x1)
			{
                                y = f[k+1]; break;
                        }
			k++;
                }
                if (rangeFound == 0)
                {
                        if (x[k] == x1 || x[k+1] == x1)
                        {
                                sqErr += pow(y-ref_f[j], 2);
                                fprintf(fw, "%10.4f \t %10.4f \n", x1, y);
                        }
                        continue;
                }
                y = -z[i]*pow( (x1-x[i+1]),3 )/(6.0*d[i]) +
                     z[i+1]*pow( (x1-x[i]),3 )/(6.0*d[i]) +
                    (f[i+1]/d[i] - z[i+1]*d[i]/6.0)*(x1 - x[i]) +
                    (-f[i]/d[i] + z[i]*d[i]/6.0)*(x1 - x[i+1]);

		sqErr += pow(y-ref_f[j], 2);

                fprintf(fw, "%10.4f \t %10.4f \n", x1, y);
	}

	fclose(fw);
        return sqrt(sqErr);
}

double ErrorPhononDispersion(int state) {

        int BandNum = 0, numBands = 36, i = 0;
	double normalizingFactor = 0.0;
        int regionsPerBand = 3, pointsPerRegion_vasp = 101;

        double x_vasp_data[numBands][20000], f_vasp_data[numBands][20000];
        double x_gulp_data[numBands][20000], f_gulp_data[numBands][20000];

        double x[20000], f[20000];
        double ref_x[regionsPerBand*pointsPerRegion_vasp], ref_f[regionsPerBand*pointsPerRegion_vasp];
        double Error[numBands], w[numBands];

        for (i = 0; i < numBands; i++)
	{
		if (i < 3) { w[i] = 1.0; normalizingFactor += 1.0;}
		else if (i >= 3) { w[i] = 0.1; normalizingFactor += 0.1;}
	}

        double SquaredError = 0.0;

        double min_x = 10000000000.0, max_x = -10000000000.0, min_f = 10000000000.0, max_f = -10000000000.0;

        struct Range *lim = (struct Range*)malloc(sizeof(struct Range));

        ReadData(numBands, x_vasp_data, f_vasp_data, x_gulp_data, f_gulp_data, state);
  
        for (BandNum=0; BandNum < numBands; BandNum++)
        {
        	for(i = 0; i < 500; i++)
        	{
        		x[i] = x_gulp_data[BandNum][i];
        		f[i] = f_gulp_data[BandNum][i];
        		if (x[i] <= min_x) min_x = x[i];
                        if (x[i] >= max_x) max_x = x[i];
			if (f[i] <= min_f) min_f = f[i];
                        if (f[i] >= max_f) max_f = f[i];
                }

                for(i = 0; i < regionsPerBand*pointsPerRegion_vasp; i++)
                {
                        ref_x[i] = x_vasp_data[BandNum][i]*1000.0;
                        ref_f[i] = f_vasp_data[BandNum][i]*33.35641;
                }

                lim->min_x = min_x; lim->min_f = min_f; lim->max_x = max_x; lim->max_f = max_f;

                Error[BandNum] = CubicSplineInterp(530, x, f, lim, regionsPerBand*pointsPerRegion_vasp, ref_x, ref_f, BandNum);

                SquaredError += w[BandNum]*Error[BandNum];
	}
        //for (BandNum=0; BandNum < numBands; BandNum++) printf("%10.6f \t ", Error[BandNum];
	return SquaredError;
	//return (SquaredError/(normalizingFactor*regionsPerBand*pointsPerRegion_vasp));
}

void initialize_population(int population_num) {
    
    FILE *population;
    
    population = fopen("ga.in", "w");
    
    double frac_perturb = 0.10;
    double rand_frac;
    fprintf(population, "%d\n", population_num);
    
    //initial guesses for the 3 As, 3 rhos, 3 Bs, 2 lambdas, 2 gamma_3b
    //double variables[14] = {2.0, 6.0, 3.0, 1.6, 0.4, 0.4, 29.0, 12.0, 33.0, 66.0, 23.0, 1.70, 5.70, 80.0};
    //double variables[14] = {7.311202,3.182849,1.751833,0.263154,0.276447,34.792870,7.942051,27.888579,87.012161,17.268118,2.098628,5.993827,77.305379}
    double variables[14] = {2.3,5.710990,2.890970,1.6533,0.524863,0.308964,37.596895,13.340301,33.898837,63.438181,18.979162,2.007693,5.657942,78.174734};
    double rand_var[14]; // array for holding random perturbations of variables
    
    int i,j,k;
    
    srand(time(NULL));
    
    for (i=0;i<population_num;i++) {
        for (j=0;j<14;j++) {
	    if (j <= 2) rand_frac = rand() / (double)(RAND_MAX)*(2*0.03) - 0.03;
            else if (j > 2 && j <= 12 ) rand_frac = rand() / (double)(RAND_MAX)*(2*frac_perturb) - frac_perturb;
	    else if (j == 13) rand_frac = rand() / (double)(RAND_MAX)*(2*0.05) - 0.05;
            rand_var[j] = (1-rand_frac)*variables[j];
            fprintf(population,"%lf ",rand_var[j]);
        }
        fprintf(population,"\n");
    }
    
    fclose(population);
}

void processor_run(int id, int iter_num, char folder[], char inputs[]) 
{
    
    //initialize variables (make global only if using private copies in this subroutine (OpenMPI))
    double A[3], rho[3], B[3], lambda[2], gamma_3b[2], theta; //variables modified in the force field file
    double err_a, elast, chi_sq[3]; 
    char file[200], path[200], c,dc[20];
    int bands;

    double weight_a = 1.0, weight_c11 = 1.0;

    mkdir(folder,S_IRWXU);
    // copy cell, forcefield, etc. into folder
    
    const char *gulp_in[6] = {"UTIL/cell1", "UTIL/cell2", "UTIL/cell3"};
    const char *states[3] = {"Compressed", "ZeroStrain", "Expanded"};
    const char *segment[3] = {"G_M","M_K","K_G"};
    char StateDir[200], CopyTo[200];
    int i,j,k;    
   
    sscanf(inputs,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &A[0],&A[1],&A[2],&rho[0],&rho[1],&rho[2],&B[0],&B[1],&B[2],&lambda[0],&lambda[1],&gamma_3b[0],&gamma_3b[1], &theta);

    for (i = 0; i < 3; i++)
    {
	sprintf(path, "%s/%s", folder, states[i]);
        for (j = 0; j < 2; j++)
        {
		if (j == 0)
		{
			strcpy(file, gulp_in[i]);	// Copy gulp file to file string
			mkdir(path, S_IRWXU);
			sprintf(CopyTo, "%s/%s", path, "cell");
		}
		else if (j==1)
		{
			strcpy(file, "UTIL/forcefield");
			sprintf(CopyTo, "%s/%s", path, "forcefield");
		}
		file_copy(file, CopyTo, 0, i);	// 0 as last argument is just a place holder for creating a cell file at given path. It will overwritten for different regions
	}
    
    	chdir(path);   //enter folder

    	modify_input(A,rho,B,lambda,gamma_3b, theta); // modifies "forcefield" file

    	FILE *gulp_input;

    	gulp_input = fopen("in.gulp","w"); // create gulp input file
    
    	if (i == 0) fprintf(gulp_input, "phon nofreq\n");
	else if (i == 1) fprintf(gulp_input, "optim relax conp comp phon nofreq\n");
	else if (i == 2) fprintf(gulp_input, "phon nofreq\n");

    	fclose(gulp_input);
    
    	char file1[100], file2[100], path1[200], path2[200];
   
    	strcpy(file1,"cell");
    	strcpy(path1,"in.gulp");
    	bands = 0;
    	file_copy(file1,path1,bands,i); // appends "cell" to in.gulp
    
        // Modify this file 
        strcpy(file2,"forcefield");
        file_copy(file2,path1,bands,i); // appends "forcefield" to in.gulp

       char file0[3][100], path0[3][100];
 
       for(k=0;k<3;k++) 
       {
        
		strcpy(dc,segment[k]);
        	mkdir(dc,S_IRWXU); // creates directory for segment
        
        	strcpy(file0[k],"in.gulp");
        	strcpy(path0[k],dc);
        	strcat(path0[k],"/");
        	strcat(path0[k],file0[k]);  // copies in.gulp into new directory
        	bands = k+1;
        
        	file_copy(file0[k],path0[k],bands,i);
        
        	chdir(dc); // enter segment directory
        
		system("/home/pv-02/hpc-23/kris658/SOFTWARE/GULP/Src/Linux/gulp < in.gulp > afterfit.out"); // runs gulp

        	chdir(".."); // exits new directory
       }

       chdir(segment[0]); // arbitrary choice of segment for objectives 1 & 2 (same for all)
       //printf("%d Extracting elastic constant \n", id);
       if (i == 1) 
       {
		read_output(&err_a,&elast); // extracts objective 1 & 2: error in lattice constant a and error in elastic constant
		err_a *= weight_a; elast *= weight_c11;
       }

       //printf("Extracting elastic constant --done\n");    

       chdir(".."); // moves back to 'folder' directory to read all segments for phonon dispersion chi squared calculation

       //chi_sq = phonon_disp(); 
       //printf("Computing error in dispersion \n");
       chi_sq[i] = ErrorPhononDispersion(i);
       chdir("../..");
    }
    chdir(folder);
    FILE *output;

    output = fopen("ga_line","w"); 
    fprintf(output,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",A[0],A[1],A[2],rho[0],rho[1],rho[2],B[0],B[1],B[2],lambda[0],lambda[1],gamma_3b[0],gamma_3b[1],theta,err_a,elast,chi_sq[0],chi_sq[1],chi_sq[2]);
    
     chdir("..");
     fclose(output);
}


