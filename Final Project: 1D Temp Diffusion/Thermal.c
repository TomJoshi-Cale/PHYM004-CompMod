/*
 Title:		PHYM004 Final Project: Project F - Thermal Diffusion in 1D
 Author:	Tom Joshi-Cale tj294 670000715
 Style:		Linux Kernel
 		https://www.kernel.org/doc/html/v4.10/process/coding-style.html
 License:	Public Domain
*/
static const char *VERSION = "2.0.0";
static const char *REV_DATE = "25-Mar-2021";
/*
 * This code, Thermal.c, is designed to model the Temperature of a 1-D system as
 * a function of time - and produce a graph of Temperature vs Radius for chosen
 * snapshots of simulation time.
 * The code runs in one of two modes, Spherical mode (chosen by running with
 * argument -e or --egg) or cylindrical mode (argument -r or --rod). 
 * The code is compiled as 
 *	$ gcc Thermal.c -lm -lgsl -lgslcblas -o Thermal
 * and is run using
 *	$ ./Thermal -r
 * or
 *	$ ./Thermal -e
 *
 * The gnuplot script provided with submission will only work for v5.0 and later
 *
 * ACKNOWLEDGEMENTS
 * The cylindrical case was written with the assistance of 'Numerical Solutions 
 * of Partial Differential Equations' by Louise Olsen-Kettle, which is available
 * at https://espace.library.uq.edu.au/eserv/UQ:239427/Lectures_Book.pdf
 *
 * The spherical case was written while consulting 'The Science of Boiling an 
 * Egg' by Charles D. H. Williams, available from 
 * https://newton.ex.ac.uk/teaching/CDHW/egg/
 *
 * The Gnu Scientific Library documentation was used when implementing matrices 
 * and vectors, and for the matrix inversion using linear algebra. The relevant
 * documentation can be found at 
 * https://www.gnu.org/software/gsl/doc/html/vectors.html and
 * https://www.gnu.org/software/gsl/doc/html/linalg.html while the matrix-vector
 * multiplication function used on line 176 is described at
 * https://www.gnu.org/software/gsl/doc/html/blas.html#level-2
 *
 * Eggs harmed in the making of this code: 1 (cracked while measuring radius)
 */

/*
 Date 		Version 	Comments
 ----		-------		--------
 18-Mar-2021	1.0.0		Merged Rod.c and Egg.c codes into 
 				this final mixed code

 20-Mar-2021	1.1.0		Customised Rod case to take variable time inputs
 				based on snapshot_percent constant, and 
 				generalised create_a() function and 
 				run_simulation() function to hold for both egg 
 				and rod

 21-Mar-2021	1.2.0		Code now runs for both Rod and Egg, and allows 
 				user to set what % of egg is yolk. Gain 
 				parameter is now treated as a vector, as are 
 				kappa, heat cap and rho.

 21-Mar-2021	1.2.1		Changed format of code so 'r'od and 'e'gg cases 
 				set up initial conditions, but run_simulation() 
 				function is called outside of cases.

 23-Mar-2021	1.2.2		Took setting initial conditions out of main.
 
 25-Mar-2021	2.0.0		Layout tweaks, and ability to use 'dumb' gnuplot
 				terminal type added to code.
*/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <errno.h>
#include <gsl/gsl_linalg.h>

/*
=============================
======== DEFINITIONS ========
=============================
*/
#define MAX_COMMAND_LENGTH 255
// [s]	seconds in a minute
#define MINS_TO_SECS 	60.0
// [s]	seconds in a (Julian) year
#define YEARS_TO_SECS 	31557600.0
#define FILE_NAME 	"temp.dat"
#define SCRIPT 		"Temp.script"
//The Gnuplot terminal type. Try "wxt" if you're having problems, or try 'dumb'! 
#define GNU_TERM_TYPE	"qt"	//Default: qt

#define LENGTH(array) (int) (sizeof(array) / sizeof(array[0]))

typedef enum object { EGG, ROD } Object; //Used as a flag throughout code

enum ERROR_CODES {
	SUCCESS = 0,
	BAD_ARGUMENTS = -1,
	BAD_RUNTIME = -2,
	BAD_SNAPSHOT = -3,
	BAD_PROPERTIES = -4,
	OPEN_ERROR = -5,
	BAD_TYPE = -6
};
/*
=============================
========= CONSTANTS =========
=============================
*/
// USER-DEFINED INFORMATION - change to change behaviour of code
// Which percentages of simulation time to plot (max 100)
static const double snapshot_percent[] = {1, 10.0, 50.0, 100.0};
//Number of spacial grid points (default 99)
static const double SPACIAL_RESOLUTION = 99;
//Number of time steps (default 1000)
static const double TIME_RESOLUTION = 1000;

// Rod Information
// [years] - Simulation time for the rod.
static const double rod_run_time = 100.0;
// [%] - radius of rod as a % of simulation length
static const double rod_r_percent = 25.0;
// [%] - boundaries of material properties as a % of the rod radius
static const double rod_boundaries[] = {100.0};
/* [cm^2 /s] - Thermal Diffusivity of Rod - number of elements must match those
 * in rod_boundaries
 */
static const double kappa_rod[] = {0.634};

//Egg Information
// [minutes] - Simulation time for the egg.
static const double egg_run_time = 15.0;
// [%] - boundaries of material properties of the egg as a % of egg radius
static const double egg_boundaries[] = {33.3, 100.0};
//[cm^2 /s] - Thermal Diffusivity of egg yolk, then egg white
static const double kappa_egg[] = {1.22e-3, 1.40e-3};

//ROD CONSTANTS
//[K] - Environmental Temp for Rod
static const double environment_temp = 300;
//[] - Number of source conditions required i.e. T_rod, rod_radius and tau_0
static const int n_source_conditions = 3;

//EGG CONSTANTS
//[K] - 100C (water boiling point)
static const double water_temp = 100.0;
//[K] - 4 C (Egg fresh from fridge)
static const double egg_temp = 4.0;


/*
=============================
========= FUNCTIONS =========
=============================
*/

// Prints the version and revision date to stdout at runtime.
static void info()
{
	printf("#VERSION: %s\n#REVISION DATE: %s\n", VERSION, REV_DATE);
}

// Check that the user-defined information is valid
static void validate_runtime(double run_time, int n_kappa , int n_boundaries)
{
	if (run_time <= 0) {
		fprintf(stderr, 
			"Run-time must be positive.\nInputted runtime: %lg\n",
			run_time);
		exit(BAD_RUNTIME);
	}
	if (snapshot_percent[0] > 100 || snapshot_percent[0] < 0 ) {
		fprintf(stderr,
			"All snapshot percents must be less than 100%% and greater than 0%%\n");
		exit(BAD_RUNTIME);
	}
	for (int i = 1; i < LENGTH(snapshot_percent); i++) {
		if (snapshot_percent[i] > 100) {
			fprintf(stderr,
				"All snapshot percents must be less than 100%%\n");
			exit(BAD_RUNTIME);
		}
		if (snapshot_percent[i] <= snapshot_percent[i-1]) {
			fprintf(stderr,
				"Ensure snapshot_percent array is in ascending order, and has no duplicate elements\n");
			exit(BAD_RUNTIME);
		}
	}

	if (n_kappa != n_boundaries) {
		fprintf(stderr, 
			"Diffusivity does not have the same number of elements at boundaries.\nPlease ensure the properties match\n");
		exit(BAD_PROPERTIES);
	}
}

// Opens a file, and checks it opened successfully
static FILE *open_and_check(char *filename, char *opentype) {
	FILE *f = fopen(filename, opentype);
	if (!f) {
		fprintf(stderr, "Error: Could not open file '%s'/\n", filename);
		exit(OPEN_ERROR);
	}
	return f;
}

// Calculates the gain parameter (denoted s in Olsen-Kettle) and returns it
static gsl_vector *gain_parameter(double n, double rc, double dt, gsl_vector *r,
				   gsl_vector *kappa, const double *type_kappa)
{
	gsl_vector *s = gsl_vector_calloc(n);
	int i_region = 0;
	double dr = rc / (n + 1);

	for (int i = 0; i < n; i++) {
		gsl_vector_set(r, i, dr + (i * (((rc - (2 * dr)) / (n - 1)))));
		if (gsl_vector_get(r, i) < ((rod_boundaries[i_region] / 100) * rc)) {
			gsl_vector_set(kappa, i, type_kappa[i_region]);
		} else {
			i_region++;
			gsl_vector_set(kappa, i, type_kappa[i_region]);
		}
		gsl_vector_set(s, i, 
			(gsl_vector_get(kappa, i) * dt) / (dr * dr));
	}
	return s;
}

/*
 * allocates memory for and populates the matrix which is denoted A in 
 * Olsen-Kettle.
 */
static gsl_matrix *create_a(int n, gsl_vector *s, int factor)
{
	gsl_matrix *a = gsl_matrix_calloc(n, n);

	for (int i=1; i <= n; i++) {
		for (int j=1; j <= n; j++) {
			if (i==j) {
				gsl_matrix_set(a, i-1, j-1,
					1 + 2 * gsl_vector_get(s, i-1));
			} else if (i == (j - 1)) {
				gsl_matrix_set(a, i-1, j-1,
					 -gsl_vector_get(s, i-1) - 
				 	gsl_vector_get(s, i-1) / (factor * i));
			} else if (i == (j + 1)) {
				gsl_matrix_set(a, i-1, j-1,
					-gsl_vector_get(s, i-1) + 
					gsl_vector_get(s, i-1) / (factor * i));
			} else {
				gsl_matrix_set(a, i-1, j-1, 0);
			}
		}
	}

	return a;
}

// Sets the boundary conditions for a new system, returns them as gsl_vector b
static gsl_vector *set_boundary_conditions(gsl_matrix *a, gsl_vector *s, 
						double n, int type)
{
	gsl_vector *b = gsl_vector_calloc(n);

	switch (type) {
	case ROD:
		/* Implement Neumann r=0 boundary condition:
		 * dT/dr = 0 at r=0, so T(r=0, t) = T(r=1, t)
		 * Which changes the 0,0th element of matrix A
		 */ 
		gsl_matrix_set(a, 0, 0, 1 + gsl_vector_get(s, 0) + 
			(gsl_vector_get(s, 0) / 2));

		/* Implement Dirichlet r=rc boundary condition:
		 * T(r=rc, t) = T_environment
		 */
		gsl_vector_set(b, n-1, -(-gsl_vector_get(s, n-1) - 
			gsl_vector_get(s, n-1) / (2 * n)) * environment_temp);
		break;

	case EGG:
		/* Implement Dirichlet boundary condition
		 * T(r = rc) = water_temp at all times
		 */
		gsl_vector_set(b, n - 1, -((-gsl_vector_get(s, n-1) - 
			(gsl_vector_get(s, n-1) / (n))) * water_temp));

		break;

	default:	
		fprintf(stderr,
			"In function 'set_boundary_conditions', type not recognised\n");
		exit(BAD_TYPE);
	}

	return b;
}

// Sets the initial temperature distrubution for each system
static gsl_vector *set_initial_temperature(double n, gsl_vector *r, double rc, 
						int type)
{
	gsl_vector *T_tk = gsl_vector_calloc(n);
	switch (type) {
	case ROD:
		gsl_vector_set_all(T_tk, environment_temp);
		break;

	case EGG:
		for (int i = 0; i < n; i++) {
			if (gsl_vector_get(r, i) < rc) {
				gsl_vector_set(T_tk, i, egg_temp);
			} else {
				gsl_vector_set(T_tk, i, water_temp);
			}
		}
		break;

	default:
		fprintf(stderr, 
			"In function 'set_intial_temperature', type not recognised\n");
		exit(BAD_TYPE);
	}

	return T_tk;
}

/*
 * Allocates memory for and calculates the inverse of a matrix using gsl
 * implementation of LU Decomposition, then frees the memory the original matrix
 * is occupying.
 */
static gsl_matrix *invert(gsl_matrix *input)
{
	int signum;

	gsl_matrix *inverse = gsl_matrix_calloc(input->size1, input->size2);
	gsl_permutation *perm = gsl_permutation_alloc(input->size1);
	gsl_linalg_LU_decomp(input, perm, &signum);
	gsl_linalg_LU_invert(input, perm, inverse);
	gsl_permutation_free(perm);
	gsl_matrix_free(input);

	return inverse;
}
/*
 * Takes in the initial conditions for the problem, and iterates the solver over
 * time. Uses the Backwards Euler method as outlined in Olsen-Kettle to solve 
 * the problem.
 */
static void run_simulation(double t, double m, double dt, gsl_matrix *inv_a, 
				gsl_vector *source_vector, gsl_vector *T_tk, 
				gsl_vector *b, gsl_vector *kappa, 
				double *source_con, gsl_matrix *solutions) 
{
	double n = inv_a->size1;
	int isolutions = 0; //current solutions

	gsl_vector *source = gsl_vector_calloc(n);
	gsl_vector *v = gsl_vector_calloc(n);
	gsl_vector *T_tk_1 = gsl_vector_calloc(n);

	for (int k = 0; k <= m; k++) {
		t = t + dt;

		for (int i = 0; i < n; i++) {
			gsl_vector_set(source, i, 
				(source_con[1] * expf( -t / source_con[2]) / 
				(source_con[0] * source_con[0])) * 
				gsl_vector_get(source_vector, i));
			gsl_vector_set(v, i, 
				gsl_vector_get(T_tk, i) + gsl_vector_get(b, i) +
				(gsl_vector_get(kappa, i) * dt *
				gsl_vector_get(source, i)));
		}

	//Multiplies matrix inv_a with vector v, and stores result in T_tk_1
		gsl_blas_dgemv(CblasNoTrans, 1, inv_a, v, 0, T_tk_1); 

		for (int i = 0; i < n; i++) {
			gsl_vector_set(T_tk, i, gsl_vector_get(T_tk_1, i));
		}

		double check = floor(m / (100 / snapshot_percent[isolutions]));
		if (k == check) {
			for (int i = 0; i < n; i++) {
				gsl_matrix_set(solutions, i, isolutions,
					gsl_vector_get(T_tk_1, i));
			}
			isolutions++;
		}
	}
	gsl_vector_free(source);
	gsl_vector_free(v);
	gsl_vector_free(T_tk_1);
}

static void write_file(int n, double tf, int n_snaps, gsl_vector *r, 
			gsl_matrix *solutions, int type)
{
	FILE *f = open_and_check(FILE_NAME, "w");
	
	//Write the header line
	fprintf(f, "r");
	for (int i = 0; i < n_snaps; i++) {
		if (type == ROD) {
			fprintf(f, "\t \"%lg years\"",
				(snapshot_percent[i] / 100) * tf);
		}
		else {
			fprintf(f, "\t \"%lg mins\"",
				(snapshot_percent[i] / 100) * tf);
		}
	}
	fprintf(f, "\n");

	//Write the data to the file
	for (int i=0; i < n; i++) {
		fprintf(f, "%lg", gsl_vector_get(r, i));
		for (int j = 0; j < n_snaps; j++) {
			fprintf(f, "\t %.7lg", gsl_matrix_get(solutions, i, j));
		}
		fprintf(f, "\n");
	}
	fclose(f);
}
/*
========================
========= MAIN =========
========================
*/

int main(int argc, char **argv)
{
	info();
	char command[MAX_COMMAND_LENGTH];
	int c;
	if (LENGTH(snapshot_percent) == 0) {
		fprintf(stderr, 
			"No snapshot times detected.\nSnapshots can be defined at the top of the code\n");
		exit(BAD_SNAPSHOT);
	}
	int type = -1; //gets set to 0 for Egg or 1 for Rod.
	double source_conditions[n_source_conditions];

	double rc, n; //spacial variables
	double t, tf, m, dt; //temporal variables
	gsl_matrix *a;
	gsl_vector *kappa, *s, *b, *r, *source_vector, *T_tk;

	while (1)
	{
		static struct option long_options[] = 
		{
			{"egg",		no_argument, 0, 'e'},
			{"rod", 	no_argument, 0, 'r'},
			{0, 0, 0, 0}
		};

		int option_index = 0;

		c = getopt_long (argc, argv, "er",
				 long_options, &option_index);

		if (c == -1)
			break;

		switch (c) {
		case 0: 
			if (long_options[option_index].flag != 0)
				break;
			printf("option %s", long_options[option_index].name);
			if (optarg)
				printf(" with arg %s", optarg);
			printf("\n");
			break;

		case 'r':
			type = ROD;
			validate_runtime(rod_run_time, LENGTH(kappa_rod),
					 LENGTH(rod_boundaries));

			/* ==========================
			 * === INITIAL CONDITIONS ===
			 * ==========================
			 * Comment format //[Units] - Physical Meaning - Olsen-Kettle
			 */
			rc = 100; //[cm] - length simulated - r_c
			n = SPACIAL_RESOLUTION; //[] - number of grid points - n

			t = 0; //[s] - Starting time
			tf = rod_run_time; //[years] - Simulation time - T_f
			m = TIME_RESOLUTION; //[] - Time resolution - m
			dt = (tf * YEARS_TO_SECS) / m; //[secs] - time spacing - \Delta t

			r = gsl_vector_calloc(n);
			kappa = gsl_vector_calloc(n);

			s = gain_parameter(n, rc, dt, r, kappa, kappa_rod); //

			a = create_a(n, s, 2);

			b = set_boundary_conditions(a, s, n, type);

			gsl_vector_free(s);

			
			/* ==========================
			 * ===== SOURCE HEATING =====
			 * ==========================
			 */			

			// Define Source Heating term
			double rod_radius = (rod_r_percent / 100) * rc; //[cm] - Radius of Nuclear rod - rc
			double T_rod = 1.0; //[K] - Initial Temp of Rod above the environment - T_{rod}
			double tau_0 = 100.0 * YEARS_TO_SECS; //[secs]

			source_conditions[0] = rod_radius; 	//[cm]
			source_conditions[1] = T_rod;		//[K]
			source_conditions[2] = tau_0;		//[secs]

			//Source_vector is 1 inside rod and 0 elsewhere, as source heating only occurs inside rod
			source_vector = gsl_vector_calloc(n); 

			for (int i = 0; i < n; i++) {
				if (gsl_vector_get(r, i) < rod_radius)
					gsl_vector_set(source_vector, i, 1);
			}

			/* ===========================
			 * === INITIAL TEMPERATURE ===
			 * ===========================
			 */			

			/* Implement Initial Temperature - at t=0 all points at T_env:
			 * T(r, 0) = T_env
			 */
			T_tk = set_initial_temperature(n, r, rc, type);

			break;

		case 'e':
			type = EGG;
			validate_runtime(egg_run_time, LENGTH(kappa_egg),
					 LENGTH(egg_boundaries));

			/* ==========================
			 * === INITIAL CONDITIONS ===
			 * ==========================
			 * Comment format //[Units] - Physical Meaning
			 */

			//spacial
			rc = 2.5; //[cm] - radius of Spherical Egg
			n = SPACIAL_RESOLUTION; //[] - Number of grid points


			//temporal
			t = 0; //[s] - Initial time
			tf = egg_run_time; //[minutes] - Simulation time
			m = TIME_RESOLUTION; //[] - time resolution
			dt = (tf * MINS_TO_SECS) / m; //[s] - time separation

			r = gsl_vector_calloc(n);
			kappa = gsl_vector_calloc(n);

			
			s = gain_parameter(n, rc, dt, r, kappa, kappa_egg);
			
			a = create_a(n, s, 1);

			b = set_boundary_conditions(a, s, n, type);

			gsl_vector_free(s);

			/* ==========================
			 * ===== SOURCE HEATING =====
			 * ==========================
			 */

			//Egg has no internal source heating, so Source Vector stays at 0
			source_vector = gsl_vector_calloc(n);
			for (int i = 0; i < n_source_conditions; i++) {
				source_conditions[i] = 1; //These get multiplied by zero, so have no physical meaning
			}

			/* ===========================
			 * === INITIAL TEMPERATURE ===
			 * ===========================
			 */

			 /* Within the egg, initial temp is egg_temp
			 * and outside the egg, initial temp is water_temp
			 */
			T_tk = set_initial_temperature(n, r, rc, type);

			break;
			
		case '?':
			fprintf(stderr,
				"Allowed arguments are '-e' for Egg and '-r' for Rod\n");
			exit(BAD_ARGUMENTS);
			break;

		default: 
			abort();
		}
	}
	if (type == -1) {
		errno = -1;
		fprintf(stderr, "Input '-e' for Egg and '-r' for Rod\n");
		exit(BAD_ARGUMENTS);

	}
	
	gsl_matrix *inv_a = invert(a);

	// Initialise a matrix to store the required time snapshots in
	gsl_matrix *solutions = gsl_matrix_calloc(n, LENGTH(snapshot_percent));

	run_simulation(t, m, dt, inv_a, source_vector, T_tk, b, kappa, 
			source_conditions, solutions);

	gsl_vector_free(source_vector);
	gsl_vector_free(b);
	gsl_vector_free(T_tk);
	gsl_matrix_free(inv_a);
	gsl_vector_free(kappa);

	//write selected time snapshots into a text file
	write_file(n, tf, LENGTH(snapshot_percent), r, solutions, type);
	
	gsl_matrix_free(solutions);
	gsl_vector_free(r);

	snprintf(command, sizeof(command),
		 "gnuplot -e \"type=%i; filename='%s'; termtype='%s'\" %s",
		 type, FILE_NAME, GNU_TERM_TYPE, SCRIPT);
	system(command);

	if (optind < argc) {
		printf("non-option ARGV-elements: ");
		while (optind < argc)
			printf("%s ", argv[optind++]);
		putchar('\n');
	}

	exit(SUCCESS);
}
