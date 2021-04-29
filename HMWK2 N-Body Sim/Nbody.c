/*
 Title:		PHYM004 Homework 2
 Author:	Tom Joshi-Cale tj294 670000715
 Style:		Linux Kernel - https://www.kernel.org/doc/html/v4.10/process/coding-style.html
*/
static const char *VERSION = "1.1.0";
static const char *REV_DATE = "13-Dec-2020";
/*
 Date 		Version 	Comments
 ----		-------		--------
 15-Nov 	0.0.1		Began code by combining getopt_long example from
 				https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html
 				and the ReadOrbits.c code by C. D. H. Williams

 30-Nov		0.2.0		Added the ability to read in the timestep as a double after 
 				the data file.

 03-Dec		0.3.0		Added an implementation of the velocity-verlet function, and the calc_acceleration
 				function. Shooting bodies into space at the moment, I believe that the
 				acceleration is being calculated incorrectly.

 08-Dec		0.4.0		Added code to overwrite previous contents of output files so I don't
 				need to re-write them each time I run the code. I am planning to eventually
 				plot all files in one text file, making the create_output_files() function 
 				obsolete. Also, re-wrote calc_accleration to correctly calculate acceleration
 				Finally, added code to execute a gnuplot script and automatically plot
 				the calculated orbits.

 09-Dec		0.5.0		Added the energy_calc() function, which calculates the potential, kinetic and total
 				energy of each body, and outputs it into a 'bodyname_energy.dat' file, where it can
 				be plotted using gnuplot and my 'energy.script' script.
 
 11-Dec		0.6.0		Code now completes all the required tasks for Velocity-Verlet, 
 				(reading in from text file, being able to choose between VV and RK4, 
 				Calculating total energy, plotting orbits via Gnuplot), and I have used
 				it to model a figure 8 stable orbit solution for the three body problem
 				(starting conditions found in fig8.txt, and output in fig8.gif), and it
 				has been used to simulate Sun-Earth-Moon system (starting conditions:
 				3Body.txt, output: 3Body.gif), and the full solar system (starting 
 				conditions: Horizons.txt, output: Solar.gif).

 12-Dec		1.0.0		Implemented the GSL Runge Kutta 4-Step method, and removed the requirement
 				for MAX_BODIES to be defined by adding get_lines(), which counts how many
 				lines are in the text file and sets that as the max amount of bodies.

 13-Dec 	1.1.0		Correctly implemented the calc_period() function through use of the bodyvar() struct
 				Code now performs all required tasks through both methods of calling.

 14-Dec		1.1.1		Quality of Life checks and tweaks, implementing the feedback from Homework 1.

 TO-DO				- Create all submission files, and document them within the below text.
*/
/*
 * This code, Nbody.c, is designed to take an input data file and read in
 * starting conditions of N solar bodies[i]s, where N is the number of lines in the 
 * data-file, then use one of two user selected methods to solve the equations
 * of motion for the N-Body system, then plot its motions using Gnuplot.
 * After compiling with:
 *	$ gcc -Wall Nbody.c -lm -lgsl -lgslcblas -o Nbody
 * It is called as:
 * 	$ ./Nbody -[METHOD] [datafile] [size of timestep in days] [length of simulation time in years]
 * so 
 * 	$ ./Nbody -v Horizons.txt 1 30
 * will use the Velocity Verlet method to calculate the motions of bodies as described
 * in the Horizons.txt file, with a step size of 1 day and a total simulation time of 30
 * years (or just over 1 full orbit of Saturn).
 * The different methods that can be used are a written implementation of the Velocity-Verlet
 * method (see: https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet) which is
 * called using the '-v' argument, or the GNU Scientific Library implementation of the 4-Step
 * Runge Kutta method (see: https://www.gnu.org/software/gsl/doc/html/ode-initval.html for GSL,
 * and https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods for Runge Kutta), called using
 * the '-r' argument.
 *
 *
 * Also provided alongside the submission is 
 */
/*
==============================
======== HEADER FILES ========
==============================
*/
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <limits.h> //used for calling a gnuplot script.
#include <gsl/gsl_odeiv2.h> //used to import the gsl rk4 method
#include <gsl/gsl_errno.h> //used for error handling in gsl
#include <gsl/gsl_const_mksa.h> //importing value for gravitational constant
/*
=============================
======== DEFINITIONS ========
=============================
*/
#define MAX_FILE_LINE_SIZE 				  255
#define ITEMS_PER_LINE		  			    8 //Name, Mass, X, Y, Z, Vx, Vy, Vz
#define MAX_NAME_SIZE					  255
#define DIM 						    6 //Number of dimensions for the system of differential equations (default: 6)
#define SOFTENING					  4e8
#define PLOTTING 				    "gnuplot"
#define ENERGY_PLOT			      "energy.script"
#define ERROR_BOUND 					 1e-8
#define G		GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT //Change to 1 when plotting figure8 3 body solution
#define AU 				       1.495978707e11 //Change to 1 to not scale x, y, z values when writing to file
#define DAYS_TO_SECS 					86400 //Change to 1 to input timestep in seconds
#define YEAR 					       365.25 //Change to 1 to input max_time in years
#define ORBIT_PLOT				"3Body.script" //Change this to plot with a different script


//Used to handle error-management.
enum ERROR_CODES {
	SUCCESS = 0,
	ALLOC_FAILED = -1,
	FILE_OPEN_FAILED = -2,
	TIMESTEP_NOT_FOUND = -3,
	BAD_INPUT = -4
};

typedef enum _config_error error_t;

struct _errordesc {
	int code;
	char *message;
} errordesc[] = {
	{ SUCCESS, "No error" },
	{ ALLOC_FAILED, "Memory Allocation Failed" },
	{ FILE_OPEN_FAILED, "File could not be opened" },
	{ TIMESTEP_NOT_FOUND, "No inputted timestep found" },
	{ BAD_INPUT, "Timestep was greater than maximum time"}
};

//used to allow looping over each direction
typedef enum coords { X, Y, Z, N_COORDS } Coords;

//used to store information about each body in the system
typedef struct body {
	char name[MAX_NAME_SIZE];
	double mass;
	double r[N_COORDS]; //Displacement
	double v[N_COORDS]; //Velocity
	double a[N_COORDS]; //Acceleration
} Body;
//used to store information required to calculate the period, if required.
typedef struct period_variables {
	double orbit_angle;
	double initial_angle;
	double last_angle;
	double last_t;
	double half_step;
	double running_total;
	double orbitN;
	double average_period;
} PVar;
/*
=======================
======== FLAGS ========
=======================
*/
static int energy_flag = 0; //set to 1 if you wish to create a "body_energy.dat" data file
static int period_flag = 0; //set to 1 if you wish to calculate the orbital periods of the bodies

/*
===========================
======== FUNCTIONS ========
===========================
*/
/*
 * xcalloc() combines the memory allocation function calloc() with a check that the memory
 * has been allocated correctly. This is wrapped in a function for conciseness.
 */
void *xcalloc(int N, size_t bytes)
{
	void *ret_val = calloc(N, bytes);
	if (ret_val) {
		return ret_val;
	}
	exit (ALLOC_FAILED);
}
/*
 * open_and_check combines opening a file with checking it has opened correctly. This is done
 * numerous times in main and so is wrapped in a function for conciseness.
 */
FILE *open_and_check(char *filename, char *opentype) {
	FILE *f = fopen(filename, opentype);
	if (!f) {
		fprintf(stderr, "Error: Could not open file '%s'.\n", filename);
		exit(FILE_OPEN_FAILED);
	}
	return f;
}
/*
 * get_lines() extracts the number of lines in the inputted datafile
 * This is necessary to know how much memory to allocate to bodies.
 */
static int get_lines(char *filename)
{
	int count = 0;
	char line[MAX_FILE_LINE_SIZE];

	FILE *f = open_and_check(filename, "r");

	while (fgets(line, MAX_FILE_LINE_SIZE, f)) {
		if (line[0] != '#') {
			count += 1;
		}
	}

	fclose(f);
	
	return count;

}
/*
 * get_data() extracts the data from an inputted datafile, and populates
 * the bodies struct with the data.
 */
static void get_data(char *filename, Body *bodies)
{
	char line[MAX_FILE_LINE_SIZE];
	char name_buf[MAX_FILE_LINE_SIZE];
	FILE *input = open_and_check(filename, "r");

	int bodyN = 0;
	while (fgets(line, MAX_FILE_LINE_SIZE, input)) {
		if (line[0] != '#') {
			int nFound = sscanf(line, "%s %lg %lg %lg %lg %lg %lg %lg",
				name_buf, &bodies[bodyN].mass,
				&bodies[bodyN].r[X], &bodies[bodyN].r[Y], &bodies[bodyN].r[Z],
				&bodies[bodyN].v[X], &bodies[bodyN].v[Y], &bodies[bodyN].v[Z]);
			if (nFound == ITEMS_PER_LINE) {
				strncpy(bodies[bodyN].name, name_buf, MAX_NAME_SIZE);
				bodyN++;
			} else {
				fprintf(stderr, "Unknown format: %s\n", line);
			}
		}
	}
	fclose(input);
}
/*
 * calc acceleration is used in the Velocity Verlet method, and calculates the
 * acceleration on each body due to the gravitational attraction of all the
 * other bodies in the system.
 */
static void calc_acceleration(Body *bodies, int Nbodies)
{
	for (int i = 0; i < Nbodies; i++) {
		for (int j=X; j < N_COORDS; j++)
			bodies[i].a[j] = 0;
	}
	for (int i = 0; i < Nbodies; i++) {
		for (int j = 0; j < Nbodies; j++) {
			if (i == j) {
				continue;
			}
			
			double distance = sqrt(pow((bodies[i].r[X] - bodies[j].r[X]), 2) + pow((bodies[i].r[Y] - bodies[j].r[Y]), 2) + pow((bodies[i].r[Z] - bodies[j].r[Z]), 2) + SOFTENING);
			for (int l = X; l < N_COORDS; l++)
				bodies[i].a[l] += -((G * bodies[j].mass) / pow(distance, 3)) * (bodies[i].r[l] - bodies[j].r[l]);
		}
	}
}
/*
 * vel_verl() implements the velocity verlet algorithm to solve the differential
 * equations. This calculates the next position using current velocity and acceleration,
 * then calculates the new acceleration in the new position, before calculating the
 * new velocity using the old velocity, and the average of the old and new acceleration
 * multiplied by the timestep.
 */ 
static void vel_verl(Body *bodies, int Nbodies, double timestep)
{
	for (int i = 0; i < Nbodies; i++) {
		for (int j = X; j < N_COORDS; j++) {
			bodies[i].v[j] = bodies[i].v[j] + bodies[i].a[j] * timestep * 0.5;
			bodies[i].r[j] = bodies[i].r[j] + bodies[i].v[j] * timestep;
		}
	}

	calc_acceleration(bodies, Nbodies);
	
	for (int i = 0; i < Nbodies; i++) {
		for (int j=X; j<N_COORDS; j++) {
			bodies[i].v[j] = bodies[i].v[j] + 0.5 * bodies[i].a[j] * timestep;
		}
	}
}
/*
 * calc_energy() calculates the gravitation potential and kinetic energy at each timestep,
 * then sums them to find the total energy. It outputs the energy of each body into a
 * bodyname_energy.dat file, which can then be plotted using the energy.script gnuplot file.
 * This function is used for verifying the accuracy of the model.  
 */
static void calc_energy(Body *bodies, int Nbodies, int t)
{
	double potential[Nbodies];
	double kinetic[Nbodies];
	double total[Nbodies];

	for (int i=0; i<Nbodies; i++) {
		potential[i] = 0;
		kinetic[i] = 0;
		total[i] = 0;
		for (int j=0; j<Nbodies; j++) {
			if (i==j)
				continue;

			potential[i] += -(G * bodies[i].mass * bodies[j].mass) / sqrt(pow((bodies[i].r[X] - bodies[j].r[X]), 2) + pow((bodies[i].r[Y] - bodies[j].r[Y]), 2) + pow((bodies[i].r[Z] - bodies[j].r[Z]), 2));
		}
		
		kinetic[i] += 0.5 * bodies[i].mass * (pow(bodies[i].v[X], 2) + pow(bodies[i].v[Y], 2) + pow(bodies[i].v[Z], 2));
		
		total[i] = kinetic[i] + potential[i];

		char filename[255];
		strcpy(filename, bodies[i].name);
		strcat(filename, "_energy.dat");
		FILE *f = open_and_check(filename, "a");
		fprintf(f, "%i \t %lg \t %lg \t %lg\n", t, potential[i], kinetic[i], total[i]);
		fclose(f);
	}
}
/*
 * calc_period() calculates the orbital period of each body. It is used for verifying the
 * accuracy of the model, so has been written specifically for the solar system, including
 * the earth's moon. It uses the bodyvar struct to store the initial angle of orbit, then 
 * each time it is called it updates the previous angle and calculates the current angle. If 
 * the condition previous_angle < initial_angle < current_angle; then it is either half way 
 * through its orbit or has completed a full orbit, since -y/-x = y/x it cycles between -pi 
 * and pi radians, rather than 0-2pi. For this reason when it the condition is satisfies it 
 * checks the half-orbit flag. If the flag is 1, then only a half step has been completed, so
 * it sets the flag to 0. The next time the condition is satisfied the flag will now read 0, and 
 * a full-orbit will have been completed, so the running total of all the orbital periods for the
 * body is updated, along with how many orbits has been completed. It then calculates the average 
 * orbital period for that body, which updates with each full orbit. It returns the array of each 
 * bodies average period.
 */
static void calc_period(Body *bodies, double t, int Nbodies, int i, PVar *bodyvar)
{
	if (!strcmp(bodies[i].name, "Moon")) {
		bodyvar[i].orbit_angle = atan((bodies[i].r[Y] - bodies[3].r[Y]) / (bodies[i].r[X] - bodies[3].r[X]));
	} else {
		bodyvar[i].orbit_angle = atan(bodies[i].r[Y] / bodies[i].r[X]);
	}
	if (t==0) {
		bodyvar[i].last_angle = 0;
		bodyvar[i].initial_angle = bodyvar[i].orbit_angle;
		bodyvar[i].orbitN = 0;
		bodyvar[i].half_step = 1;
	}

	if (bodyvar[i].last_angle < bodyvar[i].initial_angle && bodyvar[i].orbit_angle > bodyvar[i].initial_angle) {

		if (bodyvar[i].half_step == 0) {
			if (bodyvar[i].orbitN == 0) {
				bodyvar[i].last_t = t;
				bodyvar[i].running_total = t;
			} else {
				bodyvar[i].running_total += t - bodyvar[i].last_t;
				bodyvar[i].last_t = t;

			}
			bodyvar[i].orbitN++;
			bodyvar[i].average_period = bodyvar[i].running_total / (bodyvar[i].orbitN);
			bodyvar[i].half_step = 1;
		} else if (bodyvar[i].half_step == 1) {
			bodyvar[i].half_step = 0;
		}

	}
	bodyvar[i].last_angle = bodyvar[i].orbit_angle;
}
/*
 * func() is used by the GSL ODE algorithms, and is used to set up the system
 * of differential equations to be solved. For the Nbody system we need a system
 * of DIM * Nbodies differential equations, as each body requires DIM differential
 * equations described its position in 3 directions, and velocity in three
 * directions. The GSL documentation page for its ODEIV2.h package can be found
 * here: https://www.gnu.org/software/gsl/doc/html/ode-initval.html
 */
int func (double t, const double y[], double f[], void *params)
{
	(void) t;
	double *mass = (double *)params;
	int Nbodies = mass[0];

	for (int i = 0; i < Nbodies; i++) {
		f[DIM * i + 0] = y[DIM * i + 3];
		f[DIM * i + 1] = y[DIM * i + 4];
		f[DIM * i + 2] = y[DIM * i + 5];

		double x_acc = 0;
		double y_acc = 0;
		double z_acc = 0;

		for (int j = 0; j < Nbodies; j++) {
			if (i != j) {
				double x_diff = y[DIM * i + 0] - y[DIM * j + 0];
				double y_diff = y[DIM * i + 1] - y[DIM * j + 1];
				double z_diff = y[DIM * i + 2] - y[DIM * j + 2];

				double distance = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
				
				double pre = (G * mass[j+1] / pow(distance, 3));

				x_acc -= pre * (x_diff);
				y_acc -= pre * (y_diff);
				z_acc -= pre * (z_diff);
			}
		f[DIM * i + 3] = x_acc;
		f[DIM * i + 4] = y_acc;
		f[DIM * i + 5] = z_acc;
		}
	}
	return GSL_SUCCESS;
}
/*
 * gnu_plot() is used to call a system command to run an inputted script through
 * gnuplot. In main, the scripts used are #define'd at the top of the document.
 */
void gnu_plot (char *plotwith, char *script)
{
	char command[PATH_MAX];
	snprintf(command, sizeof(command), "%s %s", plotwith, script);
	system(command);
}

/*
 * update_bodies() is used to populate the bodies struct with new positions
 * and velocities from the y array outputted by the gsl rk4 algorithm. This
 * is performed in a function as it is used numerous times in main.
 */
void update_bodies(Body *bodies, int Nbodies, double *y)
{
	for (int i = 0; i < Nbodies; i++) {
		for (int l = X; l < N_COORDS; l++) {
			bodies[i].r[l] = y[DIM * i + l];
			bodies[i].v[l] = y[DIM * i + (3 + l)];
		}
	}
}
/*
======================
======== MAIN ========
======================
 */
int main (int argc, char **argv)
{
	int c, Nbodies;
	double timestep, max_time;

	printf("#REVISION DATE: %s\n#VERSION: %s\n", REV_DATE, VERSION);
	
	while (1)
	{
		static struct option long_options[] =
		{
			{"Velocity-Verlet",     required_argument, 	0, 'v'},
			{"Runge-Kutta4",  	required_argument, 	0, 'r'},
			{0, 0, 0, 0}
		};
		// getopt_long stores the option index here.
		int option_index = 0;

		c = getopt_long (argc, argv, "v:r:",
			long_options, &option_index);

		if (c == -1)
			break;

		Nbodies = get_lines(optarg);

		Body *bodies = xcalloc(Nbodies, sizeof(Body));
		PVar *bodyvar = xcalloc(Nbodies, sizeof(PVar));

		get_data(optarg, bodies);

		double expected_period[] = {0, 87.97, 224.65, 365.256, 687.0, 4331.98, 10776.3, 30688.5, 60182, 27.321}; //The observed orbital periods of each celestial body in the Horizons.txt input file.
		
		if (energy_flag) {
			
			for (int i = 0; i < Nbodies; i++) {
				char filename[MAX_NAME_SIZE];
				strcpy(filename, bodies[i].name);
				strcat(filename, "_energy.dat");
				FILE *f = open_and_check(filename, "w+");
				fprintf(f, "#Time \t PotentialEnergy \t KineticEnergy \t Total Energy\n");
				fclose(f);
			}
		}

		FILE *outfile = open_and_check("OrbitalData.dat", "w+");
		fprintf(outfile, "#Time\t");
		for (int i=0; i<Nbodies; i++) {
			fprintf(outfile, "%sX\t%sY\t%sZ\t", bodies[i].name, bodies[i].name, bodies[i].name);
		}
		fprintf(outfile, "\n");

		if (optind < argc ) {
			timestep = atof(argv[optind]) * DAYS_TO_SECS;
			optind++;	
		} else {
			printf("Please input a timestep.\n");
			exit(TIMESTEP_NOT_FOUND);
		}
		if (optind < argc) {
			max_time = atof(argv[optind]) * YEAR * DAYS_TO_SECS;
			optind++;
		}

		if (timestep > max_time / YEAR) {
			printf("The timestep must be less than the max time elapsed\n");
			exit(BAD_INPUT);
		}

		double average_period[Nbodies];

		switch (c) {
		case 0:
			if (long_options[option_index].flag != 0)
				break;
			printf ("option %s", long_options[option_index].name);
			if (optarg)
				printf (" with arg %s", optarg);
			printf ("\n");
			break;

		case 'v':
			printf("# Velocity-Verlet method on %s\n", optarg);
			printf("# Plotting every %lg days for %lg years.\n", timestep / DAYS_TO_SECS, max_time / DAYS_TO_SECS / YEAR);
			calc_acceleration(bodies, Nbodies);

			for (double t = 0; t < max_time; t += timestep) {
				vel_verl(bodies, Nbodies, timestep);
				if (energy_flag) {
					calc_energy(bodies, Nbodies, t);
				}
				for (int i = 0; i < Nbodies; i++) {
					if (i == 0 ) {
						fprintf(outfile, "%lg\t", t / DAYS_TO_SECS);
					}
					
					fprintf(outfile, "%lg \t %lg \t %lg \t", bodies[i].r[X] / AU,
						bodies[i].r[Y] / AU, bodies[i].r[Z] / AU);
					
					if (i == Nbodies - 1) {
						fprintf(outfile, "\n");
					}
					if (period_flag) {
						calc_period(bodies, t / DAYS_TO_SECS, Nbodies, i, bodyvar);
						average_period[i] = bodyvar[i].average_period;
						if (t >= max_time - timestep) {
							if (strcmp(bodies[i].name, "Sun")) {
								if (average_period[i] > 1) {
									double percent_error = ((expected_period[i] - average_period[i]) / expected_period[i]) * 100.0;
									printf("%s average period: %lg, expected period: %lg. Accuracy: %lg%% \n", bodies[i].name, average_period[i], expected_period[i], percent_error);
								} else {
									printf("%s did not complete an orbit.\n", bodies[i].name);
								}
							}
						}
					}
				}
			}	
			fclose(outfile);
			free(bodies);
			free(bodyvar);
			break;

		case 'r':
		{
			printf("# 4-Step Runge Kutta method on %s\n", optarg);
			printf("# Plotting every %lg days for %lg years.\n", timestep / DAYS_TO_SECS, max_time / DAYS_TO_SECS / YEAR);
			
			double masses[Nbodies + 1];
			masses[0] = Nbodies; //Nbodies is set to the first element of masses,
					     //so it can be passes into func, since sizeof(array)/sizeof(array[0])
					     //doesn't work in functions.
			
			for (int i = 1; i <= Nbodies; i++) {
				masses[i] = bodies[i - 1].mass;
			}
			
			gsl_odeiv2_system sys = {func, NULL, DIM * Nbodies, masses}; //func describes the system, NULL because we don't need a Jacobian matrix,
										 //DIM*Nbodies differential equations, as there are N equations for (x, y, z), (Vx, Vy, Vz)
										 //masses contains Nbodies as first element, then the masses of all the bodies. 

			double y[DIM * Nbodies];
			
			for (int i = 0; i < Nbodies; i++) {
				for (int l = X; l < N_COORDS; l++) {
					y[i * DIM + l] = bodies[i].r[l];
					y[i * DIM + (3 + l)] = bodies[i].v[l];

				}
			}

			gsl_odeiv2_driver *d = 
				gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,
								timestep, ERROR_BOUND, ERROR_BOUND);

			double t = 0.0;
			for (double j = 0; j <= max_time; j+= timestep) {
				double ti = t + timestep;
				int s = gsl_odeiv2_driver_apply(d, &t, ti, y);
				if (s != GSL_SUCCESS) {
					printf("Error: driver returned %d\n", s);
					break;
				}

				if (energy_flag) {
					update_bodies(bodies, Nbodies, y);
					calc_energy(bodies, Nbodies, j / DAYS_TO_SECS);
				}

				if (period_flag) {

					update_bodies(bodies, Nbodies, y);

					for (int i = 0; i < Nbodies; i++) {
						calc_period(bodies, j / DAYS_TO_SECS, Nbodies, i, bodyvar);
						average_period[i] = bodyvar[i].average_period;
						if (j >= max_time - timestep) {
							if (strcmp(bodies[i].name, "Sun")) {
								if (average_period[i] > 1) {
									double percent_error = ((expected_period[i] - average_period[i]) / expected_period[i]) * 100.0;
									printf("%s average period: %lg, expected period: %lg. Accuracy: %lg%% \n", bodies[i].name, average_period[i], expected_period[i], percent_error);
								} else {
									printf("%s did not complete an orbit.\n", bodies[i].name);
								}
							}
						}
					}
				}

				fprintf(outfile, "%lg\t", j / DAYS_TO_SECS);
				for (int i = 0; i < Nbodies; i++) {
					fprintf(outfile, "%lg\t%lg\t%lg\t", y[DIM * i + 0] / AU, y[ DIM * i + 1] / AU, y[DIM * i + 2] / AU);
				}
				fprintf(outfile, "\n");
			}
			fclose(outfile);
			gsl_odeiv2_driver_free(d);
			free(bodies);
			free(bodyvar);
			break;
		}

		case '?':
		/* getopt_long already printed an error message. */
			break;

		default:
			abort ();
		}
	}
	if (optind < argc)
	{
		printf ("non-option ARGV-elements: ");
		while (optind < argc)
			printf ("%s ", argv[optind++]);
		putchar ('\n');
	}
	
	gnu_plot(PLOTTING, ORBIT_PLOT);

	if (energy_flag) {
		gnu_plot(PLOTTING, ENERGY_PLOT);
	}

	exit(SUCCESS);
}
