/*
 Title:		PHYM004 Homework 1
 Author:	Tom Joshi-Cale tj294 670000715
 Style:		Linux Kernel - https://www.kernel.org/doc/html/v4.10/process/coding-style.html
*/
static const char *VERSION = "1.0.2";
static const char *REV_DATE = "26-Oct-2020";
/*
 Date		Version 	Comments
 ----		-------		--------
 13-Oct		0.1.0		Successfully wrote a function to read from input
				file the dimensions of the input matrix. Unable to
				extract matrix data and pass it to main. Will need
				to use memory allocation.
 
 14-Oct 	0.1.1		Successfully managed to extract the matrix and pass
				it to main. Managed to successfully calculate the 
				Frobenius norm of input matrices of any size.
				Task One Complete!
 
 15-Oct		0.1.2		Changed the memory allocation of Create2D to map
				onto a 1D array as suggested in lecture P-04a. This
				stores the elements in contiguous memory, making the
				elements faster to access.
 
 20-Oct		0.1.3		Re-wrote the case(f) and do_frobenius function to 
				match the output style of mat_gen.c.
				Successfully read two input arguments into the
				matrix product case.
				Tested for square matrix in the Determinant case.
 
 20-Oct 	0.2.0		Managed to successfully calculate and print the
				Transpose matrix.
				Task 2 Complete!

 23-Oct 	0.3.0		Managed to successfully calculate and print the
				results of matrix multiplication between two
				n-dimensional matrices.

 26-Oct 	0.4.0		Managed to successfully calculate the deteminant of
				an N-Dimensional matrix, using the method laid out
				by Khamsi, M. A. at Determinant and inverse of matrices
				http://www.sosmath.com/matrix/inverse/inverse.html (1999)
				This involved the creation of a cofactor function which
				will be useful for the remaining tasks.

 26-Oct 	0.5.0		Managed to successfully calculate the adjoint of an
				N-Dimensional matrix: The adjoint is the inverse of
				the matrix of cofactors, and I already had functions
				to calculate a matrix of cofactors and inverse a
				matrix, so this was as simple as calling those 
				functions in order.

 26-Oct 	0.5.1		Changed the if(i%2=0), multiply cofactor element by 1
				else multiply by -1 to one line of code:
				  cof_mat[i] = do_determinant(minors, (dim-1)) * pow(-1, i%2);
				as it is more concise.

 26-Oct 	0.5.2		The version of the multiplying by -1 only works for
				NxN matrices where N is odd. Fixed this by replacing
				the [*pow(-1, i%2)] to [*pow(-1, el_row+el_col)]

 26-Oct 	1.0.0		Successfully calculated the inverse of an NxN matrix,
				and tested with 2x2, 3x3, 4x4 and 5x5 matrices. All
				6 tasks are now complete, so code moves to 1.0.0

 26-Oct 	1.0.1		Added a print_out() function, to print the output
				of each task to the command line. This is because
				some tasks were doing it in functions, others were
				doing it in main and I wanted it to be standard,
				which would either require it occuring in main for
				each task (6 times) or writing it into a function
				(written out once). Adding it as a function is clearly
				the better choice.
				Current Bugs: When doing -d, -a or -i for matrices of
				8x8 or greater, the compiler kills the program. This
				is apparently due to it running out of memory.

 26-Oct 	1.0.2		Managed to sort the determinant for large matrices:
				as the det(A) = sum_j^n a_{i, j}*A_{i, j} for ANY
				FIXED i OR j, there is no need to create a full NxN
				matrix of cofactors when calculating the determinant,
				only a Nx1 matrix. This is much faster, and now the
				determinant of 10x10 matrices can be computed in <1s.
				This should also help with computing the adjoint (and
				therefore inverse) matrices, as finding the determinants
				for the cofactor matrix is now much faster and easier.

 27-Oct 	1.0.3		Added error handling to read_matrix() function, so that if
 				read matrix does not match the dimensions given in the file
 				the user will be warned. Also if one of the elements of
 				the inputted matrix is Not a Number, the user will be
 				told.
 
 28-Oct		1.0.4		Moved the print_header(), get_dimensions() and read_matrix()
 				functions out of each case so that it is only written in
 				main once, to ensure my code holds to the DRY principle.
*/
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h> 	//For parsing the command line
#include <string.h> 	//For parsing the input file
#include <math.h>	// For performing matrix calculations - Compile with -lm flag
#include <errno.h>
/*
 * This code, mat_test.c, is designed to take one or more input matrices from
 * the mat_gen.c code created by C.D.H. Williams (2019) and perform various
 * actions to those matrices. The actions it should be able to perform are:
 *	1. Calculate Frobenius norm of a single matrix
 *	2. Calculate the Transpose of a single matrix
 *	3. Calculate the Matrix Product of two inputted matrices
 *	4. Calculate the determinant of a single matrix
 *	5. Calculate the adjoint of a single matrix
 *	6. Calculate the inverse of a single matrix
 * it is designed to be invoked from the command line as:
 * 	
 * ./mat_test -[flag] [input_matrix.txt] ([2nd_input_matrix.txt] - optional)
 * 
 * If required, then an output file will be created in the same format as the 
 * output of mat_gen.c in the current directory:
 *
 *	# ./mat_test -[flag] [input_matrix.txt] ([2nd_input_matrix.txt])
 *	# Version = 1.0.0, Revision date = 29-Oct-2020
 *	matrix NROWS NCOLS
 *	<DATA>
 *	end
 *
 * where NROWS and NCOLS represent the dimensions of the data, '<DATA>' represents
 * the tab-separated values of the output matrix.
 */
#define MAX_FILE_LINE_SIZE 255
/*
 * =============================================================================
 * Functions:
 * =============================================================================
 */
/*
 * The print_header() function prints the two comment lines required by the
 * output format, which consist of the command line input used to run the code
 * and the code version and latest revision date of the code version run.
 */
static void print_header(char **argv, int argc)
{
	printf("#");
	for (int i = 0; i < argc; i++)
		printf(" %s", argv[i]);
	printf("\n");
	printf("# Version = %s, Revision date = %s\n", VERSION, REV_DATE);
}
/* 
 * The dimension_get function takes in the filename of the matrix being operated on
 * and reads out of that file the number of rows and number of columns of the input
 * matrix, which are saved to 2 inputted pointers.
 */
static void dimension_get(char *filename, int *rows, int *cols)
{

	char *status;
	char line[MAX_FILE_LINE_SIZE];

	FILE *inf = fopen( filename, "r");
	if (!inf) {
		fprintf(stderr, "Error, could not open file '%s'\n", filename);
		exit(1);
	}

	status = fgets(line, sizeof(line), inf);
	while (status) {
		status = fgets(line, sizeof(line), inf);
		if (status[0] == 'm') {
			sscanf(line, "matrix %i %i", rows, cols);
			break;
		}
	}
	fclose(inf);
}
/* 
 * The Create2D function dynamically allocates memory to story an inputted files data
 * to, and returns it out of the function. When the array has finished being used, 
 * it is important to free up the memory using free().
 * Thanks for this function given to user Ryyker on StackOverflow.com, who assisted
 * in solving my memory allocation problem.
 */
static double *create2D(int cols, int rows)
{
	double *arr = calloc(rows * cols, sizeof(double)+1);	
	return arr;
}
/*
 * The read_matrix function parses the data from the inputted text file, converts
 * it into a 2D array, and then populates the dynamically allocated 2D array 
 * generated by Create2D with the data.
 */
static double *read_matrix(char *filename, int rows, int cols)
{
	char *status;
	char *endptr;
	char line[MAX_FILE_LINE_SIZE];
	int i=0;
	char *data_line;
	double *array = create2D(rows, cols);

	FILE *inf = fopen(filename, "r");
	if (!inf) {
		fprintf(stderr, "Error, could not open file '%s'\n", filename);
		exit(1);
	}
	while ((status = fgets(line, sizeof(line), inf)) != NULL) {
		if (status[0] == 'e') { //finish searching when line is 'end'
			break;
		} else if (status[0] == '#' || status[0] == 'm') { //skip over comments (#) and matrix size (m)
			continue;
		} else { //read in matrix data
			int j = 0;
			data_line = strtok(line, " \t");
			errno = 0;
			while (data_line != NULL) {
				array[cols*i + j] = strtod(data_line, &endptr);
				//If the data read in is not a number, return an error.
				if (array[cols*i + j] == 0 && (errno != 0 || endptr == data_line) && strcmp("\n", data_line)) {
					fprintf(stderr, "Error: Inputted matrix contains non-number elements on line %i\n", i+1);
					exit(1);
				}
				data_line = strtok(NULL, " \t");
				j++;
			}
			if (j-1 != cols) { //checking that the dimensions of matrix in file matches dimensions in file
				fprintf(stderr, "Error: Line %i of inputted matrix does not match inputted dimensions.\n", i+1);
				exit(1);
			}
			i += 1;
		}
	}
	fclose(inf);
	return array;
}
/*
 * The do_frobenius function finds the Frobenius norm of the inputted matrix using
 * the method outlined by:
 * Weisstein, Eric W. "Frobenius Norm." From MathWorld--A Wolfram Web Resource. https://mathworld.wolfram.com/FrobeniusNorm.html
 * The Frobenius Norm is defined as the square root of the sum of absolute 
 * squares of the elements of a matrix.
 */ 
static double *do_frobenius(char *filename, double *matrix, int rows, int cols)
{		
	double *ans= create2D(1, 1);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++)
			ans[0] = ans[0] + matrix[cols * i + j] * matrix[cols * i + j];
	}
	ans[0] = sqrt(ans[0]);
	return ans;
	free(ans);
}
/*
 * The do_transpose function creates an output matrix with the opposite dimensions to
 * the input matrix [i.e. an input of 3x4 will create an output of 4x3], and then
 * reads input_mat[i,j] into output_mat[j,i]. Since the matrices have been mapped
 * onto a 1D array, the transposition takes place as output[rows*j+i] = input[cols*i+j]
 */
static double *do_transpose(double *inpmat, int rows, int cols) 
{
	double *outmat = create2D(rows, cols);
	
	if (outmat) {
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++)
				outmat[rows * j + i] = inpmat[cols * i + j];
		}
	}
	return outmat;
	free(outmat);
}
/*
 * The do_matprod performs matrix multiplication between two matrices.
 */
static double *do_matprod(double *mat1, double *mat2, int rows1, int rows2, int cols1, int cols2)
{

	if (cols1 != rows2) {
		fprintf(stderr, "Error: Matrix 1 must have number of columns equal to Matrix 2's number of rows.\n");
		exit(1);
	}

	double *outmat = create2D(cols2, rows1);
	
	for (int i = 0; i < rows1; i++) {
		for (int j = 0; j < cols2; j++) {
			for (int n = 0; n < cols1; n++)
				outmat[cols2 * i + j] += mat1[cols1 * i + n] * mat2[cols2 * n + j];
		}
	}
	return outmat;
	free(outmat);
}
/*
 * The det_cofactor() function both references and is referenced by the do_determinant
 * function, so it is initialised before do_determinant, but written out afterwards.
 */
static double *det_cofactor(double *matrix, int dim);
/*
 * The do_determinant function finds the determinant of an NxN matrix. For a 1x1 it
 * just returns the number in the matrix. For a 2x2 it performs the simplistic a*d-b*c
 * calculation, and for all other dimensions it creates a matrix of cofactors that are
 * required in the N-dimensional calculation from M. A. Khamsi.
 */
static double *do_determinant(double *matrix, int rows)
{
	int dim = rows;
	double det = 0;

	double *cof = calloc(dim, sizeof(double));
	double *ans = create2D(1, 1);
	
	if (dim == 1) {
		det = matrix[0];
	} else if (dim == 2) {
		det = matrix[0] * matrix[3] - matrix[1] * matrix[2];
	} else {
		cof = det_cofactor(matrix, dim);
		for (int j = 0; j < dim; j++)
			det = det + (cof[j] * matrix[j]);
	}
	free(cof);
	ans[0]=det;
	return ans;
	free(ans);
}
/*
 * The det_cofactor() matrix creates an array of cofactors for the 0th row of an
 * inputted matrix. This Nx1 array of cofactors is all that is required for
 * calculating the determinant of the inputted matrix, so this function is much
 * faster than creating a NxN array of cofactors and only using one row.
 */
static double *det_cofactor(double *matrix, int dim) {
	int el_col = 1;
	int el_row = 1;
	int l = 0;

	double *cof_mat = calloc(dim, dim*sizeof(double));

	for (int i = 0; i < dim; i++) {
		double *minors = create2D(dim-1, dim-1);
		for (int j = 0; j < dim * dim; j++) {
			if ((j < el_row * dim && j >= (el_row-1) * dim) || j % dim == el_col - 1)
				continue;
			minors[l] = matrix[j];
			l += 1;
		}
		cof_mat[i] = do_determinant(minors, (dim - 1))[0] * pow(-1, el_col + 1);
		free(minors);
		el_col += 1;
		l = 0;
	}
	return cof_mat;
	free(cof_mat);
}
/*
 * The cofactor() function calculates the cofactor matrix for an inputted n-dimensional
 * matrix. This involves for each element of the matrix, calculating the resulting
 * matrix if the row and column of the element is removed, and passing that matrix back
 * to the do_determinant() function. It must the multiply this cofactor for each element
 * by -1 ^ {row + col}. This function is only called by the do_adjoint() function, but
 * is a separate function from that one as finding the matrix of cofactors for a given
 * matrix is a common task, so this function is separate in the interests of adaptability
 * and modularity.
 */
static double *cofactor(double *matrix, int dim)
{	
	int el_row = 1;
	int el_col = 1;
	int l = 0;
	int n_elements = dim*dim;

	double *cof_mat = create2D(dim, dim);

	for (int i = 0; i < n_elements; i++) {
		double *minors = create2D(dim - 1, dim - 1);
		if (el_col == dim + 1) {
			el_col = 1;
			el_row += 1;
		}
		for (int j = 0; j < n_elements; j++) {
			if ((j < el_row * dim && j >= (el_row - 1) * dim) || j % dim == el_col - 1)
				continue;
			minors[l] = matrix[j];
			l += 1;
		}
		cof_mat[i] = do_determinant(minors, (dim - 1))[0] * pow(-1, (el_row + el_col));	
		free(minors);
		el_col += 1;
		l = 0;
	}

	return cof_mat;
	free(cof_mat);
}
/*
 * The do_adjoint() function simply calculates the matrix of cofactors using
 * the cofactor() function, and then uses the do_transpose() function turn it
 * into the adjoint.
 */
static double *do_adjoint(double *matrix, int rows)
{
	double *cof = create2D(rows, rows);
	if (rows != 1) {
		cof = cofactor(matrix, rows);
	} else {
		cof[0] = 1;
	}
	return do_transpose(cof, rows, rows);
}
/*
 * The do_inverse() function multiplies the adjoint matrix by the inverse
 * of the determinant.
 */
static double *do_inverse(double *matrix, int dim, double det)
{
	if (det == 0) {
		fprintf(stderr, "Matrix is singular (det = 0)\nCannot calculate inverse.\n");
		exit(1);
	} else {
		double *adj = create2D(dim, dim);

		adj = do_adjoint(matrix, dim);

		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++)
				adj[dim * i + j] = (1 / det) * adj[dim * i + j];
		}
		return adj;
	}
}
/*
 * The print_out() function allows each use-case of mat_test.c to print its
 * output in a standard format, which would be readable by mat_test.c. For
 * outputs which would be a single number, it outputs them as a 1x1 matrix.
 * if the user wishes to print these outputs into a file, the code can be
 * piped into a txt file as so:
 * 	$ ./mat_test -m matrix1.txt matrix2.txt > output.txt 
 * which will  print the typical output into a text file named output.txt
 * The outputted precision has been chosen as 12 decimal places to match
 * the precision of mat_gen.c by C.D.H. Williams.
 */
static void print_out(double *output, int rows, int cols)
{
	printf("matrix %i %i\n", rows, cols);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++)
			printf("%.12g\t", output[rows * i + j]);
		printf("\n");
	}
	printf("end\n");
}
/*
 * =============================================================================
 * Main:
 * =============================================================================
 */
int main (int argc, char **argv)
{
	int c, rows1, cols1, rows2, cols2;
	double *ans = NULL;
	double *mat1 = NULL;
	double *mat2 = NULL;
	double *mat3 = NULL;

	rows1 = 0;
	cols1 = 0;
	rows2 = 0;
	cols2 = 0;

	while (1)
	{
		static struct option long_options[] = { 
			{"froeb",	required_argument, 0, 'f'},
			{"trans", 	required_argument, 0, 't'},
			{"mprod",	required_argument, 0, 'm'},
			{"det",		required_argument, 0, 'd'},
			{"adj",		required_argument, 0, 'a'},
			{"inv",		required_argument, 0, 'i'},
			{0, 0, 0, 0}
		};
		// getopt_long stores the option index here. 
		int option_index = 0;

		c = getopt_long (argc, argv, "f:t:m:d:a:i:",
			long_options, &option_index);

		if (c == -1)
			break;

		print_header(argv, argc);
		dimension_get(optarg, &rows1, &cols1);
		mat1 = read_matrix(optarg, rows1, cols1);

		switch (c) {
		case 0:
			if (long_options[option_index].flag != 0)
				break;
			printf("option %s", long_options[option_index].name);
			if (optarg)
				printf(" with arg %s", optarg);
			printf("\n");
			break;

		case 'f': //Calculate Frobenius Norm - COMPLETE 14-OCT 4:43PM 
			if (mat1) {
				ans = do_frobenius(optarg, mat1, rows1, cols1);
				free(mat1);
			}
			print_out(ans, 1, 1);
			break;

		case 't': //Calculate Transpose - COMPLETE 20-OCT 7:50PM
			if (mat1) {
				mat2 = do_transpose(mat1, rows1, cols1);
				free(mat1);
			}
			print_out(mat2, rows1, cols1);
			break;

		case 'm': //Calculate Matrix Product - COMPLETE 23-OCT 6:32PM
			if (argv[optind] == NULL) {
				fprintf(stderr, "Error: Please input two matrices.\n");
				exit(1);
			}
			dimension_get(argv[optind], &rows2, &cols2);
			mat2 = read_matrix(argv[optind], rows2, cols2);
			if (mat1 && mat2) {
				mat3 = do_matprod(mat1, mat2, rows1, rows2, cols1, cols2);
				free(mat1);
				free(mat2);
			}
			print_out(mat3, rows1, cols2);
			optind = optind + 1;
			break;


		case 'd': //Calculate Determinant - COMPLETE 26-OCT 11:09AM
			if (rows1 != cols1) {
				fprintf(stderr, "Error: %s is not a square matrix. Please input a square matrix.\n", optarg);
				exit(1);
			} else {
				if (mat1) {
					ans = do_determinant(mat1, rows1);
					free(mat1);	
				}
			}
			print_out(ans, 1, 1);
			break;

		case 'a': //Calculate Adjoint - COMPLETE 26-OCT 11:32AM
			if (rows1 != cols1) {
				fprintf(stderr, "Error: %s is not a square matrix. Please input a square matrix.\n", optarg);
				exit(1);
			} else {
				if (mat1) {
					mat2 = do_adjoint(mat1, rows1);
					free(mat1);
				}
				print_out(mat2, rows1, cols1);
			}
			break;

		case 'i': //Calculate Inverse - COMPLETE 26-OCT 12:47PM
			if (rows1 != cols1) {
				fprintf(stderr, "Error: %s is not a square matrix. Please input a square matrix.\n", optarg);
				exit(1);
			} else {
				ans = do_determinant(mat1, rows1);
				mat3 = do_inverse(mat1, rows1, ans[0]);
				print_out(mat3, rows1, cols1);
			}
			break;

		case '?':
			break;

		default:
			abort ();
		}
	}
	if (optind < argc) {
		printf("non-option ARGV-elements: ");
		while (optind < argc)
			printf("%s ", argv[optind++]);
		putchar ('\n');
	}
	return 1;
}
