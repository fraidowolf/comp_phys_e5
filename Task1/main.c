// this include is needed, remember to add the gsl library when you compile as in E1.
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#define PI 3.14159265359
#define a0 0.529177249 // bohr radius a0

#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>


void print_mat(gsl_matrix *matrix, size_t len){
    size_t i, j;
    double element;

    for (i = 0; i < len; ++i) {
        for (j = 0; j < len; ++j) {
            element = gsl_matrix_get(matrix, i, j);
            printf("%f ", element);
        }
        printf("\n");
    }
}

void init_mat(gsl_matrix *matrix,double h, size_t len){
	size_t i, j;
	double value = 0;
	for (i = 0; i < len; ++i) {
		for (j = 0; j < len; ++j) {
			if (i==j){
				value = -2.0/(h*h);
			}else if ((i-1==j) | (j-1==i)){
				value = 1.0/(h*h);
			}else{
				value = 0;
			}

			gsl_matrix_set(matrix, i, j, value);
		}
	}
}

double rho(double x){
	return (1.0/(PI*a0*a0*a0))*exp(-2.0*x/a0);
}

void print_vec(gsl_vector *vec, double a, double h, size_t len){
	FILE *file = fopen("data.csv", "w");
	for(int i =0; i<(double)len; i++){
		double value = gsl_vector_get(vec,i);
		//printf("%e\n",value);
		fprintf(file,"%f\t%f\n", a + h*(double) i, value/ (a + (double) i*h));
	}
	fclose(file);
}

void init_f(gsl_vector *f, double a, double h,  size_t len){

	for(size_t i=0;i<len;i++){
		double x = a + h * (double) i;
		gsl_vector_set(f,i,-4*PI*x*rho(x));
	}
	double x = a + h * (double) (len-1);
	gsl_vector_set(f,len-1,-4*PI*x*rho(x)-1/(h*h));
}
  

void solve_matrix_eq(gsl_matrix *matrix, gsl_vector *b, gsl_vector *x,size_t len){
	gsl_permutation *p = gsl_permutation_alloc(len);
	int s;

	gsl_linalg_LU_decomp(matrix, p, &s);
	gsl_linalg_LU_solve(matrix,p,b,x);
	
	gsl_permutation_free(p);
	}
			


int main() {
	
	const size_t steps = 1000;
	gsl_vector *U_gsl=gsl_vector_alloc(steps);
	gsl_vector *f_gsl=gsl_vector_alloc(steps);
	double a = 0.01;   // r_0
	double b = 20; // r_n (Å)
	double h = (b-a)/(double) steps;
	gsl_matrix *mat = gsl_matrix_alloc(steps, steps);

	init_mat(mat,h,steps);
	init_f(f_gsl,a,h,steps);

	solve_matrix_eq(mat,f_gsl,U_gsl,steps);

	//print_mat(mat,steps);
	print_vec(U_gsl,a,h,steps);	
	gsl_vector_free(U_gsl);
	gsl_vector_free(f_gsl);
	gsl_matrix_free(mat);	

	return 0;
}
