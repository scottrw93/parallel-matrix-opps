#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define MAXTHRDS 124

void initMat(int n, double mat[]){
  	int i, j;
	for (i=0;i<n;i++){
    	for (j=0;j<n;j++){
      		mat[i*n+j] = i+j;
      	}
    }
}


void printMat(int n, double mat[] , char* matName){
	int i, j;
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) {
			printf("%s[%d][%d]=%.1f ", matName, i, j, mat[i*n+j]);
		}
		printf("\n");
	}
}

int main(int argc, char *argv[]) {
	struct timeval tv1, tv2;
	struct timezone tz;
	
   	double *x, *y, *z, norm;
   	int i, j;
   	
   	int n = atoi(argv[1]);
	
   	x = malloc(n*n*sizeof(double));
   	y = malloc(n*n*sizeof(double));  
	z = malloc(n*n*sizeof(double));
   	initMat(n, x);
   	initMat(n, y);
   	
   	gettimeofday(&tv1, &tz);
   	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, x, n, y, n, 0.0, z, n);
   	
   	norm = 0.;
   
	for(i=0;i<n;i++){
    	double row_sum = 0.;
        for(j=0;j<n;j++){
        	row_sum += z[i*n+j];
        }
        if(row_sum>norm){
        	norm = row_sum;
        }
    }
   	gettimeofday(&tv2, &tz);
   	double elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;
	printf("Time = %f\n",elapsed);
	
  	//printMat(n, z, "z");
 	//printf("\nRow Sum Norm = %f\n", norm);
 	
   	free(x);   	
   	free(y);   
   	free(z);       
   	return(0);
}
