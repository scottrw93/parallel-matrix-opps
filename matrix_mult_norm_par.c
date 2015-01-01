#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <sys/time.h>

#define MAXTHRDS 124

typedef struct {
	double *x;
   	double *y;
   	double *z;
   	int n;
   	int sub_n;
}mat_mult_t;

typedef struct {
   	double *z;
   	int n;
   	int sub_n;
   	double *global_norm;
   	pthread_mutex_t *mutex;
}mat_norm_t;

void *matrix_mult(void *arg){
   	mat_mult_t *thread_mat_mult_data = arg;
   	int n = thread_mat_mult_data->n;
   	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, thread_mat_mult_data->sub_n, n, n,
   	 1.0, thread_mat_mult_data->x, n, thread_mat_mult_data->y, n, 0.0, thread_mat_mult_data->z, n);
   	pthread_exit(NULL);
}

void *matrix_norm(void *arg){
   	mat_norm_t *thread_mat_norm_data = arg;
   	int n = thread_mat_norm_data->n;
   	int i, j;
   	double norm = 0.;
   	
   	for(i=0;i<thread_mat_norm_data->sub_n;i++){
   		double row_sum = 0.;
   		for(j=0;j<n;j++){
   			row_sum += *(thread_mat_norm_data->z+i*n+j);
   		}
   		if(row_sum>norm){
   			norm = row_sum;
   		}
   	}

	pthread_mutex_lock(thread_mat_norm_data->mutex);
   	if (norm > *(thread_mat_norm_data->global_norm)){
    		*(thread_mat_norm_data->global_norm)=norm;
   	}
   	pthread_mutex_unlock(thread_mat_norm_data->mutex);	
   	pthread_exit(NULL);
}

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

	pthread_t *working_thread;
   	mat_mult_t *thread_mat_mult_data;
   	mat_norm_t *thread_mat_norm_data;
   	pthread_mutex_t *mutex;
   	double *x, *y, *z, norm;
   	int i, rows_per_thread;
   	void *status;
   	
   	int n = atoi(argv[1]); 	
	int num_of_thrds = atoi(argv[2]);// Works when this is 2, not when 4 
	
	if(n<=num_of_thrds && num_of_thrds < MAXTHRDS){
		printf("Matrix dimension must be greater than num of thrds\nand num of thrds less than 124.\n");
		return (-1);
	}
	
   	x = malloc(n*n*sizeof(double));
   	y = malloc(n*n*sizeof(double));  
	z = malloc(n*n*sizeof(double));
   	initMat(n, x);
   	initMat(n, y);
   	
   	working_thread = malloc(num_of_thrds * sizeof(pthread_t));
   	thread_mat_mult_data = malloc(num_of_thrds * sizeof(mat_mult_t));
	rows_per_thread = n/num_of_thrds;
	gettimeofday(&tv1, &tz);
   	
	for(i=0;i<num_of_thrds;i++){
      		thread_mat_mult_data[i].x = x + i * rows_per_thread * n;
      		thread_mat_mult_data[i].y = y;
      		thread_mat_mult_data[i].z = z + i * rows_per_thread * n;
      		thread_mat_mult_data[i].n = n;
      		thread_mat_mult_data[i].sub_n = (i == num_of_thrds-1) ? n-(num_of_thrds-1)*rows_per_thread : rows_per_thread;
      		pthread_create(&working_thread[i], NULL, matrix_mult, (void *)&thread_mat_mult_data[i]);
   	}
   	
   	for(i=0;i<num_of_thrds;i++){
        pthread_join(working_thread[i], NULL);
    }

    	//printMat(n, z , "z");
    
    working_thread = malloc(num_of_thrds * sizeof(pthread_t));
    thread_mat_norm_data = malloc(num_of_thrds * sizeof(mat_norm_t));
    mutex = malloc(sizeof(pthread_mutex_t));
	pthread_mutex_init(mutex, NULL);
    
	for(i=0;i<num_of_thrds;i++){
    	  	thread_mat_norm_data[i].z = z + i * rows_per_thread * n;
      		thread_mat_norm_data[i].n = n;
      		thread_mat_norm_data[i].global_norm = &norm;
      		thread_mat_norm_data[i].sub_n = (i == num_of_thrds-1) ? n-(num_of_thrds-1)*rows_per_thread : rows_per_thread;
      		thread_mat_norm_data[i].mutex = mutex;
      		pthread_create(&working_thread[i], NULL, matrix_norm, (void *)&thread_mat_norm_data[i]);
   	}

   	for(i=0;i<num_of_thrds;i++){
        pthread_join(working_thread[i], &status);
	}
    
	gettimeofday(&tv2, &tz);
 	//printf("\nRow Sum Norm = %f\n", norm);
 	double elapsed = (double) (tv2.tv_sec-tv1.tv_sec) + (double) (tv2.tv_usec-tv1.tv_usec) * 1.e-6;
	printf("Time = %f\n",elapsed);
 	
   	free(x);   	
   	free(y);   
   	free(z);       
   	free(working_thread);
   	free(thread_mat_mult_data);
   	free(thread_mat_norm_data);
   	pthread_mutex_destroy(mutex);
   	free(mutex);

   	return(0);
}
