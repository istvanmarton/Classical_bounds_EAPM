#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <omp.h>
#include "functions.h"

void calc_jmin_jmax(int* index, unsigned long long int* jMin, unsigned long long int* jMax, unsigned long long int* steps, unsigned long long int* steps_remainder){
	*jMax = (*index + 1) * *steps - 1;
	*jMin = *index * *steps;
	if(*index < *steps_remainder) *jMax += *index + 1;
	else *jMax += *steps_remainder;
	if(*index <= *steps_remainder) *jMin += *index;
	else *jMin += *steps_remainder;
}

void Ln(std::vector<std::vector<int>>& indexes_local, int_type* d_mtx_as_vec, int* d_partitions, unsigned long long int steps, unsigned long long int steps_remainder, int_type *Ln_vector, int cube_size, int n){ // This function calculates the Ld norm
	int i, l, index;
	std::vector<int> i_Data;
	int_type temp[16][16], product, Ln = 0;
	unsigned long long int jMax, jMin, ix, jx, jy, idx;

	index = omp_get_thread_num(); // Index of threads.
	calc_jmin_jmax(&index, &jMin, &jMax, &steps, &steps_remainder); // This function calculates the minimal (jMin-th) and the maximal (jMax-th) word of the binary reflected Gray code for which the calculations must be performed by a given thread.
	for(ix = jMin; ix <= jMax; ix++){
		for(i = 0; i < n*n; i++){
			for(l = 0; l < cube_size; l++){
				temp[i][l] = 0;
			}
		}
		
		for(jy = 0; jy < cube_size; jy++){
			for(jx = 0; jx < cube_size; jx++){
				idx = n * d_partitions[ix * cube_size + jx] + d_partitions[ix * cube_size + jy];
				for(i = 0; i < cube_size; i++) {temp[idx][i] += d_mtx_as_vec[i * cube_size * cube_size + jy * cube_size + jx];
				}
			}
		}
	
		product = 0;
		for(i = 0; i < n*n ; i++){
			for(l = 0; l < cube_size; l++){
				product += abs(temp[i][l]); //printf("%d ", temp[i][l]);
			}
		}

		if(product > Ln) {
		        indexes_local = std::vector<std::vector<int>>();
			Ln = product; // If the current Ln sum is greater than the previous one, it modifies both the value
		}

        if(product == Ln) {
            i_Data.push_back(Ln);
            for(l=0; l < cube_size; l++){i_Data.push_back(d_partitions[ix * cube_size + l]);}
            indexes_local.push_back(i_Data);
            i_Data.clear();
        }
	}
	Ln_vector[index] = Ln; // Every thread writes the biggest found Ln sum to the d_Ln_vector.
//	fclose(fp);
}

void print_partition(int* partition, item_calc* second) {
  int i, *temp;
  temp = (int*) realloc(second->partitions, (second->counter) * (second->cube_size) * sizeof(int));
  if(!temp){
    printf("Memory allocation failed!");
    return;
  }
  second->partitions = temp;
  for(i = 0; i < second->cube_size; i++){
    (second->partitions)[(second->counter - 1) * (second->cube_size) + i] = partition[i];
  }
}

void backtrack(int index, int max_used, int* partition, item_calc* second) {
	
	if (index == second->cube_size) {
        if (max_used == second->n) {
            (second->counter)++;
            print_partition(partition, second);
        }
        return;
    }

    for (int i = 0; i < max_used; i++) {
        partition[index] = i;
        backtrack(index + 1, max_used, partition, second);
    }

    if (max_used < second->n) {
        partition[index] = max_used;
        backtrack(index + 1, max_used + 1, partition, second);
    }
}

void calc_Parameters(item* first, item_calc* second, int* num_of_threads){
	int *partition;
	
	double bigger = second->iRows_reduced > second->iCols_reduced ? second->iRows_reduced : second->iCols_reduced;
	second->cube_size = (int) cbrt(bigger);
	printf("size of tensor: %d\n", second->cube_size);
	partition = (int*) calloc(second->cube_size, sizeof(int));
    	second->partitions = NULL;
    	second->counter = 0;
    	partition[0] = 0;
printf("iRows_reduced1: %d, iCols_reduced1: %d, second->n: %d\n", second->iRows_reduced, second->iCols_reduced, second->n);
	backtrack(1, 1, partition, second);
	first->total_num_to_calc = second->counter;
	
	second->copyNum = *num_of_threads > first->total_num_to_calc ? first->total_num_to_calc : *num_of_threads; // copyNum is the actual number of threads. It cannot be more than the number of possible L values.
	second->steps = second->copyNum == 0 ? 0 : first->total_num_to_calc/ second->copyNum;
	second->steps_remainder = second->copyNum == 0 ? 0 : first->total_num_to_calc % second->copyNum;

	second->steps = second->copyNum == 0 ? 0 : first->total_num_to_calc/ second->copyNum;
	second->steps_remainder = second->copyNum == 0 ? 0 : first->total_num_to_calc % second->copyNum;

	printf("num_ofThread: %d\n", *num_of_threads);
	if(second->iCols_reduced > length) {printf("Matrix is too large. The length variable %d should be greater or equal than %d.\n", length, second->iCols_reduced); exit(-1);}
	free(partition);
}

void calc_Lnorm(item* first, item_calc* second, int* num_of_threads){
std::string fileName = std::string("strategy_")+ first->fileName;
        
	int i, iMax;
	int_type *Ln_vector;
	Ln_vector = (int_type*)calloc(second->copyNum, sizeof(int_type)); // The code allocates memory in the host for the possible L norms.
	second->strategy = (int*) calloc(second->iRows_reduced , sizeof(int));
#pragma omp parallel num_threads(second->copyNum)
{
        std::vector<std::vector<int>> indexes_local;
	Ln(indexes_local, second->mtx_as_vec, second->partitions, second->steps, second->steps_remainder, Ln_vector, second->cube_size, second->n);
#pragma omp critical
        {
                second->indexes.insert(std::end(second->indexes), std::begin(indexes_local), std::end(indexes_local)); //ionisation_px.insert(std::end(IonisationStatist->ionisation_px), std::begin(local_IonisationStatist.ionisation_px), std::end(local_IonisationStatist.ionisation_px));
        }
}
        
//	printf("second->copyNum: %llu\n", second->copyNum);
	second->Lnorm = Ln_vector[0];
	iMax = 0; // Determining the maximal element of Ln_vector, which is the L norm, and the index of the corresponding strategy vector as well.
	for(i = 1; i < second->copyNum; i++){ if(second->Lnorm < Ln_vector[i]) {second->Lnorm = Ln_vector[i]; iMax = i;}}
	
        std::ofstream myfile;
        myfile.open(fileName);
        for(std::vector<std::vector<int>>::iterator it = second->indexes.begin(); it != second->indexes.end(); ++it){
            if( *(it->begin()) == second->Lnorm ){
                for(std::vector<int>::iterator col = it->begin()+1; col != it->end(); ++col){
			myfile << (*col) << ' ';
                }
        	myfile << '\n';
            }
        }
        myfile.close();
	
	free(Ln_vector); // Deallocates the vectors in the host memory.
}

void arguments_OMP(item* first, int* num_of_threads, int* argc, char** argv){
	FILE *fp;
	int sd, t;
	char msg[] = "Use the following command: ./L_OpenMP number_of_threads filename_of_matrix order_of_the_L_norm";
	if(*argc < 4){
		printf("Incorrect number of input arguments. %s\n", msg);
		exit(-1);
	}
	
	sd = sscanf(argv[1], "%d", &t);
	if((sd == 0) || (t < 1)){
		printf("Please make sure that the number of threads is a positive integer. %s\n", msg);
		exit(-1);
	}
	*num_of_threads = t;
	
	sprintf(first->fileName,"%s", argv[2]);
	fp = fopen(first->fileName, "r");
	if(fp == NULL){
		printf("Please make sure that the file containig the matrix exists within this directory. %s\n", msg);
	}
	fclose(fp);
	
	sd = sscanf(argv[3], "%d", &t);
	if((sd == 0) || (t < 2)){
		printf("Please make sure that the order of the L norm is a bigger than 1. %s\n", msg);
		exit(-1);
	}
	first->n_original = t;
	
	if(*argc < 5){first->stat = 'n';}
	else {first->stat = argv[4][0];}
	
	if(t > RANK_OF_NORM) {printf("The order of the L norm is too large. Please increase the RANK_OF_NORM variable in the code to %d and compile and run it again.\n", t); exit(-1);}
}

void load_parameters(item* first, item_calc* second, int* argc, char** argv, int* num_of_threads){
	arguments_OMP(first, num_of_threads, argc, argv);
	matrix_read(first);
	second->iRows_reduced = first->iRows;
	second->iCols_reduced = first->iCols;
	second->n = first->n_original;
	convert_mtx_to_vec_noTranspose(first, second);
	calc_Parameters(first, second, num_of_threads); // Calculates the necessary parameters for the calculation.
}

int main(int argc, char *argv[]){
	int num_of_threads;
	item first;
	item_calc second;

	load_parameters(&first, &second, &argc, argv, &num_of_threads);
	
	if(second.iRows_reduced == 0) {printf("This is a zero matrix.\n"); second.Lnorm = 0;}
	else {calc_Lnorm(&first, &second, &num_of_threads);} // The function 'calc_Lnorm' calculates the L norm of order n of the input matrix.
	print_results(&first, &second);
	
	free_first(&first);
	if(second.iRows_reduced > 0) free_second(&second);
	return 0;
}
