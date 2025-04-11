#include<stdio.h>
#include<vector>
#include<fstream>
#define length 4096
#define RANK_OF_NORM 20
#define NUM_OF_BITS 8 * sizeof(unsigned long long int)
typedef int int_type;

struct data_k{
	int n;
	int iRows_reduced;
	int iCols_reduced;
	int_type Lnorm;
	unsigned long long int steps;
	unsigned long long int steps_remainder;
	unsigned long long int copyNum;
	int_type *mtx_as_vec;
	int *strategy;
	unsigned long long int counter;
	int cube_size;
	int *partitions;
	std::vector<std::vector<int>> indexes;
};

typedef struct data_k item_calc;

struct data_s{
	int_type **matrix;
	int n_original;
	int iRows;
	int iCols;
	char fileName[1024];
	char stat;
	int original_length;
	unsigned long long int total_num_to_calc;
};

typedef struct data_s item;

void matrix_read(item* first){
	int i = 0, j = 0, k = 0;
	int_type value, *row;
	first->matrix = NULL;
	row = NULL;

	char g, cNum[256];

	FILE *fp;
	fp = fopen(first->fileName,"r");
	
	do{
		g = fgetc(fp);	
		if((((g - '0') < 10) && ((g - '0') >= 0)) || (g == 'e') || ( g == 'E') || (g == '.') || (g == '+') || (g == '-')) {cNum[i] = g; i++;}
		else {
			cNum[i] = '\0'; 
			if(cNum[0] != '\0') {sscanf(cNum, "%d", &value); j++; i = 0;  row = (int_type*) realloc(row, j * sizeof(int_type)); row[j-1] = value;}
			if( ((g == '\n') || (g == EOF)) && (j > 0)){first->iCols = j; j = 0; k++; first->matrix = (int_type**) realloc(first->matrix, k * sizeof(int_type*)); first->matrix[k-1] = row; row = NULL;}
		}
	}while(!feof(fp));
	first->iRows = k;
	printf("rows: %d, cols: %d\n",first->iRows, first->iCols); 
	fclose(fp);
}

void convert_mtx_to_vec_noTranspose(item* first, item_calc* second){
	int i, j;
	second->mtx_as_vec = (int_type*) calloc(second->iRows_reduced * second->iCols_reduced, sizeof(int_type));
	for(i = 0; i < second->iRows_reduced; i++){
		for(j = 0; j < second->iCols_reduced; j++){
			second->mtx_as_vec[i * second->iCols_reduced + j] = first->matrix[i][j];
		}
	}
	first->original_length = first->iRows;
}

void mtx_free(item* first){
	int i;
	for(i = 0; i < first->iRows; i++){
		free(first->matrix[i]);
	}
	free(first->matrix);
}

void free_first(item* first){
	mtx_free(first);
}

void free_second(item_calc* second){
	free(second->strategy);
	free(second->mtx_as_vec);
}

int* calc_Pattern(int* n){
	int i, *iPattern;
	iPattern = (int*) calloc(2 * *n, sizeof(int));
	for(i = 0; i < *n; i++){
		iPattern[i] = i;
		iPattern[2 * *n-i-1]=i;
	}
	return iPattern;
}

void print_results(item* first, item_calc* second){
	int i, j, val;
	char fileOutput[1024]; // The variable 'fileOutput' is the name of the file to which the strategy vector found to be optimal is written
	FILE *fp;
		 
	printf("L%d is: %d\n", first->n_original, second->Lnorm); // Write out the value of the L norm to the screen.
}
