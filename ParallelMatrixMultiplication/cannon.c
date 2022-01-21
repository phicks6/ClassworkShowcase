/*
	Philip Hicks
	Implements Cannon's algorithm for matrix multiplication
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

//Prints a matrix
void printMatrix(int * M, int dim){
	int i;
	int j;
	for(j = 0; j < dim; j++){
		for(i = 0; i < dim; i++){
			printf("%5d ",M[i + j*dim]);
		}
		printf("\n");
	}
}

//Prints a matrix to a file 
void fprintMatrix(FILE * fout, int * M, int dim){
	int i;
	int j;
	for(j = 0; j < dim; j++){
		for(i = 0; i < dim; i++){
			fprintf(fout,"%3d ",M[i + j*dim]);
		}
		fprintf(fout, "\n");
	}
}

//Preforms the local matrix multiplication
void matrixMultiply(int *a, int *b, int *c, int dim ) {
	int i;
	int j;
	int k;
	for(j = 0; j < dim; j++){
		for(i = 0; i < dim; i++){
			for(k = 0; k < dim; k++){
				c[i+j*dim] += a[k+(j)*dim] * b[i+(k)*dim];
			}
		}
	}
}


int main ( int argc, char *argv[] ){
	int rank;
	int numP;
	
	double startComm, endComm, startComp, endComp, startTotal, endTotal;
	double commTime = 0.0;
	double compTime = 0.0;
	int N;
	N = atoi(argv[1]);
	
	//Initialize mpi and processors
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numP);
	
	
	int numPinRow = sqrt(numP);
	int blockDim = N/numPinRow;
	int blockSize = N*N/numP;
	

	//Map processors into 2D grid
	int dims[2] = {0, 0};
    MPI_Dims_create(numP, 2, dims);
	int periods[2] = {1, 1};
	MPI_Comm cartComm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cartComm);
	

	//Get coords of processor
	int coords[2];
	MPI_Cart_coords(cartComm, rank, 2, coords);

	//Creates room for the matrixes
	int * subMatrixA = malloc(sizeof(int)*blockSize);
	int * subMatrixB = malloc(sizeof(int)*blockSize);
	int * resultMatrix = malloc(sizeof(int)*blockSize);
	
	int i;
	for(i = 0; i < blockSize; i++){
		resultMatrix[i] = 0;
	}

	int* sendCounts = (int*)malloc(sizeof(int) * numP);
	int* displacements = (int*)malloc(sizeof(int) * numP);
	int * fullMatrix;
	int * fullResultMatrix;

	//Creates a MPI_Datatype so that the MPI_Scatterv can be used to properly distribute the matrix from the root processor
	MPI_Datatype type, resizedtype;
	int sizes[2]    = {N,N}; 
	int subsizes[2] = {N/numPinRow,N/numPinRow}; 
	int starts[2]   = {0,0};  

	MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &type);
	MPI_Type_create_resized(type, 0, blockDim*sizeof(int), &resizedtype);
	MPI_Type_commit(&resizedtype);
    
	//Rank 0 creates the NxN matrix
	if(rank == 0){
		fullMatrix = malloc(sizeof(int)*N*N);
		fullResultMatrix = malloc(sizeof(int)*N*N);
		int i;
		int j;
		for(i = 0; i < N*N; i++){
			fullMatrix[i] = i%25; //can be used to show it works on non just 1's matrix
			//fullMatrix[i] = 1; //The result matrix will just be N in every slot
		}
		//printMatrix(fullMatrix,N);


		for (i = 0; i < numP; i++) {
			sendCounts[i] = 1;
		}

		//Caculate the displacements so we can use Scatterv to effciently send the submatrixes to the right processors
		int disp = 0;
		for(j = 0; j < numPinRow; j++){
			
			for(i = 0; i < numPinRow; i++){
				displacements[i+ (j*numPinRow)] = disp;
				disp += 1;
			}
			disp += (blockDim - 1)* numPinRow;
		}
	}


	//Start timing
	startTotal = MPI_Wtime();
	startComm = MPI_Wtime();


	//Scatter matrix twice to all processors
	MPI_Scatterv(fullMatrix, sendCounts, displacements, resizedtype, subMatrixA, blockDim*blockDim, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(fullMatrix, sendCounts, displacements, resizedtype, subMatrixB, blockDim*blockDim, MPI_INT, 0, MPI_COMM_WORLD);
	

	//Do initializing shift because we can't do caculations with given subMatrix 
	int left, right, up, down;
	MPI_Cart_shift(cartComm, 1, coords[0], &left, &right);
	MPI_Sendrecv_replace(subMatrixA, blockSize, MPI_INT, left, 1, right, 1, cartComm, MPI_STATUS_IGNORE);
	MPI_Cart_shift(cartComm, 0, coords[1], &up, &down);
	MPI_Sendrecv_replace(subMatrixB, blockSize, MPI_INT, up, 1, down, 1, cartComm, MPI_STATUS_IGNORE);

	endComm = MPI_Wtime();
	commTime+= endComm-startComm;

	startComp = MPI_Wtime();
	//Do first computation pass
	if(rank == -1){ //Set rank == 0 to get a set by set printout from one processor
			printf("\n\nSub matrix A:\n");
			printMatrix(subMatrixA,blockDim);
			printf("\nSub matrix B:\n");
			printMatrix(subMatrixB,blockDim);
			printf("\n\n");
			matrixMultiplyVerbose(subMatrixA,subMatrixB,resultMatrix,blockDim);
	}else{
		matrixMultiply(subMatrixA,subMatrixB,resultMatrix,blockDim);
	}

	endComp = MPI_Wtime();
	compTime+= endComp-startComp;


	//Do sqrt(P)-1 shifts then computations
	for(i = 0; i < numPinRow-1; i++){
		//Do shift
		startComm = MPI_Wtime();
		MPI_Cart_shift(cartComm, 1, 1, &left, &right);
		MPI_Cart_shift(cartComm, 0, 1, &up, &down);
		MPI_Sendrecv_replace(subMatrixA, blockSize, MPI_INT, left, 1, right, 1, cartComm, MPI_STATUS_IGNORE);
		MPI_Sendrecv_replace(subMatrixB, blockSize, MPI_INT, up, 1, down, 1, cartComm, MPI_STATUS_IGNORE);
		endComm = MPI_Wtime();
		commTime+= endComm-startComm;

		startComp = MPI_Wtime();

		//Do local matrix multiplication
		if(rank == -1){
			printf("\n\nSub matrix A:\n");
			printMatrix(subMatrixA,blockDim);
			printf("\nSub matrix B:\n");
			printMatrix(subMatrixB,blockDim);
			printf("\n\n");
			matrixMultiplyVerbose(subMatrixA,subMatrixB,resultMatrix,blockDim);
		}else{
			matrixMultiply(subMatrixA,subMatrixB,resultMatrix,blockDim);
		}
		endComp = MPI_Wtime();
		compTime+= endComp-startComp;
	}
	
	startComm = MPI_Wtime();

	//Gather results so the whole matrix is on every processor
	MPI_Gatherv(resultMatrix, blockSize, MPI_INT,
		fullResultMatrix, sendCounts, displacements, resizedtype,
		0, MPI_COMM_WORLD);
	endComm = MPI_Wtime();
	commTime+= endComm-startComm;

	endTotal = MPI_Wtime();

	//Have root processor print result and statistics
	if(rank == 0){
		printf("Result is:\n");
		printMatrix(fullResultMatrix, N);
		
		printf("%3d processors COMMed for : %.10f\n",numP,commTime);
		printf("               COMPed for : %.10f\n",compTime);
		printf("           for a total of : %.10f\n",endTotal-startTotal);


		FILE * fout = fopen("ResultMatrix.txt","w");
		fprintMatrix(fout,fullResultMatrix, N);

		FILE * tout = fopen("cannonTiming.txt","a");
		fprintf(tout,"%d processors\n %.10f, %.10f, %.10f \n",numP, commTime, compTime, endTotal-startTotal);
		fclose(fout);
		fclose(tout);
	}

	
	MPI_Finalize();
    return 1;
}

