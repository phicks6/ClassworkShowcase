/*
	Philip Hicks
	Implements DNS algorithm for matrix multiplication
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
	
	
	int numPinRow = (int)cbrt(numP);
	int numPinGrid = numPinRow*numPinRow;
	int blockDim = (N/numPinRow);
	int blockSize = blockDim*blockDim;
	
	//Map processors into 3D grid
	int dims[3] = {0, 0, 0};
    MPI_Dims_create(numP, 3, dims);
	int periods[3] = {1, 1 , 1};
	MPI_Comm cartComm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &cartComm);
	
	
	//Make a bunch of comms and groups for sub rings and sub grid
	MPI_Comm kgridComm, columnComm, iComm, jComm;
	int gridDims[3] = {0, 1, 1};
    MPI_Cart_sub(cartComm, gridDims, &kgridComm);

	int columnDims[3] = {1, 0, 0};
	MPI_Cart_sub(cartComm, columnDims, &columnComm);
	int iDims[3] = {0, 0, 1};
	MPI_Cart_sub(cartComm, iDims, &iComm);
	int jDims[3] = {0, 1, 0};
	MPI_Cart_sub(cartComm, jDims, &jComm);

	MPI_Group worldGroup,columnGroup,iGroup, jGroup, cartGroup;
	MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);
	MPI_Comm_group(columnComm, &columnGroup);
	MPI_Comm_group(iComm, &iGroup);
	MPI_Comm_group(jComm, &jGroup);
	MPI_Comm_group(cartComm, &cartGroup);
	
	
	//Get coords of processor
	int coord3d[3];
	int coord2d[2];
	MPI_Cart_coords(cartComm, rank, 3, coord3d);
	
	//Creates room for the matrixes
	int * subMatrixA = malloc(sizeof(int)*blockSize);
	int * subMatrixB = malloc(sizeof(int)*blockSize);
	int * resultMatrix = malloc(sizeof(int)*blockSize);
	
	int i;
	for(i = 0; i < blockSize; i++){
		resultMatrix[i] = 0;
	}

	int* sendCounts = (int*)malloc(sizeof(int) * numPinGrid);
	int* displacements = (int*)malloc(sizeof(int) * numPinGrid);
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

		
		for (i = 0; i < numPinGrid; i++) {
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


	//Scatter matrix twice along bottom grid of processors
	if(rank < numPinGrid){ 
		MPI_Scatterv(fullMatrix, sendCounts, displacements, resizedtype, subMatrixA, blockDim*blockDim, MPI_INT, 0, kgridComm);
		MPI_Scatterv(fullMatrix, sendCounts, displacements, resizedtype, subMatrixB, blockDim*blockDim, MPI_INT, 0, kgridComm);
	}
	MPI_Barrier(MPI_COMM_WORLD);


	//Send B up the matrix to correct processors
	int columnRank = 0;
	int senderRank;
	int receiverRank;
	int senderCoords[3];
	MPI_Group_translate_ranks(columnGroup, 1, &columnRank, worldGroup, &senderRank);
	MPI_Cart_coords(cartComm, senderRank, 3, senderCoords);
	MPI_Group_translate_ranks(columnGroup, 1, &senderCoords[1], worldGroup, &receiverRank);
	if(senderRank!=receiverRank){
		if(rank == senderRank){
			MPI_Send(subMatrixB, blockSize, MPI_INT, receiverRank, 0, MPI_COMM_WORLD);
		}
		if(rank == receiverRank){
			MPI_Recv(subMatrixB, blockSize, MPI_INT, senderRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}

	//Send A up the matrix to correct processors
	MPI_Group_translate_ranks(columnGroup, 1, &senderCoords[2], worldGroup, &receiverRank);
	if(senderRank!=receiverRank){
		if(rank == senderRank){
			MPI_Send(subMatrixA, blockSize, MPI_INT, receiverRank, 0, MPI_COMM_WORLD);
		}
		if(rank == receiverRank){
			MPI_Recv(subMatrixA, blockSize, MPI_INT, senderRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}

	

	
	//Spread B matrix along j
	int broadCaster = rank/numPinGrid;
	MPI_Bcast(subMatrixB, blockSize, MPI_INT, broadCaster, jComm);
	MPI_Barrier(MPI_COMM_WORLD);

	
	//Spread A matrix along i
	broadCaster = rank/numPinGrid;
	MPI_Bcast(subMatrixA, blockSize, MPI_INT, broadCaster, iComm);
	MPI_Barrier(MPI_COMM_WORLD);
	
	endComm = MPI_Wtime();
	commTime+= endComm-startComm;

	startComp = MPI_Wtime();

	//Do the single local matrix multiplication
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
	
	//The results are still distributed along the collums of our Cube of processors
	//So we need to add them together and put it on the bottom grid of processors
	int * addedResultMatrix;
	if(rank < numPinGrid){
		addedResultMatrix = malloc(sizeof(int)*blockSize);
	}
	MPI_Comm_rank(columnComm, &columnRank);

	startComp = MPI_Wtime();
	//Reduce results to bottom grid
	MPI_Reduce(resultMatrix,addedResultMatrix,blockSize,MPI_INT, MPI_SUM, 0, columnComm);
	endComp = MPI_Wtime();
	compTime+= endComp-startComp;
	MPI_Barrier(MPI_COMM_WORLD);
	

	startComm = MPI_Wtime();
	if(rank < numPinGrid){
		//Gather results so the whole matrix is on every processor of the bottom grid
		MPI_Gatherv(addedResultMatrix, blockSize, MPI_INT,
			fullResultMatrix, sendCounts, displacements, resizedtype,
			0, kgridComm);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	endComm = MPI_Wtime();
	commTime+= endComm-startComm;
	endTotal = MPI_Wtime();

	//Have root processor print result and statistics
	if(rank == 0){
		printf("Result is:\n");
		printMatrix(fullResultMatrix, N);
		
		printf("%3d processors communicated for : %.10f\n",numP,commTime);
		printf("               computated for : %.10f\n",compTime);
		printf("           for a total of : %.10f\n",endTotal-startTotal);

		
		FILE * fout = fopen("ResultMatrix.txt","w");
		fprintMatrix(fout,fullResultMatrix, N);
		
		FILE * tout = fopen("/nics/b/home/phicks6/PA3/dnsTiming.txt","a");
		fprintf(tout,"%d processors\n %.10f, %.10f, %.10f \n",numP, commTime, compTime, endTotal-startTotal);
		fclose(fout);
		fclose(tout);
	}
	
	MPI_Finalize();
    return 1;
}
