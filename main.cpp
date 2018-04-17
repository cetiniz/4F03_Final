#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <mpi.h>
#include <math.h>

#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#define epsilon 0.000000000000000222
#define g pow(6.673*10, -11)

double **contigArrayGenerator(int row, int col){
	double **contigarray = (double **)malloc(row*sizeof(double));
	double *pointer = (double *)malloc(row*col*sizeof(double));
	for(int i = 0; i < row; i++){
		contigarray[i] = &(pointer[col*i]);
	}
	return contigarray;
}

double computeForce(double x1, double y1, int w1, double x2, double y2, int w2, int axisToCompute){
	double force = 0;
	int dx = x1-x2; // compute position change in x axis
	int dy = y1-y2; // compute position change in y axis
	int distance = sqrt(dx*dx+dy*dy); // compute distance between particles
	int distance3 = distance*distance*distance; // compute distance cubed between particles

	/* COMPUTE FORCE */
	if(axisToCompute == 0){
		force = -g*w1*w1/distance3*dx;
	} else if (axisToCompute == 1){
		force = -g*w1*w1/distance3*dy;
	}
	return force;
}

int main(int argc, char* argv[]){
	
	if( argc != 10){
		printf("Usage: %s numParticlesLight numParticlesMedium numParticlesHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrex\n", argv[0]);
		return 1;
	}

	MPI_Init(&argc,&argv);

	int p, my_rank;
	int i,j,m,offset;

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	//variables
	int numParticlesLight = std::stoi(argv[1]);
	int numParticlesMedium = std::stoi(argv[2]);
	int numParticlesHeavy = std::stoi(argv[3]);
	int numParticlesTotal = numParticlesLight + numParticlesMedium + numParticlesHeavy; //total number of particles is sum of light, medium, heavy particle numbers
	int frameTotal = 10;
 	int * w = (int *) malloc(sizeof(int) * numParticlesTotal); //array to store weight of particles
 	double * s_x = (double *) malloc(sizeof(double) * numParticlesTotal); //array to store positions of particles in x dimension
 	double * s_y = (double *) malloc(sizeof(double) * numParticlesTotal); //array to store positions of particles in y dimension
 	double * v_x = (double *) malloc(sizeof(double) * numParticlesTotal); //array to store velocities of particles in x dimesion
 	double * v_y = (double *) malloc(sizeof(double) * numParticlesTotal); //array to store velocities of particles in y dimesion
 	double **f_x = contigArrayGenerator(numParticlesTotal,numParticlesTotal); //matrix to store forces of particles in x dimension
 	double **f_y = contigArrayGenerator(numParticlesTotal,numParticlesTotal); //matrix to store forces of particles in y dimension

 	int h = std::stoi(argv[6]);

 	int imageWidth = std::stoi(argv[7]);
 	int imageHeight = std::stoi(argv[8]);

 	int numSteps = 0;
 	int subSteps = 0;
 	double timeSubStep;

 	int width, height;
 	int particlesToReceive;
 	int particlesPerProcessor = numParticlesTotal/p;
 	int particlesRemaining = numParticlesTotal%p;

 	int * pointerForOriginalArray;
 	int * particlesToCompute_s_x;
 	int * particlesToCompute_s_y;
 	int * particleWeights;
 	int * particlesToCompute_v_x;
 	int * particlesToCompute_v_y;


 	int * pointerForLocalArray;
 	int * pointerForTempArray;
 	int * tempWeights;
 	int * localWeights;
 	int * localArray_s_x;
 	int * localArray_s_y;
 	int * localArray_f_x;
 	int * localArray_f_y;
 	int * tempArray_s_x;
 	int * tempArray_s_y;
 	int * tempArray_f_x;
 	int * tempArray_f_y;

 	int * weights;
 	int * forces_x;
 	int * forces_y;

 	unsigned char* image;

 	MPI_Status status;
 	int source = status.MPI_SOURCE;

 /***************************** MASTER TASK ***************************/
 	if(my_rank == 0){

 	/******** Allocate particle weight, position, and velocity to array ********/
 		for(i = 0; i < numParticlesTotal; i++){
 			if(i < numParticlesLight){
 				w[i] = 1;
 				s_x[i] = drand48()*imageWidth;
 				s_y[i] = drand48()*imageHeight;
 				v_x[i] = drand48();
 				v_y[i] = drand48();
 			} else if(i >= numParticlesLight && i < (numParticlesLight+numParticlesMedium)){
 				w[i] = 2;
 				s_x[i] = drand48()*imageWidth;
 				s_y[i] = drand48()*imageHeight;
 				v_x[i] = drand48();
 				v_y[i] = drand48();
 			} else{
 				w[i] = 3;
 				s_x[i] = drand48()*imageWidth;
 				s_y[i] = drand48()*imageHeight;
 				v_x[i] = drand48();
 				v_y[i] = drand48();
 			}
 		}

 		for (int frameNum = 0; frameNum < frameTotal; frameNum++) {
 			printf("My thread number is %d and my loop (frame) is %d\n", my_rank,frameNum);
 			
 			for (int dest = 1; dest < p; dest++){
 				printf("My thread number is %d and my loop (masterSetup) is %d\n", my_rank, dest);

 				/******* STEP 1: ALLOCATE NUMBER OF PARTICLES TO EACH PROCESSOR *******/
 				particlesToReceive = (dest < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;
 				/******* STEP 2: CREATE ARRAYS TO STORE PARTICLE VALUES & LOCATION IN ORIGINAL ARRAY (Particle number) *******/
 				particlesToCompute_s_x = (int *) malloc(sizeof(int) * particlesToReceive); 
 				particlesToCompute_s_y = (int *) malloc(sizeof(int) * particlesToReceive); 
 				particleWeights = (int *) malloc(sizeof(int) * particlesToReceive); 
 				pointerForOriginalArray = (int *) malloc(sizeof(int) * particlesToReceive); //contains pointers that store location of original matrix location (when the Master gathers everything back at the end)

 		/******* STEP 3: DISTRIBUTE PARTICLES FROM ORIGINAL ARRAYS TO NEW ARRAYS & MARK LOCATION IN ORIGINAL ARRAYS *******/
 				m=0;
 				offset = dest;
 				for(i = offset; i < numParticlesTotal; i+=p){
 					particlesToCompute_s_x[m] = s_x[i];
 					particlesToCompute_s_y[m] = s_y[i];
 					*particleWeights = w[i];
 					pointerForOriginalArray[m] = i;
 					m++;
 				}
				/******* SEND ARRAYS TO SLAVE PROCESSORS *******/
 				MPI_Send(&(particleWeights[0]), particlesToReceive, MPI_INT, dest, 0, MPI_COMM_WORLD);
 				MPI_Send(&(particlesToCompute_s_x[0]), particlesToReceive, MPI_INT, dest, 0, MPI_COMM_WORLD);
 				MPI_Send(&(particlesToCompute_s_y[0]), particlesToReceive, MPI_INT, dest, 0, MPI_COMM_WORLD);
 				MPI_Send(&(pointerForOriginalArray[0]), particlesToReceive, MPI_INT, dest, 0, MPI_COMM_WORLD);

 			}
			// /************** RING LOOP WILL GO HERE***************/
 			for(int ringNumber = 0; ringNumber < p - 1; ringNumber++){
 				printf("My thread number is %d and my loop (ringNumMaster) is %d\n", my_rank, ringNumber);

				//Send to dest AND receive from source
 				MPI_Sendrecv(&(particleWeights[0]), particlesToReceive, MPI_INT, (my_rank-1+p)%p, 0, &(tempWeights[0]), particlesToReceive, MPI_INT, (my_rank+1)%p, 0, MPI_COMM_WORLD, &status);
 				MPI_Sendrecv(&(localArray_f_x[0]),  particlesToReceive, MPI_INT, (my_rank-1+p)%p, 0, &(tempArray_s_x[0]), particlesToReceive, MPI_INT, (my_rank+1)%p, 0, MPI_COMM_WORLD, &status);
 				MPI_Sendrecv(&(localArray_f_y[0]),  particlesToReceive, MPI_INT, (my_rank-1+p)%p, 0, &(tempArray_s_y[0]), particlesToReceive, MPI_INT, (my_rank+1)%p, 0, MPI_COMM_WORLD, &status);
 				MPI_Sendrecv(&(pointerForOriginalArray[0]),  particlesToReceive, MPI_INT, (my_rank-1+p)%p, 0, &(pointerForTempArray[0]), particlesToReceive, MPI_INT, (my_rank+1)%p, 0, MPI_COMM_WORLD, &status);

				//Calculate forces


/**************************** MASTER RECIEVES TASKS FROM SLAVES **********************/
for(i = 1; i < p; i++){
	source = i;
	weights = (int *) malloc(particlesToReceive * sizeof(int));
	MPI_Recv(&(weights[0]), particlesToReceive, MPI_INT, source, 0, MPI_COMM_WORLD, &status);
	forces_x = (int *) malloc(particlesToReceive * sizeof(int));
	MPI_Recv(&(forces_x[0]), particlesToReceive, MPI_INT, source, 0, MPI_COMM_WORLD, &status);
	forces_y = (int *) malloc(particlesToReceive * sizeof(int));
	MPI_Recv(&(forces_y[0]), particlesToReceive, MPI_INT, source, 0, MPI_COMM_WORLD, &status);
	pointerForOriginalArray = (int *) malloc(particlesToReceive * sizeof(int)); //allocate memory for pointer array
	MPI_Recv(&(pointerForOriginalArray[0]), particlesToReceive, MPI_INT, source, 0, MPI_COMM_WORLD, &status);

	for(j = 0; j < particlesToReceive; j++){
		v_x[pointerForOriginalArray[j]] += h*forces_x[j]/weights[j];
		v_y[pointerForOriginalArray[j]] += h*forces_y[j]/weights[j];
	}

	for(j = 0; j < particlesToReceive; j++){
		s_x[pointerForOriginalArray[j]] += h*v_x[pointerForOriginalArray[j]];
		s_y[pointerForOriginalArray[j]] += h*v_y[pointerForOriginalArray[j]];
	}
}

saveBMP(argv[9], image, width, height);
}

 	}
 /*************************** SLAVE TASKS **********************************/
 	if(my_rank > 0){
 		localWeights = (int *) malloc(sizeof(int) * particlesToReceive); 
 		localArray_s_x = (int *) malloc(sizeof(int) * particlesToReceive); 
 		localArray_f_x = (int *) malloc(sizeof(int) * particlesToReceive); 
 		localArray_s_y = (int *) malloc(sizeof(int) * particlesToReceive); 
 		localArray_f_y = (int *) malloc(sizeof(int) * particlesToReceive); 
 		pointerForLocalArray = (int *) malloc(sizeof(int) * particlesToReceive); 
 		
 		tempArray_s_x = (int *) malloc(sizeof(int) * particlesToReceive); 
 		tempArray_f_x = (int *) malloc(sizeof(int) * particlesToReceive); 
 		tempArray_s_y = (int *) malloc(sizeof(int) * particlesToReceive); 
 		tempArray_f_y = (int *) malloc(sizeof(int) * particlesToReceive); 

 	/******* Recieve particles from MASTER *******/
 		for (int numFrames = 0; numFrames < frameTotal; numFrames++){
 			printf("My thread number is %d and my loop (slaveFrame) is %d\n", my_rank,numFrames);
 			MPI_Recv(&(particleWeights[0]), particlesToReceive, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
 			MPI_Recv(&(localArray_s_x[0]), particlesToReceive, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
 			MPI_Recv(&(localArray_s_y[0]), particlesToReceive, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
 			MPI_Recv(&(pointerForLocalArray[0]), particlesToReceive, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

 			for(i = 0; i < particlesToReceive; i++){
 				tempArray_s_x[i] = localArray_s_x[i];
 				tempArray_s_y[i] = localArray_s_y[i];
 			}


			// RING LOOP GOES HERE
 			for(int ringNumber = 0; ringNumber < p - 1; ringNumber++){
 				printf("My thread number is %d and my loop (slaveRingNumber) is %d\n", my_rank,ringNumber);
				/******* Send & Recieve particles from another SLAVE *******/
 				MPI_Sendrecv(&(particleWeights[0]), particlesToReceive, MPI_INT, (my_rank-1+p)%p, 0, &(tempWeights[0]), particlesToReceive, MPI_INT, (my_rank+1)%p, 0, MPI_COMM_WORLD, &status);
 				MPI_Sendrecv(&(localArray_s_x[0]),  particlesToReceive, MPI_INT, (my_rank-1+p)%p, 0, &(tempArray_s_x[0]), particlesToReceive, MPI_INT, (my_rank+1)%p, 0, MPI_COMM_WORLD, &status);
 				MPI_Sendrecv(&(localArray_s_y[0]),  particlesToReceive, MPI_INT, (my_rank-1+p)%p, 0, &(tempArray_s_y[0]), particlesToReceive, MPI_INT, (my_rank+1)%p, 0, MPI_COMM_WORLD, &status);
 				MPI_Sendrecv(&(pointerForOriginalArray[0]),  particlesToReceive, MPI_INT, (my_rank-1+p)%p, 0, &(pointerForTempArray[0]), particlesToReceive, MPI_INT, (my_rank+1)%p, 0, MPI_COMM_WORLD, &status);

				//Calculate forces
 				j=0;
				for(i = 0; i < particlesToReceive; i++){
				 	while(pointerForLocalArray[i] >= pointerForTempArray[j] && j < particlesToReceive){ //find index where particle number in tempArray is greater than localArray
				 		j++;
				 	}
				 	if(pointerForLocalArray[i] < pointerForTempArray[j]){
				 		localArray_f_x[i] += computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[i], tempArray_s_y[i], tempWeights[i], 0);
				 		localArray_f_y[i] += computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[i], tempArray_s_y[i], tempWeights[i], 1);
				 		tempArray_f_x[j] -= computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[i], tempArray_s_y[i], tempWeights[i], 0);
				 		tempArray_f_y[j] -= computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[i], tempArray_s_y[i], tempWeights[i], 1);
				 	}
				 }




		}
		// FINAL SEND GOES HERE
		MPI_Send(&(particleWeights[0]), particlesToReceive, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&(localArray_f_x[0]), particlesToReceive, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&(localArray_f_y[0]), particlesToReceive, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&(pointerForOriginalArray[0]), particlesToReceive, MPI_INT, 0, 0, MPI_COMM_WORLD);

	}
}

free(image);

MPI_Finalize();
return 0;
}
