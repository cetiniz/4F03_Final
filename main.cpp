#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <mpi.h>
#include <math.h>

#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#define _XOPEN_SOURCE
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

double computeForce(int x1, int y1, int w1, int x2, int y2, int w2, int axisToCompute){
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
void printArray(int * arr, int length){
	printf("PRINTING OUT ARRAY\n");
	for (int i = 0; i < length; i++) {
		printf("%d ",arr[i]);
	}
	printf("\n");
}

int main(int argc, char* argv[]){
	
	if( argc != 10){
		printf("Usage: %s numParticlesLight numParticlesMedium numParticlesHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrefix\n", argv[0]);
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
	int numSteps = std::stoi(argv[4]);
	int numSubSteps = std::stoi(argv[5]);
	int timeSubSteps = std::stoi(argv[6]);
	int numParticlesTotal = numParticlesLight + numParticlesMedium + numParticlesHeavy; //total number of particles is sum of light, medium, heavy particle numbers
 	int * w = (int *) malloc(sizeof(int) * numParticlesTotal); //array to store weight of particles
 	int * s_x = (int *) malloc(sizeof(int) * numParticlesTotal); //array to store positions of particles in x dimension
 	int * s_y = (int *) malloc(sizeof(int) * numParticlesTotal); //array to store positions of particles in y dimension
 	double * v_x = (double *) malloc(sizeof(double) * numParticlesTotal); //array to store velocities of particles in x dimesion
 	double * v_y = (double *) malloc(sizeof(double) * numParticlesTotal); //array to store velocities of particles in y dimesion
 	double **f_x = contigArrayGenerator(numParticlesTotal,numParticlesTotal); //matrix to store forces of particles in x dimension
 	double **f_y = contigArrayGenerator(numParticlesTotal,numParticlesTotal); //matrix to store forces of particles in y dimension

 	int imageWidth = std::stoi(argv[7]);
 	int imageHeight = std::stoi(argv[8]);

 	double timeSubStep;

 	int width, height;
 	int particlesToReceive;
 	int particlesPerProcessor = numParticlesTotal/p;
 	

 	int * pointerForOriginalArray;
 	int * particlesToCompute_s_x;
 	int * particlesToCompute_s_y;
 	int * particleWeights;
 	double * particlesToCompute_v_x;
 	double * particlesToCompute_v_y;

	int particlesRemaining = numParticlesTotal%p;	

 	int * pointerForLocalArray;
 	int * pointerForTempArray;
 	int * tempWeights;
 	int * localWeights;
 	int * localArray_s_x;
 	int * localArray_s_y;
 	double * localArray_f_x;
 	double * localArray_f_y;
 	int * tempArray_s_x;
 	int * tempArray_s_y;
 	double * tempArray_f_x;
 	double * tempArray_f_y;

 	int * weights;
 	double * forces_x;
 	double * forces_y;


 	MPI_Status status;
 	int source = status.MPI_SOURCE;

 /***************************** MASTER TASK ***************************/
 	if(my_rank == 0){
 		for(i = 0; i < numParticlesTotal; i++){
 			if(numParticlesLight > 0){
 				w[i] = drand48() * (massLightMax-massLightMin+1) + 1;
 				printf("Particle L initial weight: %d\n", w[i]);
 				s_x[i] = drand48() * imageWidth;
 				s_y[i] = drand48() * imageHeight;
 				printf("Particle L initial positions: x: %d / y: %d\n", s_x[i],s_y[i]);
 				v_x[i] = drand48() * (velocityLightMax-velocityLightMin+1) + velocityLightMin;
 				v_y[i] = drand48() * (velocityLightMax-velocityLightMin+1) + velocityLightMin;
 				printf("Particle L initial velocities: x: %f / y: %f\n", v_x[i],v_y[i]);
 				numParticlesLight--;
 			} else if(numParticlesMedium > 0){
 				w[i] = drand48() * (massMediumMax-massMediumMin+1) + massMediumMin;
 				printf("Particle M initial weight: %d\n", w[i]);
  				s_x[i] = drand48()*imageWidth;
 				s_y[i] = drand48()*imageHeight;
 				printf("Particle M initial positions: x: %d / y: %d\n", s_x[i],s_y[i]);
 				v_x[i] = drand48() * (velocityMediumMax-velocityMediumMin+1) + velocityMediumMax;
 				v_y[i] = drand48() * (velocityMediumMax-velocityMediumMin+1) + velocityMediumMax;
 				printf("Particle M initial velocities: x: %f / y: %f\n", v_x[i],v_y[i]);
 				numParticlesMedium--;
 			} else if(numParticlesHeavy > 0){
 				w[i] = drand48() * (massHeavyMax-massHeavyMin+1) + massHeavyMin;
 				printf("Particle H initial weight: %d\n", w[i]);
 				s_x[i] = drand48()*imageWidth;
 				s_y[i] = drand48()*imageHeight;
 				printf("Particle H initial positions: x: %d / y: %d\n", s_x[i],s_y[i]);
 				v_x[i] = drand48() * (velocityHeavyMax-velocityHeavyMin+1) + velocityHeavyMin;
 				v_y[i] = drand48() * (velocityHeavyMax-velocityHeavyMin+1) + velocityHeavyMin;
 				printf("Particle H initial velocities: x: %f / y: %f\n", v_x[i],v_y[i]);
 				numParticlesHeavy--;
 			}
 		}
		for (int frameNum = 0; frameNum < numSteps * numSubSteps; frameNum++) {	
			// ************** ALLOCATED FOR MASTER *************** //
			int * masterWeights;
			int * masterArray_s_x;
			double * masterArray_f_x;
			int * masterArray_s_y;
			double * masterArray_f_y;
			int * masterPointerForLocalArray;
			particlesToReceive = (my_rank < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;
			localWeights = (int *) malloc(sizeof(int) * particlesToReceive); 
			localArray_s_x = (int *) malloc(sizeof(int) * particlesToReceive); 
			localArray_f_x = (double *) malloc(sizeof(double) * particlesToReceive); 
			localArray_s_y = (int *) malloc(sizeof(int) * particlesToReceive); 
			localArray_f_y = (double *) malloc(sizeof(double) * particlesToReceive); 
			pointerForLocalArray = (int *) malloc(sizeof(int) * particlesToReceive); 

			int * tempWeights = (int *) malloc(sizeof(int) * particlesToReceive); 
			int * tempArray_s_x = (int *) malloc(sizeof(int) * particlesToReceive); 
			double * tempArray_f_x = (double *) malloc(sizeof(double) * particlesToReceive); 
			int * tempArray_s_y = (int *) malloc(sizeof(int) * particlesToReceive); 
			double * tempArray_f_y = (double *) malloc(sizeof(double) * particlesToReceive); 
			int * pointerForTempArray = (int *) malloc(sizeof(int) * particlesToReceive); 

			// ************** ALLOCATED FOR SLAVES *************** //
			for (int dest = 0; dest < p; dest++){
				/******* STEP 1: ALLOCATE NUMBER OF PARTICLES TO EACH PROCESSOR *******/
 				particlesToReceive = (dest < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;
 				/******* STEP 2: CREATE ARRAYS TO STORE PARTICLE VALUES & LOCATION IN ORIGINAL ARRAY (Particle number) *******/
				particlesToCompute_s_x = (int *) malloc(sizeof(int) * particlesToReceive); 
	 			particlesToCompute_s_y = (int *) malloc(sizeof(int) * particlesToReceive); 
	 			particleWeights = (int *) malloc(sizeof(int) * particlesToReceive); 
				pointerForOriginalArray = (int *) malloc(sizeof(int) * particlesToReceive);				
 				/******* STEP 3: DISTRIBUTE PARTICLES FROM ORIGINAL ARRAYS TO NEW ARRAYS & MARK LOCATION IN ORIGINAL ARRAYS *******/
				m = 0;
				offset = dest;
				// 
				for(i = offset; i < numParticlesTotal; i+=p){
					particlesToCompute_s_x[m] = s_x[i];
					particlesToCompute_s_y[m] = s_y[i];
					particleWeights[m] = w[i];
					pointerForOriginalArray[m] = i;
					m++;
				}
				/******* MASTER VS PROCESSORS *******/
				if (dest == 0) {
					masterWeights = particleWeights;
	 				masterArray_s_x = particlesToCompute_s_x;
					masterArray_f_x = (double *) malloc(sizeof(double) * particlesToReceive);
					masterArray_s_y = particlesToCompute_s_y;
					masterArray_f_y = (double *) malloc(sizeof(double) * particlesToReceive);
					masterPointerForLocalArray = pointerForOriginalArray;
					for(i = 0; i < particlesToReceive; i++){
						tempArray_s_x[i] = masterArray_s_x[i];
						tempArray_s_y[i] = masterArray_s_y[i];
					}
				}
				else {
					MPI_Send(&(particleWeights[0]), particlesToReceive, MPI_INT, dest, 7, MPI_COMM_WORLD);
					MPI_Send(&(particlesToCompute_s_x[0]), particlesToReceive, MPI_INT, dest, 7, MPI_COMM_WORLD);
					MPI_Send(&(particlesToCompute_s_y[0]), particlesToReceive, MPI_INT, dest, 7, MPI_COMM_WORLD);
					MPI_Send(&(pointerForOriginalArray[0]), particlesToReceive, MPI_INT, dest, 7, MPI_COMM_WORLD);
				}
			}
			// /************** RING LOOP WILL GO HERE***************/
			for(int ringNumber = 0; ringNumber < p - 1; ringNumber++){
				printf("My thread number is %d and my loop (ringNumMaster) is %d\n", my_rank, ringNumber);

				//Send to dest AND receive from source
				int nextRank = (my_rank-1+p)%p;;
				/*if(my_rank == p-1){
					nextRank = 0;
				} else {
					nextRank = my_rank + 1;
				} */
				printf("*******p: %d\n", p);
				printf("*******nextRank: %d\n", nextRank);

				int prevRank = (my_rank+1)%p;
				/* if(my_rank == 0){
					prevRank = p-1;
				} else {
					prevRank = my_rank - 1;
				} */


				//Send to dest AND receive from source
				//MPI_Send(&(localWeights[0]), particlesToReceive, MPI_INT, nextRank, 1, MPI_COMM_WORLD);
				//MPI_Recv(&(tempWeights[0]), particlesToReceive, MPI_INT, prevRank, 1, MPI_COMM_WORLD, &status);
				MPI_Sendrecv(&(localWeights[0]), particlesToReceive, MPI_INT, nextRank, 1, &(tempWeights[0]), particlesToReceive, MPI_INT, prevRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_s_x[0]),  particlesToReceive, MPI_INT, nextRank, 2, &(tempArray_s_x[0]), particlesToReceive, MPI_INT, prevRank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_s_y[0]),  particlesToReceive, MPI_INT, nextRank, 3, &(tempArray_s_y[0]), particlesToReceive, MPI_INT, prevRank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_f_x[0]),  particlesToReceive, MPI_DOUBLE, nextRank, 4, &(tempArray_f_x[0]), particlesToReceive, MPI_DOUBLE, prevRank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				//MPI_Sendrecv(&(tempArray_f_y[0]),  particlesToReceive, MPI_DOUBLE, nextRank, 5, &(tempArray_f_y[0]), particlesToReceive, MPI_DOUBLE, prevRank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				//MPI_Sendrecv(&(pointerForLocalArray[0]),  particlesToReceive, MPI_INT, nextRank, 6, &(pointerForTempArray[0]), particlesToReceive, MPI_INT, prevRank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				j = 0;
				for(i = 0; i < particlesToReceive; i++){
				 	while(pointerForLocalArray[i] >= pointerForTempArray[j] && j < particlesToReceive){ //find index where particle number in tempArray is greater than localArray
				 		j++;
				 	}
				 	if(pointerForLocalArray[i] < pointerForTempArray[j]){
				 		masterArray_f_x[i] += computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[i], tempArray_s_y[i], tempWeights[i], 0);
				 		masterArray_f_y[i] += computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[i], tempArray_s_y[i], tempWeights[i], 1);
				 		tempArray_f_x[j] -= computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[i], tempArray_s_y[i], tempWeights[i], 0);
				 		tempArray_f_y[j] -= computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[i], tempArray_s_y[i], tempWeights[i], 1);
				 	}
				 }
			}

/**************************** MASTER RECIEVES TASKS FROM SLAVES **********************/
			for(int dest = 0; dest < p; dest++) {
				particlesToReceive = (dest < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;

				if (dest == 0) {
					weights = masterWeights;
					forces_x = masterArray_f_x;
					forces_y = masterArray_f_y;
					pointerForOriginalArray = masterPointerForLocalArray;
				}
				else {
					MPI_Recv(&(weights[0]), particlesToReceive, MPI_INT, dest, 8, MPI_COMM_WORLD, &status);
					MPI_Recv(&(forces_x[0]), particlesToReceive, MPI_DOUBLE, dest, 8, MPI_COMM_WORLD, &status);
					MPI_Recv(&(forces_y[0]), particlesToReceive, MPI_DOUBLE, dest, 8, MPI_COMM_WORLD, &status);
					MPI_Recv(&(pointerForOriginalArray[0]), particlesToReceive, MPI_INT, dest, 8, MPI_COMM_WORLD, &status);
				}
				
				for(j = 0; j < particlesToReceive; j++) {
					v_x[pointerForOriginalArray[j]] += timeSubSteps * forces_x[j]/weights[j];
					v_y[pointerForOriginalArray[j]] += timeSubSteps * forces_y[j]/weights[j];
				}

				for(j = 0; j < particlesToReceive; j++) {
					s_x[pointerForOriginalArray[j]] += timeSubSteps * v_x[pointerForOriginalArray[j]];
					s_y[pointerForOriginalArray[j]] += timeSubSteps * v_y[pointerForOriginalArray[j]];
					printf("!!!!!!!I AM FROM DEST %d\n", dest);
					printf("Next positions, x: %d, y: %d\n",s_x[pointerForOriginalArray[j]],s_y[pointerForOriginalArray[j]]);
				}
			}

			unsigned char* image = (unsigned char *) calloc(3*imageWidth*imageHeight, sizeof(unsigned char));
			// distribute particle colours at given position to array to create image
			printArray(w,numParticlesTotal);
			for(i = 0; i < numParticlesTotal; i++){
				int index = (s_y[i] * imageWidth + s_x[i])*3;
				printf("S_X: %d ", s_x[i]);
				printf("\n");
				printf("S_Y: %d ", s_y[i]);
				printf("\n");
				printf("INDEX: %d ", index);
				printf("\n");
				if (index < (sizeof(unsigned char) *3*imageWidth*imageHeight) && index >= 0){
					if(w[i] >= 1 && w[i] <= 5){
						image[index] = 68;
						image[index+1] = 214;
						image[index+2] = 44;
					} else if(w[i] >= massMediumMin && w[i] <= massMediumMax){
						image[index] = 206;
						image[index+1] = 0;
						image[index+2] = 86;

					} else{
						image[index] = 116;
						image[index+1] = 209;
						image[index+2] = 234;
					}
				}
			}
			// Make sure LOGIC HERE IS SOUND
			if (frameNum % numSubSteps == 0) {
				int frameOut = frameNum / numSubSteps;
				char result[50];
				sprintf(result, "%d", frameOut);
				strcat(result, argv[9]);
				saveBMP(result,image, imageWidth, imageHeight);
			}
			free(image);
		}
		
	}



 /*************************** SLAVE TASKS **********************************/
	if(my_rank > 0){
 		particlesToReceive = (my_rank < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;
 		printf("My thread number is %d and my SlaveParticlesToReceive is %d\n", my_rank,particlesToReceive);
		localWeights = (int *) malloc(sizeof(int) * particlesToReceive); 
		localArray_s_x = (int *) malloc(sizeof(int) * particlesToReceive); 
		localArray_f_x = (double *) malloc(sizeof(double) * particlesToReceive); 
		localArray_s_y = (int *) malloc(sizeof(int) * particlesToReceive); 
		localArray_f_y = (double *) malloc(sizeof(double) * particlesToReceive); 
		pointerForLocalArray = (int *) malloc(sizeof(int) * particlesToReceive); 

		tempWeights = (int *) malloc(sizeof(int) * particlesToReceive);
		tempArray_s_x = (int *) malloc(sizeof(int) * particlesToReceive); 
		tempArray_f_x = (double *) malloc(sizeof(double) * particlesToReceive); 
		tempArray_s_y = (int *) malloc(sizeof(int) * particlesToReceive); 
		tempArray_f_y = (double *) malloc(sizeof(double) * particlesToReceive);
		pointerForTempArray = (int *) malloc(sizeof(int) * particlesToReceive); 

 	/******* Recieve particles from MASTER *******/
		for (int numFrames = 0; numFrames < numSteps * numSubSteps; numFrames++){
			printf("My thread number is %d and my loop (slaveFrame) is %d\n", my_rank,numFrames);
			MPI_Recv(&(localWeights[0]), particlesToReceive, MPI_INT, 0, 7, MPI_COMM_WORLD, &status);
			MPI_Recv(&(localArray_s_x[0]), particlesToReceive, MPI_INT, 0, 7, MPI_COMM_WORLD, &status);
			MPI_Recv(&(localArray_s_y[0]), particlesToReceive, MPI_INT, 0, 7, MPI_COMM_WORLD, &status);
			MPI_Recv(&(pointerForLocalArray[0]), particlesToReceive, MPI_INT, 0, 7, MPI_COMM_WORLD, &status);
			printf("SLAVE RECIEVES WEIGHTS\n");
			printArray(localWeights, particlesToReceive);
			for(i = 0; i < particlesToReceive; i++){
				tempArray_s_x[i] = localArray_s_x[i];
				tempArray_s_y[i] = localArray_s_y[i];
				printf("The value in tempArray_s_x is: %d\n", tempArray_s_x[i]);
				printf("The value in tempArray_s_y is: %d\n", tempArray_s_y[i]);
			}
			// RING LOOP GOES HERE
			for(int ringNumber = 0; ringNumber < p - 1; ringNumber++){
				printf("My thread number is %d and my loop (slaveRingNumber) is %d\n", my_rank,ringNumber);
				/******* Send & Recieve particles from another SLAVE *******/
				int nextRank = (my_rank-1+p)%p;
				//if(my_rank == p-1){
					//nextRank = 1;
				//} else {
					//nextRank = my_rank + 1;
				//}

				printf("*******p: %d\n", p);
				printf("*******nextRank: %d\n", nextRank);
				int prevRank = (my_rank+1)%p;
				//if(my_rank == 1){
					//prevRank = p-1;
				//} else {
					//prevRank = my_rank - 1;
				//}

				//CHANGE BASED ON EVEN OR ODD
				//IF EVEN
				// if(my_rank % 2 == 0){
				// 	MPI_Send(&(localWeights[0]), particlesToReceive, MPI_INT, nextRank, 1, MPI_COMM_WORLD);
				// 	MPI_Recv(&(tempWeights[0]), particlesToReceive, MPI_INT, prevRank, 1, MPI_COMM_WORLD, &status);
				// }
				// else {
				// 	MPI_Recv(&(tempWeights[0]), particlesToReceive, MPI_INT, prevRank, 1, MPI_COMM_WORLD, &status);
				// 	MPI_Send(&(localWeights[0]), particlesToReceive, MPI_INT, nextRank, 1, MPI_COMM_WORLD);
				// }

				MPI_Sendrecv(&(localWeights[0]), particlesToReceive, MPI_INT, nextRank, 1, &(tempWeights[0]), particlesToReceive, MPI_INT, prevRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_s_x[0]),  particlesToReceive, MPI_INT, nextRank, 2, &(tempArray_s_x[0]), particlesToReceive, MPI_INT, prevRank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_s_y[0]),  particlesToReceive, MPI_INT, nextRank, 3, &(tempArray_s_y[0]), particlesToReceive, MPI_INT, prevRank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_f_x[0]),  particlesToReceive, MPI_DOUBLE, nextRank, 4, &(tempArray_f_x[0]), particlesToReceive, MPI_DOUBLE, prevRank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				//MPI_Sendrecv(&(tempArray_f_y[0]),  particlesToReceive, MPI_DOUBLE, nextRank, 5, &(tempArray_f_y[0]), particlesToReceive, MPI_DOUBLE, prevRank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				//MPI_Sendrecv(&(pointerForLocalArray[0]),  particlesToReceive, MPI_INT, nextRank, 6, &(pointerForTempArray[0]), particlesToReceive, MPI_INT, prevRank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				
				//Calculate forces
				j = 0;
				for(i = 0; i < particlesToReceive; i++){
				 	while(pointerForLocalArray[i] >= pointerForTempArray[j] && j < particlesToReceive) { //find index where particle number in tempArray is greater than localArray
				 		j++;
				 	}
				 	if(pointerForLocalArray[i] < pointerForTempArray[j]){
				 		printf("Local particle at position %d is interacting with temp particle at position %d\n", localArray_s_x[i], tempArray_s_x[i]);
				 		localArray_f_x[i] += computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[i], tempArray_s_y[i], tempWeights[i], 0);
				 		localArray_f_y[i] += computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[i], tempArray_s_y[i], tempWeights[i], 1);
				 		tempArray_f_x[j] -= computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[i], tempArray_s_y[i], tempWeights[i], 0);
				 		tempArray_f_y[j] -= computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[i], tempArray_s_y[i], tempWeights[i], 1);
				 	}
				}
			}
		// FINAL SEND GOES HERE
		MPI_Send(&(localWeights[0]), particlesToReceive, MPI_INT, 0, 8, MPI_COMM_WORLD);
		MPI_Send(&(localArray_f_x[0]), particlesToReceive, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD);
		MPI_Send(&(localArray_f_y[0]), particlesToReceive, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD);
		MPI_Send(&(pointerForLocalArray[0]), particlesToReceive, MPI_INT, 0, 8, MPI_COMM_WORLD);
		}
	}
	MPI_Finalize();
	return 0;
}

