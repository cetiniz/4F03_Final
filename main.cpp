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
//#define g pow(6.673*10, -11)
#define g 0.001

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
	//printf("X1:%d,Y1:%d,W1:%d,X2:%d, Y2:%d, W2: %d\n", x1,y1,w1,x2,y2,w2);
	int dx = x1-x2; // compute position change in x axis
	int dy = y1-y2; // compute position change in y axis
	
	double distance = sqrt(dx*dx+dy*dy); // compute distance between particles
	double distance3 = distance*distance*distance; // compute distance cubed between particles
	
	/* COMPUTE FORCE */
	if(axisToCompute == 0){
		force = -g*w1*w2/distance3*dx;
		//printf("%f\n",0.0000000000667);
	} else if (axisToCompute == 1){
		force = -g*w1*w2/distance3*dy;
	}
	return force;
}
void printArray(int * arr, int length){
	for (int i = 0; i < length; i++) {
		printf("%d ",arr[i]);
	}
	printf("\n");
}

void printArrayD(double * arr, int length){
	for (int i = 0; i < length; i++) {
		printf("%f ",arr[i]);
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

 	int imageWidth = std::stoi(argv[7]);
 	int imageHeight = std::stoi(argv[8]);
	
 	int particlesToReceive;
 	int tempParticlesToReceive;
 	int particlesPerProcessor = numParticlesTotal/p;
 	

 	int * pointerForOriginalArray;
 	int * particlesToCompute_s_x;
 	int * particlesToCompute_s_y;
 	int * particleWeights;

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

 /***************************** MASTER TASK ***************************/
 	if(my_rank == 0){
 		for(i = 0; i < numParticlesTotal; i++){
 			if(numParticlesLight > 0){
 				w[i] = drand48() * (massLightMax-massLightMin+1) + massLightMin;
 				s_x[i] = drand48() * imageWidth;
 				s_y[i] = drand48() * imageHeight;
 				if(i%2 ==0){
 					v_x[i] = drand48() * (velocityLightMax-velocityLightMin+1) + velocityLightMin;
 					v_y[i] = drand48() * (velocityLightMax-velocityLightMin+1) + velocityLightMin;
 				} else{
 					v_x[i] = -1*(drand48() * (velocityLightMax-velocityLightMin+1) + velocityLightMin);
 					v_y[i] = -1*(drand48() * (velocityLightMax-velocityLightMin+1) + velocityLightMin);
 				}
 				numParticlesLight--;
 			} else if(numParticlesMedium > 0){
 				w[i] = drand48() * (massMediumMax-massMediumMin+1) + massMediumMin;
  				s_x[i] = drand48()*imageWidth;
 				s_y[i] = drand48()*imageHeight;
 				if(i%2 ==0){
 					v_x[i] = drand48() * (velocityMediumMax-velocityMediumMin+1) + velocityMediumMin;
 					v_y[i] = drand48() * (velocityMediumMax-velocityMediumMin+1) + velocityMediumMin;
 				} else{
 					v_x[i] = -1*(drand48() * (velocityMediumMax-velocityMediumMin+1) + velocityMediumMin);
 					v_y[i] = -1*(drand48() * (velocityMediumMax-velocityMediumMin+1) + velocityMediumMin);
 				}
 				numParticlesMedium--;
 			} else if(numParticlesHeavy > 0){
 				w[i] = drand48() * (massHeavyMax-massHeavyMin+1) + massHeavyMin;
 				s_x[i] = drand48()*imageWidth;
 				s_y[i] = drand48()*imageHeight;
 				if(i%2 ==0){
 					v_x[i] = drand48() * (velocityHeavyMax-velocityHeavyMin+1) + velocityHeavyMin;
 					v_y[i] = drand48() * (velocityHeavyMax-velocityHeavyMin+1) + velocityHeavyMin;
 				} else{
 					v_x[i] = -1*(drand48() * (velocityHeavyMax-velocityHeavyMin+1) + velocityHeavyMin);
 					v_y[i] = -1*(drand48() * (velocityHeavyMax-velocityHeavyMin+1) + velocityHeavyMin);
 				}
 				numParticlesHeavy--;
 			}
 		}
		for (int frameNum = 0; frameNum < numSteps * numSubSteps; frameNum++) {	
			// ************** ALLOCATED FOR MASTER *************** //
			int masterParticlesToReceive = (0 < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;
			int * masterWeights;
			int * masterArray_s_x;
			double * masterArray_f_x = (double *) calloc(masterParticlesToReceive,sizeof(double)); 
			int * masterArray_s_y;
			double * masterArray_f_y = (double *) calloc(masterParticlesToReceive,sizeof(double)); 
			int * masterPointerForLocalArray;			

			int * tempWeights = (int *) malloc(sizeof(int) * masterParticlesToReceive); 
			int * tempArray_s_x = (int *) malloc(sizeof(int) * masterParticlesToReceive); 
			double * tempArray_f_x = (double *) calloc(masterParticlesToReceive,sizeof(double)); 
			int * tempArray_s_y = (int *) malloc(sizeof(int) * masterParticlesToReceive); 
			double * tempArray_f_y = (double *) calloc(masterParticlesToReceive,sizeof(double)); 
			int * pointerForTempMasterArray = (int *) malloc(sizeof(int) * masterParticlesToReceive); 

			// ************** ALLOCATED FOR SLAVES *************** //
			for (int dest = 0; dest < p; dest++){
				if (dest == 0) {
					particlesToReceive = (dest < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;
					masterArray_s_x = (int *) malloc(sizeof(int) * particlesToReceive); 
		 			masterArray_s_y = (int *) malloc(sizeof(int) * particlesToReceive); 
		 			masterWeights = (int *) malloc(sizeof(int) * particlesToReceive); 
					masterPointerForLocalArray = (int *) malloc(sizeof(int) * particlesToReceive);	

					m = 0;
					offset = dest;
					
					for(i = offset; i < numParticlesTotal; i+=p){
						tempArray_s_x[m] = s_x[i];
						masterArray_s_x[m] = s_x[i];

						tempArray_s_y[m] = s_y[i];
						masterArray_s_y[m] = s_y[i];

						tempWeights[m] = w[i];
						masterWeights[m] = w[i];

						pointerForTempMasterArray[m] = i;
						masterPointerForLocalArray[m] = i;
						m++;
					}
				}

				
				else {
					particlesToReceive = (dest < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;
					particlesToCompute_s_x = (int *) malloc(sizeof(int) * particlesToReceive); 
		 			particlesToCompute_s_y = (int *) malloc(sizeof(int) * particlesToReceive); 
		 			particleWeights = (int *) malloc(sizeof(int) * particlesToReceive); 
					pointerForOriginalArray = (int *) malloc(sizeof(int) * particlesToReceive);	

					m = 0;
					offset = dest;
					
					for(i = offset; i < numParticlesTotal; i+=p){
						particlesToCompute_s_x[m] = s_x[i];
						particlesToCompute_s_y[m] = s_y[i];
						particleWeights[m] = w[i];
						pointerForOriginalArray[m] = i;
						m++;
					}

					MPI_Send(&(particleWeights[0]), particlesToReceive, MPI_INT, dest, 7, MPI_COMM_WORLD);
					MPI_Send(&(particlesToCompute_s_x[0]), particlesToReceive, MPI_INT, dest, 7, MPI_COMM_WORLD);
					MPI_Send(&(particlesToCompute_s_y[0]), particlesToReceive, MPI_INT, dest, 7, MPI_COMM_WORLD);
					MPI_Send(&(pointerForOriginalArray[0]), particlesToReceive, MPI_INT, dest, 7, MPI_COMM_WORLD);		
				}

			}
			// /************** RING LOOP WILL GO HERE***************/
			for(int ringNumber = 0; ringNumber < (p); ringNumber++){				
				//Send to dest AND receive from source
				int nextRank = (my_rank-1+p)%p;
				int prevRank = (my_rank+1)%p;

				
				//Send to dest AND receive from source
				
				/*printf("POINTER(master)\n");
				printArray(masterPointerForLocalArray, particlesToReceive);
				printf("Weight(Master)\n");
				printArray(masterWeights, particlesToReceive);
				printf("FORCE(master)\n");
				printArrayD(masterArray_f_x, particlesToReceive);*/

				//printArray(masterWeights,1);
				if (ringNumber == 0) {
					for(i = 0; i < particlesToReceive; i++){
					 	for(j = i+1; j < particlesToReceive; j++){ //MAKE SURE PARTICLES TO RECEIVE ARE DIFFERENT NUMBERS!!!!!!!!!!!!!
					 		if(masterPointerForLocalArray[i] < pointerForTempMasterArray[j]){
					 			//printf("Local particle at position %d is interacting with temp particle at position %d\n", localArray_s_x[i], tempArray_s_x[i]);
					 			masterArray_f_x[i] += computeForce(masterArray_s_x[i], masterArray_s_y[i], masterWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 0);
					 			masterArray_f_y[i] += computeForce(masterArray_s_x[i], masterArray_s_y[i], masterWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 1);
					 			tempArray_f_x[j] -= computeForce(masterArray_s_x[i], masterArray_s_y[i], masterWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 0);
					 			tempArray_f_y[j] -= computeForce(masterArray_s_x[i], masterArray_s_y[i], masterWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 1);
					 		}
					 	}
					}
				}

				int particledReceived;
				if(ringNumber ==0){
					particledReceived = (my_rank < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;
				} else{
					particledReceived = (prevRank < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;
				}
				int incomingParticles = (nextRank < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;

				//printf("PARTICLES RECIEVED: %d\n", particledReceived);
				//printf("INCOMING PARTICLES: %d\n", incomingParticles);

				MPI_Sendrecv(&(tempWeights[0]), particledReceived, MPI_INT, nextRank, 1, &(tempWeights[0]), incomingParticles, MPI_INT, prevRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_s_x[0]),  particledReceived, MPI_INT, nextRank, 2, &(tempArray_s_x[0]), incomingParticles, MPI_INT, prevRank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_s_y[0]),  particledReceived, MPI_INT, nextRank, 3, &(tempArray_s_y[0]), incomingParticles, MPI_INT, prevRank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_f_x[0]),  particledReceived, MPI_DOUBLE, nextRank, 4, &(tempArray_f_x[0]), incomingParticles, MPI_DOUBLE, prevRank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_f_y[0]),  particledReceived, MPI_DOUBLE, nextRank, 5, &(tempArray_f_y[0]), incomingParticles, MPI_DOUBLE, prevRank, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(pointerForTempMasterArray[0]), particledReceived, MPI_INT, nextRank, 6, &(pointerForTempMasterArray[0]), incomingParticles, MPI_INT, prevRank, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				

				for(i = 0; i < particlesToReceive; i++){
				 	for(j = 0; j < incomingParticles; j++){ //MAKE SURE PARTICLES TO RECEIVE ARE DIFFERENT NUMBERS!!!!!!!!!!!!!
				 		if(masterPointerForLocalArray[i] < pointerForTempMasterArray[j]){
				 			//printf("Local particle at position %d is interacting with temp particle at position %d\n", localArray_s_x[i], tempArray_s_x[i]);
				 			masterArray_f_x[i] += computeForce(masterArray_s_x[i], masterArray_s_y[i], masterWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 0);
				 			masterArray_f_y[i] += computeForce(masterArray_s_x[i], masterArray_s_y[i], masterWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 1);
				 			tempArray_f_x[j] -= computeForce(masterArray_s_x[i], masterArray_s_y[i], masterWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 0);
				 			tempArray_f_y[j] -= computeForce(masterArray_s_x[i], masterArray_s_y[i], masterWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 1);
				 		}
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
					s_x[pointerForOriginalArray[j]] += timeSubSteps * v_x[pointerForOriginalArray[j]];
					s_y[pointerForOriginalArray[j]] += timeSubSteps * v_y[pointerForOriginalArray[j]];
					//printf("!!!!!!!X");
					//printArray(s_x,particlesToReceive);
					//printf("!!!!!!!Y");
					//printArray(s_y,particlesToReceive);
					//printf("!!!!!!!I AM FROM DEST %d\n", dest);
					printf("Next positions, x: %d, y: %d\n",s_x[pointerForOriginalArray[j]],s_y[pointerForOriginalArray[j]]);
				}
			}

			unsigned char* image = (unsigned char *) calloc(3*imageWidth*imageHeight, sizeof(unsigned char));
			// distribute particle colours at given position to array to create image
			for(i = 0; i < numParticlesTotal; i++){
				int index = (s_y[i] * imageWidth + s_x[i])*3;
				/*printf("S_X: %d ", s_x[i]);
				printf("\n");
				printf("S_Y: %d ", s_y[i]);
				printf("\n");
				printf("INDEX: %d ", index);
				printf("\n"); */
				if (index < (sizeof(unsigned char) *3*imageWidth*imageHeight) && index >= 0){
					if(w[i] >= 1 && w[i] <= 5){
						image[index] = 68;
						image[index+1] = 214;
						image[index+2] = 44;
					} else if(w[i] >= massMediumMin && w[i] <= massMediumMax){
						image[index] = 135;
						image[index+1] = 206;
						image[index+2] = 250;

					} else{
						image[index] = 255;
						image[index+1] = 20;
						image[index+2] = 147;
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
 		tempParticlesToReceive = particlesToReceive;
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
			MPI_Recv(&(localWeights[0]), particlesToReceive, MPI_INT, 0, 7, MPI_COMM_WORLD, &status);
			MPI_Recv(&(localArray_s_x[0]), particlesToReceive, MPI_INT, 0, 7, MPI_COMM_WORLD, &status);
			MPI_Recv(&(localArray_s_y[0]), particlesToReceive, MPI_INT, 0, 7, MPI_COMM_WORLD, &status);
			MPI_Recv(&(pointerForLocalArray[0]), particlesToReceive, MPI_INT, 0, 7, MPI_COMM_WORLD, &status);
			
			for(i = 0; i < particlesToReceive; i++){
				tempArray_s_x[i] = localArray_s_x[i];
				tempArray_s_y[i] = localArray_s_y[i];
				tempArray_f_x[i] = 0.00;
				tempArray_f_y[i] = 0.00;
				pointerForTempArray[i] = pointerForLocalArray[i];
			}
			printArray(tempArray_s_x,particlesToReceive);
			// RING LOOP GOES HERE
			for(int ringNumber = 0; ringNumber < p; ringNumber++){
				//printf("My thread number is %d and my loop (slaveRingNumber) is %d\n", my_rank,ringNumber);
				/******* Send & Recieve particles from another SLAVE *******/
				int nextRank = (my_rank-1+p)%p;

				int prevRank = (my_rank+1)%p;

				//Calculate forces
				if (ringNumber == 0) {
					for(i = 0; i < particlesToReceive; i++){
					 	for(j = i+1; j < particlesToReceive; j++){ //MAKE SURE PARTICLES TO RECEIVE ARE DIFFERENT NUMBERS!!!!!!!!!!!!!
					 		if(pointerForLocalArray[i] < pointerForTempArray[j]){
					 			//printf("Local particle at position %d is interacting with temp particle at position %d\n", localArray_s_x[i], tempArray_s_x[i]);
					 			localArray_f_x[i] += computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 0);
					 			localArray_f_y[i] += computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 1);
					 			tempArray_f_x[j] -= computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 0);
					 			tempArray_f_y[j] -= computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 1);
					 		}
					 	}
					}
				}

				int particledReceived;
				if(ringNumber ==0){
					particledReceived = (my_rank < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;
				} else{
					particledReceived = (prevRank < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;
				}
				int incomingParticles = (nextRank < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;

				//printf("PARTICLES RECIEVED: %d\n", particledReceived);
				//printf("INCOMING PARTICLES: %d\n", incomingParticles);

				MPI_Sendrecv(&(tempWeights[0]), particledReceived, MPI_INT, nextRank, 1, &(tempWeights[0]), incomingParticles, MPI_INT, prevRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_s_x[0]),  particledReceived, MPI_INT, nextRank, 2, &(tempArray_s_x[0]), incomingParticles, MPI_INT, prevRank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_s_y[0]),  particledReceived, MPI_INT, nextRank, 3, &(tempArray_s_y[0]), incomingParticles, MPI_INT, prevRank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_f_x[0]),  particledReceived, MPI_DOUBLE, nextRank, 4, &(tempArray_f_x[0]), incomingParticles, MPI_DOUBLE, prevRank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_f_y[0]),  particledReceived, MPI_DOUBLE, nextRank, 5, &(tempArray_f_y[0]), incomingParticles, MPI_DOUBLE, prevRank, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(pointerForTempArray[0]),  particledReceived, MPI_INT, nextRank, 6, &(pointerForTempArray[0]), incomingParticles, MPI_INT, prevRank, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				

				for(i = 0; i < particlesToReceive; i++){
				 	for(j = 0; j < incomingParticles; j++){ //MAKE SURE PARTICLES TO RECEIVE ARE DIFFERENT NUMBERS!!!!!!!!!!!!!
				 		//printf("localPointer(slave): %d  ;   remote Pointer(slave):%d\n",pointerForLocalArray[i],pointerForTempArray[j]);
				 		if(pointerForLocalArray[i] < pointerForTempArray[j]){
				 			localArray_f_x[i] += computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 0);
				 			localArray_f_y[i] += computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 1);
				 			tempArray_f_x[j] -= computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 0);
				 			tempArray_f_y[j] -= computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 1);
				 		}
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
