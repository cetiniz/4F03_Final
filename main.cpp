#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <mpi.h>
#include <math.h>

#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#include <omp.h>

#define _XOPEN_SOURCE
#define epsilon 0.000000000000000222
#define g 1.5 //gravitational constant

//Function to make contiguous array
double **contigArrayGenerator(int row, int col){
	double **contigarray = (double **)malloc(row*sizeof(double));
	double *pointer = (double *)malloc(row*col*sizeof(double));
	for(int i = 0; i < row; i++){
		contigarray[i] = &(pointer[col*i]);
	}
	return contigarray;
}

//function to calculate force between two particles
double computeForce(double x1, double y1, int w1, double x2, double y2, int w2, int axisToCompute, bool isTrue = false){
	//if particles have similar x or y values, don't compute since this skews velocity by A LOT
	if ((abs(x1-x2) < 3 &&  abs(y1-y2) < 3) || (abs(x1-x2) > 100000 &&  abs(y1-y2) > 100000)) {
		return 0.0;
	}
	double force = 0;
	double dx = x1-x2; // compute position change in x axis
	double dy = y1-y2; // compute position change in y axis

	double distance = sqrt(dx*dx+dy*dy); // compute distance between particles
	double distance3 = distance*distance*distance; // compute distance cubed between particles
	/* COMPUTE FORCE */
	if(axisToCompute == 0){ //force calculation for x axis
		force = g*w1*w2/distance3*dx;
	} else if (axisToCompute == 1){ //force calculation for y axis
		force = g*w1*w2/distance3*dy;
	}
	return force;
}

//function to print array of ints
void printArray(int * arr, int length){
	for (int i = 0; i < length; i++) {
		printf("%d ",arr[i]);
	}
	printf("\n");
}
//function to print array of doubles
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
	omp_set_dynamic(0);

	/******* INITIALIZE GLOBAL VARIABLES *******/
	//initializing command line arguments:
	int numParticlesLight = std::stoi(argv[1]);
	int numParticlesMedium = std::stoi(argv[2]);
	int numParticlesHeavy = std::stoi(argv[3]);
	int numSteps = std::stoi(argv[4]);
	int numSubSteps = std::stoi(argv[5]);
	double timeSubSteps = std::stof(argv[6]);
	int imageWidth = std::stoi(argv[7]);
 	int imageHeight = std::stoi(argv[8]);

	int numParticlesTotal = numParticlesLight + numParticlesMedium + numParticlesHeavy; //total number of particles is sum of light, medium, heavy particle numbers
	int particlesToReceive;
 	int particlesPerProcessor = numParticlesTotal/p;
 	int particlesRemaining = numParticlesTotal%p; //leftover particles that need to be distributed to avaiable processors

 	//master arrays to track all particle weights, positions, and velocities
 	int * w = (int *) malloc(sizeof(int) * numParticlesTotal); //array to store weight of particles
 	double * s_x = (double *) malloc(sizeof(double) * numParticlesTotal); //array to store positions of particles in x dimension
 	double * s_y = (double *) malloc(sizeof(double) * numParticlesTotal); //array to store positions of particles in y dimension
 	double * v_x = (double *) malloc(sizeof(double) * numParticlesTotal); //array to store velocities of particles in x dimesion
 	double * v_y = (double *) malloc(sizeof(double) * numParticlesTotal); //array to store velocities of particles in y dimesion

 	int * pointerForOriginalArray; //array to track particle number
 	double * particlesToCompute_s_x;
 	double * particlesToCompute_s_y;
 	int * particleWeights;

 	//variables for ring pass communication
 	int * pointerForLocalArray;
 	int * pointerForTempArray;
 	int * tempWeights;
 	int * localWeights;
 	double * localArray_s_x;
 	double * localArray_s_y;
 	double * localArray_f_x;
 	double * localArray_f_y;
 	double * tempArray_s_x;
 	double * tempArray_s_y;
 	double * tempArray_f_x;
 	double * tempArray_f_y;
 	int * weights;
 	double * forces_x;
 	double * forces_y;

 	//variables to measure execution time
 	double minTime, maxTime, avgTime;
 	int counter;

 	MPI_Status status;

 /***************************** MASTER TASK ***************************/
 	if(my_rank == 0){
 		/******* INITIALIZE PARTICLE WEIGHTS, POSITIONS, and VELOCITIES IN MASTER ARRAYS *******/
 		#pragma omp for
 		for(i = 0; i < numParticlesTotal; i++){
 			//initialize light particles
 			if(numParticlesLight > 0){
 				w[i] = drand48() * (massLightMax-massLightMin+1) + massLightMin; //weight and velocity ranges can be found in properties.h
 				s_x[i] = drand48() * imageWidth;
 				s_y[i] = drand48() * imageHeight;
 				//randomize velocities in +x, -x, +y, and -y directions to get random distribution
 				if(i%4==0){
 					v_x[i] = drand48() * (velocityLightMax-velocityLightMin+1) + velocityLightMin;
 					v_y[i] = drand48() * (velocityLightMax-velocityLightMin+1) + velocityLightMin;
 				} else if(i%3==0){
 					v_x[i] = -1*(drand48() * (velocityLightMax-velocityLightMin+1) + velocityLightMin);
 					v_y[i] = -1*(drand48() * (velocityLightMax-velocityLightMin+1) + velocityLightMin);
 				} else if(i%2==0){
 					v_x[i] = -1*(drand48() * (velocityLightMax-velocityLightMin+1) + velocityLightMin);
 					v_y[i] = drand48() * (velocityLightMax-velocityLightMin+1) + velocityLightMin;
 				} else{
 					v_x[i] = drand48() * (velocityLightMax-velocityLightMin+1) + velocityLightMin;
 					v_y[i] = -1*(drand48() * (velocityLightMax-velocityLightMin+1) + velocityLightMin);
 				}
 				#pragma omp atomic 
 				numParticlesLight--;
 			} else if(numParticlesMedium > 0){ //initialize medium particles
 				w[i] = drand48() * (massMediumMax-massMediumMin+1) + massMediumMin;
 				s_x[i] = drand48()*imageWidth;
 				s_y[i] = drand48()*imageHeight;
 				//randomize velocities in +x, -x, +y, and -y directions to get random distribution
 				if(i%4==0){
 					v_x[i] = drand48() * (velocityMediumMax-velocityMediumMin+1) + velocityMediumMin;
 					v_y[i] = drand48() * (velocityMediumMax-velocityMediumMin+1) + velocityMediumMin;
 				} else if(i%3==0){
 					v_x[i] = -1*(drand48() * (velocityMediumMax-velocityMediumMin+1) + velocityMediumMin);
 					v_y[i] = -1*(drand48() * (velocityMediumMax-velocityMediumMin+1) + velocityMediumMin);
 				} else if(i%2==0){
 					v_x[i] = -1*(drand48() * (velocityMediumMax-velocityMediumMin+1) + velocityMediumMin);
 					v_y[i] = drand48() * (velocityMediumMax-velocityMediumMin+1) + velocityMediumMin;
 				} else{
 					v_x[i] = drand48() * (velocityMediumMax-velocityMediumMin+1) + velocityMediumMin;
 					v_y[i] = -1*(drand48() * (velocityMediumMax-velocityMediumMin+1) + velocityMediumMin);
 				}
 				#pragma omp atomic 
 				numParticlesMedium--;
 			} else if(numParticlesHeavy > 0){ //initialize heavy particles
 				w[i] = drand48() * (massHeavyMax-massHeavyMin+1) + massHeavyMin;
 				s_x[i] = drand48()*imageWidth;
 				s_y[i] = drand48()*imageHeight;
 				//randomize velocities in +x, -x, +y, and -y directions to get random distribution
 				if(i%4==0){
 					v_x[i] = drand48() * (velocityHeavyMax-velocityHeavyMin+1) + velocityHeavyMin;
 					v_y[i] = drand48() * (velocityHeavyMax-velocityHeavyMin+1) + velocityHeavyMin;
 				} else if(i%3==0){
 					v_x[i] = -1*(drand48() * (velocityHeavyMax-velocityHeavyMin+1) + velocityHeavyMin);
 					v_y[i] = -1*(drand48() * (velocityHeavyMax-velocityHeavyMin+1) + velocityHeavyMin);
 				} else if(i%2==0){
 					v_x[i] = -1*(drand48() * (velocityHeavyMax-velocityHeavyMin+1) + velocityHeavyMin);
 					v_y[i] = drand48() * (velocityHeavyMax-velocityHeavyMin+1) + velocityHeavyMin;
 				} else{
 					v_x[i] = drand48() * (velocityHeavyMax-velocityHeavyMin+1) + velocityHeavyMin;
 					v_y[i] = -1*(drand48() * (velocityHeavyMax-velocityHeavyMin+1) + velocityHeavyMin);
 				}
 				#pragma omp atomic 
 				numParticlesHeavy--;
 			}
 		}
 		int masterParticlesToReceive = (0 < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor; //assign master with the number of particles it will receive
		//initialize master's local arrays for ring pass communication
		int * masterWeights;
		double * masterArray_s_x;
		double * masterArray_f_x = (double *) calloc(masterParticlesToReceive,sizeof(double)); 
		double * masterArray_s_y;
		double * masterArray_f_y = (double *) calloc(masterParticlesToReceive,sizeof(double)); 
		int * masterPointerForLocalArray;			

		//initialize temp arrays that master will store incoming arrays in for ring pass communication
		int * tempWeights = (int *) malloc(sizeof(int) * particlesPerProcessor+1); 
		double * tempArray_s_x = (double *) malloc(sizeof(double) * particlesPerProcessor+1); 
		double * tempArray_f_x = (double *) calloc(particlesPerProcessor+1,sizeof(double)); 
		double * tempArray_s_y = (double *) malloc(sizeof(double) * particlesPerProcessor+1); 
		double * tempArray_f_y = (double *) calloc(particlesPerProcessor+1,sizeof(double)); 
		int * pointerForTempMasterArray = (int *) malloc(sizeof(int) * particlesPerProcessor+1); 

		//intialize timer variables
		double start, end;
		minTime = 10000000;
 		maxTime = 0;
 		avgTime = 0;
 		counter = 0;
  			
 		// ************** FOR-LOOP GOING THROUGH ALL SIMULATION STEPS (# frames * # substeps) *************** //
		for (int frameNum = 0; frameNum < (numSteps * numSubSteps); frameNum++) {	
			start = MPI_Wtime(); //start timer

			// ************** INITIALIZE VARIABLES & SET ARRAY VALUES FOR MASTER & SLAVES *************** //
			for (int dest = 0; dest < p; dest++){
				if (dest == 0) {
					particlesToReceive = (dest < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;
					masterArray_s_x = (double *) malloc(sizeof(double) * particlesToReceive); 
		 			masterArray_s_y = (double *) malloc(sizeof(double) * particlesToReceive); 
		 			masterWeights = (int *) malloc(sizeof(int) * particlesToReceive); 
					masterPointerForLocalArray = (int *) malloc(sizeof(int) * particlesToReceive);	

					m = 0;
					offset = dest;
					
					//cycle through master array and set master local & temp arrays according to position determined by cyclic distribution & number of particles that are expected to be receieved
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
				} else { //SLAVES:
					particlesToReceive = (dest < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;
					particlesToCompute_s_x = (double *) malloc(sizeof(double) * particlesToReceive); 
		 			particlesToCompute_s_y = (double *) malloc(sizeof(double) * particlesToReceive); 
		 			particleWeights = (int *) malloc(sizeof(int) * particlesToReceive); 
					pointerForOriginalArray = (int *) malloc(sizeof(int) * particlesToReceive);	

					m = 0;
					offset = dest;
					//cycle through slave local array and set them according to positions determined by cyclic distribution & number of particles that are expected to be receieved
					for(i = offset; i < numParticlesTotal; i+=p){
						particlesToCompute_s_x[m] = s_x[i];
						particlesToCompute_s_y[m] = s_y[i];
						particleWeights[m] = w[i];
						pointerForOriginalArray[m] = i;
						m++;
					}
					//SEND these arrays to slaves!
					MPI_Send(&(particleWeights[0]), particlesToReceive, MPI_INT, dest, 7, MPI_COMM_WORLD);
					MPI_Send(&(particlesToCompute_s_x[0]), particlesToReceive, MPI_DOUBLE, dest, 7, MPI_COMM_WORLD);
					MPI_Send(&(particlesToCompute_s_y[0]), particlesToReceive, MPI_DOUBLE, dest, 7, MPI_COMM_WORLD);
					MPI_Send(&(pointerForOriginalArray[0]), particlesToReceive, MPI_INT, dest, 7, MPI_COMM_WORLD);	

					//free memory
					free(particleWeights);
					free(particlesToCompute_s_x);
					free(particlesToCompute_s_y);
					free(pointerForOriginalArray);
				}
			}
			// /************** RING PASS COMMUNICATION  FOR MASTER ***************/
			int currentParticlesSize;
			int nextParticlesSize;
			for(int ringNumber = 0; ringNumber < (p); ringNumber++){

				particlesToReceive = (my_rank < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;
				
				int nextRank = (my_rank-1+p)%p; //processor that current processor will send temp array to in next ring pass
				int prevRank = (my_rank+1)%p; //processor that current processor will receieve temp array from in next ring pass

				if (ringNumber == 0) {
					currentParticlesSize = particlesToReceive;
				}
				else {
					currentParticlesSize = nextParticlesSize;
				}
				nextParticlesSize = (((nextRank + ringNumber)%p) < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;

				//first, compute forces between particles in the master's local array
				if (ringNumber == 0) {
					for(i = 0; i < particlesToReceive; i++){
					 	for(j = 0; j < currentParticlesSize; j++){ 
					 		if(masterPointerForLocalArray[i] < pointerForTempMasterArray[j]){
					 			//add forces in x & y direction to master's local array for smaller particle index
					 			masterArray_f_x[i] += computeForce(masterArray_s_x[i], masterArray_s_y[i], masterWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 0);
					 			masterArray_f_y[i] += computeForce(masterArray_s_x[i], masterArray_s_y[i], masterWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 1);
					 			//subtract forces in x & y direction to master's local array for larger particle index
					 			tempArray_f_x[j] -= computeForce(masterArray_s_x[i], masterArray_s_y[i], masterWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 0);
					 			tempArray_f_y[j] -= computeForce(masterArray_s_x[i], masterArray_s_y[i], masterWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 1);
					 		}
					 	}
					}	
				}

				//Send current temp array possessed by master to next processor & receive new temp array from previous processor
				MPI_Sendrecv(&(tempWeights[0]), particlesPerProcessor+1, MPI_INT, nextRank, 1, &(tempWeights[0]), particlesPerProcessor+1, MPI_INT, prevRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_s_x[0]),  particlesPerProcessor+1, MPI_DOUBLE, nextRank, 2, &(tempArray_s_x[0]), particlesPerProcessor+1, MPI_DOUBLE, prevRank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_s_y[0]),  particlesPerProcessor+1, MPI_DOUBLE, nextRank, 3, &(tempArray_s_y[0]), particlesPerProcessor+1, MPI_DOUBLE, prevRank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_f_x[0]),  particlesPerProcessor+1, MPI_DOUBLE, nextRank, 4, &(tempArray_f_x[0]), particlesPerProcessor+1, MPI_DOUBLE, prevRank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_f_y[0]),  particlesPerProcessor+1, MPI_DOUBLE, nextRank, 5, &(tempArray_f_y[0]), particlesPerProcessor+1, MPI_DOUBLE, prevRank, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(pointerForTempMasterArray[0]), particlesPerProcessor+1, MPI_INT, nextRank, 6, &(pointerForTempMasterArray[0]), particlesPerProcessor+1, MPI_INT, prevRank, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				
				//if on last step of ring pass communication, add all forces in temp array to local master array
				if (ringNumber == (p-1)){
					for(i = 0; i < particlesToReceive; i++) {
						masterArray_f_x[i] += tempArray_f_x[i];
						masterArray_f_y[i] += tempArray_f_y[i];
					}
				} else { //for all other steps in ring, compute forces between new temp array from MPI_Recv and master's local array
					for(i = 0; i < particlesToReceive; i++){
					 	for(j = 0; j < currentParticlesSize; j++){ 
					 		if(masterPointerForLocalArray[i] < pointerForTempMasterArray[j]){
					 			masterArray_f_x[i] += computeForce(masterArray_s_x[i], masterArray_s_y[i], masterWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 0);
					 			masterArray_f_y[i] += computeForce(masterArray_s_x[i], masterArray_s_y[i], masterWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 1);
					 			tempArray_f_x[j] -= computeForce(masterArray_s_x[i], masterArray_s_y[i], masterWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 0);
					 			tempArray_f_y[j] -= computeForce(masterArray_s_x[i], masterArray_s_y[i], masterWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 1);
					 		}
					 	}
					}
				}
			}

			/**************************** MASTER RECIEVES TASKS FROM SLAVES **********************/
			for(int dest = 0; dest < p; dest++) {
				particlesToReceive = (dest < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;

				//set variables to compute new velocities & positions to arrays that are incoming from end of ring pass communication
				if (dest == 0) { 
					weights = masterWeights;
					forces_x = masterArray_f_x;
					forces_y = masterArray_f_y;
					pointerForOriginalArray = masterPointerForLocalArray;
					
				}
				else {
					weights = (int *) malloc(sizeof(int) * particlesToReceive);
					forces_x = (double *) malloc(sizeof(double) * particlesToReceive);
					forces_y = (double *) malloc(sizeof(double) * particlesToReceive);
					pointerForOriginalArray = (int *) malloc(sizeof(int) * particlesToReceive);
					MPI_Recv(&(weights[0]), particlesToReceive, MPI_INT, dest, 8, MPI_COMM_WORLD, &status);
					MPI_Recv(&(forces_x[0]), particlesToReceive, MPI_DOUBLE, dest, 8, MPI_COMM_WORLD, &status);
					MPI_Recv(&(forces_y[0]), particlesToReceive, MPI_DOUBLE, dest, 8, MPI_COMM_WORLD, &status);
					MPI_Recv(&(pointerForOriginalArray[0]), particlesToReceive, MPI_INT, dest, 8, MPI_COMM_WORLD, &status);
					
				}

				for(j = 0; j < particlesToReceive; j++) {
					// ****If particle is OUT OF BOUNDS, change direction of velocity to REBOUND the particle ****//
					if((s_x[pointerForOriginalArray[j]] >= imageWidth)){
						v_x[pointerForOriginalArray[j]] *= -1.0;
					}
					if ((s_x[pointerForOriginalArray[j]] <= 0)) {
						v_x[pointerForOriginalArray[j]] *= -1.0;
					}
					if((s_y[pointerForOriginalArray[j]] >= imageHeight)){
						v_y[pointerForOriginalArray[j]] *= -1.0;
					}
					if ((s_y[pointerForOriginalArray[j]] <= 0)) {
						v_y[pointerForOriginalArray[j]] *= -1.0;
					}

					// ******** MASTER CALCULATES NEW VELOCITIES AND POSITIONS OF ALL PARTICLES ******** //
					v_x[pointerForOriginalArray[j]] += timeSubSteps * forces_x[j]/weights[j];
					v_y[pointerForOriginalArray[j]] += timeSubSteps * forces_y[j]/weights[j];
					s_x[pointerForOriginalArray[j]] += timeSubSteps * v_x[pointerForOriginalArray[j]];
					s_y[pointerForOriginalArray[j]] += timeSubSteps * v_y[pointerForOriginalArray[j]];
				}
			}

			//initialize image array & clear it to 0 using calloc
			unsigned char* image = (unsigned char *) calloc(3*imageWidth*imageHeight, sizeof(unsigned char));

			// distribute particle colours at given position to array to create image
			for(i = 0; i < numParticlesTotal; i++){
				int x = (int) s_x[i];
				int y = (int) s_y[i];
				int index = (y * imageWidth + x)*3; //index of particle in contiguous array

				if (index < (3*imageWidth*imageHeight) && index >= 0){ //ensure pixel is not out of bounds
						//set particle colours
						if(w[i] >= massLightMin && w[i] <= massLightMax){
							image[index] = 255;
							image[index+1] = 160;
							image[index+2] = 0;
						} else if(w[i] >= massMediumMin && w[i] <= massMediumMax){
							image[index] = 0;
							image[index+1] = 255;
							image[index+2] = 255;

						} else{
							image[index] = 255;
							image[index+1] = 0;
							image[index+2] = 255;
						}
				}
			}

			end = MPI_Wtime(); //end timer
			double time = end-start;

			if(time <= minTime){
				minTime = time; //allocate new min time for substep
			} 

			if(time >= maxTime){
				maxTime = time; //allocate new max time for substep
			}

			avgTime += time; 
			counter++;

			//if we are at substep that should be a frame number, save frame as image!
			if (frameNum % numSubSteps == 0) {
				int frameOut = frameNum / numSubSteps;
				char result[50];
				sprintf(result, "%s_%04i.bmp",argv[9], frameOut);
				saveBMP(result,image, imageWidth, imageHeight);
			}
			free(image);
		}
		printf("%f %f %f\n", minTime, maxTime, avgTime/counter); //print out times
	}

 /*************************** SLAVE TASKS **********************************/
	if(my_rank > 0){
		//initialize slaves's local arrays for ring pass communication & temp arrays to pass on
 		particlesToReceive = (my_rank < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;
		localWeights = (int *) calloc(particlesToReceive,sizeof(int)); 
		localArray_s_x = (double *) calloc(particlesToReceive,sizeof(double));  
		localArray_f_x = (double *) calloc(particlesToReceive,sizeof(double));  
		localArray_s_y = (double *) calloc(particlesToReceive,sizeof(double)); 
		localArray_f_y = (double *) calloc(particlesToReceive,sizeof(double));  
		pointerForLocalArray = (int *) calloc(particlesToReceive,sizeof(int));  

		tempWeights = (int *) calloc(particlesPerProcessor+1,sizeof(int)); 
		tempArray_s_x = (double *) calloc(particlesPerProcessor+1,sizeof(double)); 
		tempArray_f_x = (double *) calloc(particlesPerProcessor+1,sizeof(double)); 
		tempArray_s_y = (double *) calloc(particlesPerProcessor+1,sizeof(double));  
		tempArray_f_y = (double *) calloc(particlesPerProcessor+1,sizeof(double)); 
		pointerForTempArray = (int *) calloc(particlesPerProcessor+1,sizeof(int)); 

 	/******* Recieve particles from MASTER *******/
		for (int numFrames = 0; numFrames < numSteps * numSubSteps; numFrames++){
			MPI_Recv(&(localWeights[0]), particlesToReceive, MPI_INT, 0, 7, MPI_COMM_WORLD, &status);
			MPI_Recv(&(localArray_s_x[0]), particlesToReceive, MPI_DOUBLE, 0, 7, MPI_COMM_WORLD, &status);
			MPI_Recv(&(localArray_s_y[0]), particlesToReceive, MPI_DOUBLE, 0, 7, MPI_COMM_WORLD, &status);
			MPI_Recv(&(pointerForLocalArray[0]), particlesToReceive, MPI_INT, 0, 7, MPI_COMM_WORLD, &status);
			
			//allocate arrays received from master to temp array
			for(i = 0; i < particlesToReceive; i++){
				tempArray_s_x[i] = localArray_s_x[i];
				tempArray_s_y[i] = localArray_s_y[i];
				tempArray_f_x[i] = 0.00;
				tempArray_f_y[i] = 0.00;
				pointerForTempArray[i] = pointerForLocalArray[i];
			}

			int currentParticlesSize;
			int nextParticlesSize;
			// /************** RING PASS COMMUNICATION  FOR SLAVES ***************/
			for(int ringNumber = 0; ringNumber < p; ringNumber++){
				/******* Send & Recieve particles from another SLAVE *******/
				int nextRank = (my_rank-1+p)%p; //processor slave will send temp array to
				int prevRank = (my_rank+1)%p; //processor slave will receive temp array from
				
				if (ringNumber == 0) {
					currentParticlesSize = particlesToReceive;
				}
				else {
					currentParticlesSize = nextParticlesSize;
				}
				nextParticlesSize = (((nextRank + ringNumber)%p) < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;
			
				//first, compute interactions between particles in slave's own array
				if (ringNumber == 0) {
					for(i = 0; i < particlesToReceive; i++){
					 	for(j = 0; j < currentParticlesSize; j++){
					 		if(pointerForLocalArray[i] < pointerForTempArray[j]){
					 			//for smaller particle numbers, add force between 2 interacting particles to local & and subtract from larger particle number in temp array
					 			localArray_f_x[i] += computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 0);
					 			localArray_f_y[i] += computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 1);
					 			tempArray_f_x[j] -= computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 0);
					 			tempArray_f_y[j] -= computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 1);
					 		}
					 	}
					}
				}
				
				//Send current temp array possessed by master to next processor & receive new temp array from previous processor
				MPI_Sendrecv(&(tempWeights[0]), particlesPerProcessor+1, MPI_INT, nextRank, 1, &(tempWeights[0]), particlesPerProcessor+1, MPI_INT, prevRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_s_x[0]),  particlesPerProcessor+1, MPI_DOUBLE, nextRank, 2, &(tempArray_s_x[0]), particlesPerProcessor+1, MPI_DOUBLE, prevRank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_s_y[0]),  particlesPerProcessor+1, MPI_DOUBLE, nextRank, 3, &(tempArray_s_y[0]), particlesPerProcessor+1, MPI_DOUBLE, prevRank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_f_x[0]),  particlesPerProcessor+1, MPI_DOUBLE, nextRank, 4, &(tempArray_f_x[0]), particlesPerProcessor+1, MPI_DOUBLE, prevRank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(tempArray_f_y[0]),  particlesPerProcessor+1, MPI_DOUBLE, nextRank, 5, &(tempArray_f_y[0]), particlesPerProcessor+1, MPI_DOUBLE, prevRank, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&(pointerForTempArray[0]),  particlesPerProcessor+1, MPI_INT, nextRank, 6, &(pointerForTempArray[0]), particlesPerProcessor+1, MPI_INT, prevRank, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				//if on last step of ring pass communication, add all forces in temp array to local master array
				if (ringNumber == (p-1)){
					for(i = 0; i < particlesToReceive; i++) {
						localArray_f_x[i] += tempArray_f_x[i];
						localArray_f_x[i] += tempArray_f_y[i];
					}
				}
				else {  //for all other steps in ring, compute forces between new temp array from MPI_Recv and master's local array
					for(i = 0; i < particlesToReceive; i++){
					 	for(j = 0; j < nextParticlesSize; j++){ 
					 		if(pointerForLocalArray[i] < pointerForTempArray[j]){
					 			localArray_f_x[i] += computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 0);
					 			localArray_f_y[i] += computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 1);
					 			tempArray_f_x[j] += computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 0);
					 			tempArray_f_y[j] += computeForce(localArray_s_x[i], localArray_s_y[i], localWeights[i], tempArray_s_x[j], tempArray_s_y[j], tempWeights[j], 1);
					 		}

					 	}
					}
				}
			}
		// FINAL SEND TO MASTER GOES HERE
		MPI_Send(&(localWeights[0]), particlesToReceive, MPI_INT, 0, 8, MPI_COMM_WORLD);
		MPI_Send(&(localArray_f_x[0]), particlesToReceive, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD);
		MPI_Send(&(localArray_f_y[0]), particlesToReceive, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD);
		MPI_Send(&(pointerForLocalArray[0]), particlesToReceive, MPI_INT, 0, 8, MPI_COMM_WORLD);
		}

		//free all memory
		free(localWeights);
		free(localArray_s_x);
		free(localArray_f_x); 
		free(localArray_s_y);
		free(localArray_f_y);
		free(pointerForLocalArray); 

		free(tempWeights);
		free(tempArray_s_x); 
		free(tempArray_f_x);
		free(tempArray_s_y);
		free(tempArray_f_y);
		free(pointerForTempArray);
	}
	MPI_Finalize();
	return 0; 
}
