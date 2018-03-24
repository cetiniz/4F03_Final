#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <mpi.h>
#include <math.h>

#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#define epsilon 0.000000000000000222

double **contigArrayGenerator(int row, int col){
	double **contigarray = (double **)malloc(row*sizeof(double));
	double *pointer = (double *)malloc(row*col*sizeof(double));
	for(int i = 0; i < row; i++){
		contigarray[i] = &(pointer[col*i]);
	}
	return contigarray;
}

int main(int argc, char* argv[]){
	
	if( argc != 10){
		printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrex\n", argv[0]);
	}

	MPI_Init(&argc,&argv);

	int p, my_rank;
	int i,j,m,offset;

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	//variables
	int numParticlesLight = std::stoi(argv[1]);
	int numParticleMedium = std::stoi(argv[2]);
	int numParticleHeavy = std::stoi(argv[3]);
	int numParticlesTotal = numParticlesLight + numParticlesMedium + numParticlesHeavy; //total number of particles is sum of light, medium, heavy particle numbers

 	int * w = (int *) malloc(sizeof(int) * numParticlesTotal); //array to store weight of particles
 	double * s_x = (double *) malloc(sizeof(double) * numParticlesTotal); //array to store positions of particles in x dimension
 	double * s_y = (double *) malloc(sizeof(double) * numParticlesTotal); //array to store positions of particles in y dimension
 	double * v_x = (double *) malloc(sizeof(double) * numParticlesTotal); //array to store velocities of particles in x dimesion
 	double * v_y = (double *) malloc(sizeof(double) * numParticlesTotal); //array to store velocities of particles in y dimesion
 	double **f_x = contigArrayGenerator(numParticlesTotal,numParticlesTotal); //matrix to store forces of particles in x dimension
 	double **f_y = contigArrayGenerator(numParticlesTotal,numParticlesTotal); //matrix to store forces of particles in y dimension

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
 	int * particlesToCompute_v_x;
 	int * particlesToCompute_v_y;

 	unsigned char* image;

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
 		} else if(i >= numParticlesLight && i < (numParticlesLight+numParticleMedium)){
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


 	for (dest = 0; dest < p; dest++){

 		/******* STEP 1: ALLOCATE NUMBER OF ROWS TO EACH PROCESSOR *******/
 		particlesToReceive = (dest < particlesRemaining) ? particlesPerProcessor+1 : particlesPerProcessor;
 		particlesToCompute_s_x = (int *) malloc(sizeof(int) * particlesToReceive); 
 		particlesToCompute_s_y = (int *) malloc(sizeof(int) * particlesToReceive); 
 		pointerForOriginalArray = (int *) malloc(sizeof(int) * particlesToReceive); //contains pointers that store location of original matrix location (when the Master gathers everything back at the end)

 		m=0;
      	offset = dest;
      	for(i = offset; i < row; i+=p){
          	particlesToCompute_s_x[m] = s_x[i];
          	particlesToCompute_s_y[m] = s_y[i];
        	pointerForOriginalArray[m] = i;
        	m++;
      	}

 	}

 	 if (dest > 0){ //for all other destinations than master processor, send array containing values to compute, pointer array, matrix B
      MPI_Send(&(particlesToCompute_s_x[0]), particlesToReceive, MPI_INT, dest, 0, MPI_COMM_WORLD);
      MPI_Send(&(particlesToCompute_s_y[0]), particlesToReceive, MPI_INT, dest, 0, MPI_COMM_WORLD);
      MPI_Send(&(pointerForOriginalArray[0]), particlesToReceive, MPI_INT, dest, 0, MPI_COMM_WORLD);
    }

 	saveBMP(argv[9], image, width, height);
 }

 else{

 }

 free(image);

 MPI_Finalize();
 return 0;
}