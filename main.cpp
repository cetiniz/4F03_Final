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
	int i,j;

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	//variables
	int numParticlesLight = std::stoi(argv[1]);
	int numParticleMedium = std::stoi(argv[2]);
	int numParticleHeavy = std::stoi(argv[3]);
	int numParticlesTotal = numParticlesLight + numParticlesMedium + numParticlesHeavy; //total number of particles is sum of light, medium, heavy particle numbers

 	int * w = (int *) malloc(sizeof(int) * numParticlesTotal); //array to store weight of particles
 	double * s_x = (double *) malloc(sizeof(double) * numParticlesTotal); //matrix to store positions of particles
 	double * s_y = (double *) malloc(sizeof(double) * numParticlesTotal); //matrix to store positions of particles
 	double * v_x = (double *) malloc(sizeof(double) * numParticlesTotal); //matrix to store velocities of particles
 	double * v_y = (double *) malloc(sizeof(double) * numParticlesTotal);
 	double **f = contigArrayGenerator(numParticlesTotal,numParticlesTotal); //matrix to store forces of particles


 	int imageWidth = std::stoi(argv[7]);
 	int imageHeight = std::stoi(argv[8]);

 	int numSteps = 0;
 	int subSteps = 0;
 	double timeSubStep;

 	int width, height;

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
 		
 	}

 	saveBMP(argv[9], image, width, height);
 }

 else{

 }

 free(image);

 MPI_Finalize();
 return 0;
}