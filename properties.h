#ifndef properties_h
#define properties_h

#include "vector3d.h"

//light particles are the fastest
int velocityLightMin = 11;
int velocityLightMax = 15;

int velocityMediumMin = 6;
int velocityMediumMax = 10;

//heavy particles are the slowest
int velocityHeavyMin = 1;
int velocityHeavyMax = 5;

//mass
int massLightMin = 1;
int massLightMax = 5;

int massMediumMin = 6;
int massMediumMax = 10;

int massHeavyMin = 11;
int massHeavyMax = 15;


//colours
vec3 colourLight = vec3(0,0,1);
vec3 colourMedium = vec3(0,1,0);
vec3 colourHeavy = vec3(1,0,0);

#endif