/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include <stdio.h>
#include <stdlib.h>
#include "common_math.h"
#include "trajectories.h"
#include "common_functions.h"
#include "generate_3D_shapes.h"


bool ballsOnCircle(dataType ** image3DPtr, dataType radius, size_t xDim, size_t yDim,
	size_t zDim, dataType xCordctr, dataType yCordctr, dataType bRadius, dataType fillValue)
{
	if (image3DPtr == NULL)
		return false;

	dataType t;
	const dataType pi = 2.0 * (dataType)M_PI;
	const dataType stepSize = (dataType)0.06;
	Point3D ballCenter;

	for (t = 0; t <= pi; t += stepSize)
	{
		ballCenter.x = (dataType)(xCordctr + radius * cos(t));
		ballCenter.y = (dataType)(yCordctr + radius * sin(t));
		ballCenter.z = (dataType)(zDim / 2);
		fillBall3D(image3DPtr, zDim, xDim, yDim, bRadius, ballCenter, fillValue);
	}
	return true;
}