#ifndef LINEPOSITION_HPP
#define LINEPOSITION_HPP

#include "glmath.hpp"
#include <vector>

    /* 3D objects positions: */
    //glm::vec3 linePositions[] = 
    std::vector <glm::vec3> linePositions
    {
        glm::vec3( 0.0f,  5.0f,  0.0f), glm::vec3( 0.0f,  -5.0f, 0.0f),
        glm::vec3( 5.0f,  0.0f,  0.0f), glm::vec3( -5.0f,  0.0f, 0.0f) 
    };

    /* the number of objects: */
    int linesNum = 2;

    /* transformation (scaling) of objects: */
    //glm::vec3 linesScaling[] = 
    std::vector <glm::vec3> linesScaling
    {
        glm::vec3( 5.0f, 0.0f, 0.0f),
        glm::vec3( 7.0f, 0.0f, 0.0f)
    };

    /* rotation angle of each object and mode of rotation,
       (1.0 - continuous rotation, 0.0 - one time rotation): */
    //float linesRotAngle[0];
    std::vector <float> linesRotAngle;

    /* rotation axcess of each object: */
    //glm::vec3 linesRotation[0];
    std::vector <glm::vec3> linesRotation;

    /* RGB colour of each object: */
    //glm::vec3 linesColour[] = 
    std::vector <glm::vec3> linesColour
    {
        glm::vec3( 0.2f, 1.0f, 0.2f),
        glm::vec3( 0.5f, 0.50f, 0.5f)
    };

#endif