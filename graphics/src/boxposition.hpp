#ifndef BOXPOSITION_HPP
#define BOXPOSITION_HPP

#include "glmath.hpp"
#include <vector>

    /* 3D objects positions: */
    //glm::vec3 boxPositions[] = 
    std::vector <glm::vec3> boxPositions
    {
        glm::vec3( 0.0f,  5.0f,  0.0f), 
        glm::vec3( 2.0f,  -5.0f, -1.0f) 
    };

    /* the number of objects: */
    int boxesNum = 2;

    /* transformation (scaling) of objects: */
    //glm::vec3 boxesScaling[] = 
    std::vector <glm::vec3> boxesScaling
    {
        glm::vec3( 3.0f,  3.0f,  3.0f),
        glm::vec3( 2.0f,  2.0f,  2.0f)
    };

    /* rotation angle of each object and mode of rotation,
       (1.0 - continuous rotation, 0.0 - one time rotation): */
    //float boxesRotAngle[] = 
    std::vector <float> boxesRotAngle
    {
        20.0f, 0.0f,
        20.0f, 1.0f
    };

    /* rotation axcess of each object: */
    //glm::vec3 boxesRotation[] = 
    std::vector <glm::vec3> boxesRotation
    {
        glm::vec3( 1.0f, 0.3f, 0.5f),
        glm::vec3( 1.0f, 0.3f, 0.5f)
    };

    /* RGB colour of each object: */
    //glm::vec3 boxesColour[] = 
    std::vector <glm::vec3> boxesColour
    {
        glm::vec3( 0.10f, 0.3f, 0.5f),
        glm::vec3( 0.50f, 0.3f, 0.1f)
    };

#endif