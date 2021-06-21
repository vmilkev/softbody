#ifndef PLANEPOSITION_HPP
#define PLANEPOSITION_HPP

#include "glmath.hpp"
#include <vector>

    /* 3D objects positions: */
    //glm::vec3 planePositions[] = 
    std::vector <glm::vec3> planePositions
    {
        glm::vec3( 0.0f,  10.0f,  0.0f), 
        glm::vec3( 0.0f,  -10.0f, 0.0f),
        glm::vec3( -10.0f,  0.0f, 0.0f),
        glm::vec3( 10.0f,  0.0f, 0.0f) 
    };

    /* the number of objects: */
    int planesNum = 4;

    /* transformation (scaling) of objects: */
    //glm::vec3 planesScaling[] = 
    std::vector <glm::vec3> planesScaling
    {
        glm::vec3( 3.0f,  3.0f,  3.0f),
        glm::vec3( 3.0f,  3.0f,  3.0f),
        glm::vec3( 6.0f,  6.0f,  6.0f),
        glm::vec3( 6.0f,  6.0f,  6.0f)
    };

    /* rotation angle of each object and mode of rotation,
       (1.0 - continuous rotation, 0.0 - one time rotation): */
    //float planesRotAngle[] = 
    std::vector <float> planesRotAngle
    {
        0.0f, 0.0f,
        0.0f, 1.0f,
        90.0f, 0.0f,
        90.0f, 0.0f
    };

    /* rotation axcess of each object: */
    //glm::vec3 planesRotation[] = 
    std::vector <glm::vec3> planesRotation
    {
        glm::vec3( 1.0f, 0.3f, 0.5f),
        glm::vec3( 1.0f, 0.3f, 0.5f),
        glm::vec3( 0.0f, 0.0f, 1.0f),
        glm::vec3( 0.0f, 0.0f, 1.0f)
    };

    /* RGB colour of each object: */
    //glm::vec3 planesColour[] = 
    std::vector <glm::vec3> planesColour
    {
        glm::vec3( 0.50f, 0.5f, 0.2f),
        glm::vec3( 0.20f, 0.7f, 0.5f),
        glm::vec3( 0.90f, 0.7f, 0.5f),
        glm::vec3( 0.90f, 0.7f, 0.5f)
    };

#endif