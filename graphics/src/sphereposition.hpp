#ifndef SPHEREPOSITION_HPP
#define SPHEREPOSITION_HPP

#include "glmath.hpp"
#include <vector>

    /* 3D objects positions: */
    //glm::vec3 spherePositions[] = 
    std::vector <glm::vec3> spherePositions
    {
        glm::vec3( 1.0f,  1.0f,  1.0f),
        glm::vec3( 1.0f,  2.0f,  1.0f),
        glm::vec3( 2.0f,  1.0f,  1.0f),
        glm::vec3( 2.0f,  2.0f,  1.0f)
    };

    /* the number of objects: */
    int spheresNum = 4;

    /* transformation (scaling) of objects: */
    //glm::vec3 spheresScaling[] = 
    std::vector <glm::vec3> spheresScaling
    {
        glm::vec3( 2.0f,  2.0f,  2.0f),
        glm::vec3( 2.0f,  2.0f,  2.0f),
        glm::vec3( 2.0f,  2.0f,  2.0f),
        glm::vec3( 2.0f,  2.0f,  2.0f)
    };

    /* rotation angle of each object and mode of rotation,
       (1.0 - continuous rotation, 0.0 - one time rotation): */
    //float spheresRotAngle[] = 
    std::vector <float> spheresRotAngle
    {
        0.0f, 0.0f,
        0.0f, 0.0f,
        0.0f, 0.0f,
        0.0f, 0.0f
    };

    /* rotation axcess of each object: */
    //glm::vec3 spheresRotation[] = 
    std::vector <glm::vec3> spheresRotation
    {
        glm::vec3( 0.4f, 0.3f, 0.5f),
        glm::vec3( 0.4f, 0.3f, 0.5f),
        glm::vec3( 0.4f, 0.3f, 0.5f),
        glm::vec3( 0.4f, 0.3f, 0.5f)
    };

    /* RGB colour of each object: */
    //glm::vec3 spheresColour[] = 
    std::vector <glm::vec3> spheresColour
    {
        glm::vec3( 0.50f, 0.0f, 0.2f),
        glm::vec3( 0.70f, 0.3f, 0.1f),
        glm::vec3( 0.10f, 0.7f, 0.5f),
        glm::vec3( 0.30f, 0.1f, 0.5f)
    };

#endif