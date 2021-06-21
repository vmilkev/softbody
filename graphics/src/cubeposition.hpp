#ifndef CUBEPOSITION_HPP
#define CUBEPOSITION_HPP

#include "glmath.hpp"
#include <vector>

    /* 3D cube positions: */
    //glm::vec3 cubePositions[] = 
    std::vector <glm::vec3> cubePositions
    {
        glm::vec3( 0.0f,  0.0f,  0.0f), 
        glm::vec3( 2.0f,  5.0f, -15.0f), 
        glm::vec3(-1.5f, -2.2f, -2.5f),  
        glm::vec3(-3.8f, -2.0f, -12.3f),  
        glm::vec3( 2.4f, -0.4f, -3.5f),  
        glm::vec3(-1.7f,  3.0f, -7.5f),  
        glm::vec3( 1.3f, -2.0f, -2.5f),  
        glm::vec3( 1.5f,  2.0f, -2.5f), 
        glm::vec3( 1.5f,  0.2f, -1.5f), 
        glm::vec3(-1.3f,  1.0f, -1.5f)  
    };

    /* the number of objects: */
    int cubesNum = 10;

    /* transformation (scaling) of objects: */
    //glm::vec3 cubesScaling[] = 
    std::vector <glm::vec3> cubesScaling
    {
        glm::vec3( 1.0f,  1.0f,  1.0f),
        glm::vec3( 2.0f,  2.0f,  2.0f),
        glm::vec3( 3.0f,  3.0f,  3.0f),
        glm::vec3( 1.0f,  1.0f,  1.0f),
        glm::vec3( 1.0f,  1.0f,  1.0f),
        glm::vec3( 1.0f,  1.0f,  1.0f),
        glm::vec3( 1.0f,  1.0f,  1.0f),
        glm::vec3( 1.0f,  1.0f,  1.0f),
        glm::vec3( 1.0f,  1.0f,  1.0f),
        glm::vec3( 1.0f,  1.0f,  1.0f)
    };

    /* rotation angle of each object and mode of rotation,
       (1.0 - continuous rotation, 0.0 - one time rotation): */
    //float cubesRotAngle[] = 
    std::vector <float> cubesRotAngle
    {
        20.0f, 0.0f,
        20.0f, 0.0f,
        20.0f, 0.0f,
        20.0f, 0.0f,
        20.0f, 0.0f,
        20.0f, 0.0f,
        20.0f, 0.0f,
        20.0f, 0.0f,
        20.0f, 0.0f,
        20.0f, 0.0f
    };

    /* rotation axcess of each object: */
    //glm::vec3 cubesRotation[] = 
    std::vector <glm::vec3> cubesRotation
    {
        glm::vec3( 1.0f, 0.3f, 0.5f),
        glm::vec3( 1.0f, 0.3f, 0.5f),
        glm::vec3( 1.0f, 0.3f, 0.5f),
        glm::vec3( 1.0f, 0.3f, 0.5f),
        glm::vec3( 1.0f, 0.3f, 0.5f),
        glm::vec3( 1.0f, 0.3f, 0.5f),
        glm::vec3( 1.0f,  1.0f,  1.0f),
        glm::vec3( 1.0f,  1.0f,  1.0f),
        glm::vec3( 1.0f,  1.0f,  1.0f),
        glm::vec3( 1.0f,  1.0f,  1.0f)
    };

    /* RGB colour of each object: */
    //glm::vec3 cubesColour[] = 
    std::vector <glm::vec3> cubesColour
    {
        glm::vec3( 1.0f, 0.3f, 0.5f),
        glm::vec3( 1.0f, 0.3f, 0.5f),
        glm::vec3( 1.0f, 0.3f, 0.5f),
        glm::vec3( 1.0f, 0.3f, 0.5f),
        glm::vec3( 1.0f, 0.3f, 0.5f),
        glm::vec3( 1.0f, 0.3f, 0.5f),
        glm::vec3( 1.0f,  1.0f,  1.0f),
        glm::vec3( 1.0f,  1.0f,  1.0f),
        glm::vec3( 1.0f,  1.0f,  1.0f),
        glm::vec3( 1.0f,  1.0f,  1.0f)
    };

#endif