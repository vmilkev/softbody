#ifndef LIGHTPOSITION_HPP
#define LIGHTPOSITION_HPP

#include "glmath.hpp"
#include <vector>

    /* 3D objects positions: */
    glm::vec3 lightPositions[] = {
        glm::vec3( 10.0f,  20.0f,  0.0f)
        //glm::vec3( -10.0f,  -20.0f,  0.0f)
    };

    /* the number of objects: */
    int lightsNum = 1;

    /* transformation (scaling) of objects: */
    glm::vec3 lightsScaling[] = {
        glm::vec3( 1.0f,  1.0f,  1.0f)
    };

    /* rotation angle of each object and mode of rotation,
       (1.0 - continuous rotation, 0.0 - one time rotation): */
    float lightsRotAngle[] = {
        0.0f, 0.0f
    };

    /* rotation axcess of each object: */
    glm::vec3 lightsRotation[] = {
        glm::vec3( 1.0f, 0.3f, 0.5f)
    };

    /* RGB colour of each object: */
    glm::vec3 lightsColour[] = {
        glm::vec3( 1.0f, 1.0f, 1.0f)
    };

#endif