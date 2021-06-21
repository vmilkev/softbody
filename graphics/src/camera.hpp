#ifndef CAMERA_HPP
#define CAMERA_HPP

#include <glad/glad.h>
#include "glmath.hpp"

#include <vector>
#include <cmath>
#include <iostream>


namespace oglu
{

    // Defines several possible options for camera movement. Used as abstraction to stay away from window-system specific input methods
    enum Camera_Movement {
        FORWARD,
        BACKWARD,
        LEFT,
        RIGHT,
        UP,
        DOWN,
        CENTER,
        VIEWIN,
        VIEWOUT
    };

    // Default camera values
    const float YAW         = -90.0f;
    const float PITCH       =  0.0f;
    const float SPEED       =  3.5f;
    const float SENSITIVITY =  0.1f;
    const float ZOOM        =  35.0f;


    // An abstract camera class that processes input and calculates the corresponding Euler Angles, Vectors and Matrices for use in OpenGL
    class camera
    {
    public:
        // Camera Attributes
        glm::vec3 cameraPosition;
        glm::vec3 cameraPosition0;
        glm::vec3 lookDirection;
        glm::vec3 lookDirection0;
        glm::vec3 currentTop;
        glm::vec3 Right;
        glm::vec3 WorldTop;
        // Euler Angles
        float Yaw;
        float Pitch;
        // Camera options
        float MovementSpeed;
        float MouseSensitivity;
        float Zoom;

        
        // Constructor with vectors
        camera(glm::vec3 position = glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f), float yaw = YAW, float pitch = PITCH);

        // Constructor with scalar values
        camera(float posX, float posY, float posZ, float upX, float upY, float upZ, float yaw, float pitch);

        // Returns the view matrix calculated using Euler Angles and the LookAt Matrix
        glm::mat4 GetViewMatrix();

        // Processes input received from any keyboard-like input system. Accepts input parameter in the form of camera defined ENUM (to abstract it from windowing systems)
        void ProcessKeyboard(Camera_Movement direction, float deltaTime);

        // Processes input received from a mouse input system. Expects the offset value in both the x and y direction.
        void ProcessMouseMovement(float xoffset, float yoffset, GLboolean constrainPitch/* = true*/);

        // Processes input received from a mouse input system. Expects the offset value in both the x and y direction.
        void ProcessMouseMovementLeftBtn(float xoffset, float yoffset);

        void ProcessMouseMovementRightBtn(float xoffset, float yoffset, float deltaTime);

        // Processes input received from a mouse scroll-wheel event. Only requires input on the vertical wheel-axis
        void ProcessMouseScroll(float yoffset);

    private:
        // Calculates the front vector from the Camera's (updated) Euler Angles
        void updateCameraVectors();
    };

}

#endif