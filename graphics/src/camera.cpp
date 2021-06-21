#include "camera.hpp"

oglu::camera::camera(glm::vec3 position, glm::vec3 up, float yaw, float pitch) : lookDirection(glm::vec3(0.0f, 0.0f, -1.0f)), MovementSpeed(SPEED), MouseSensitivity(SENSITIVITY), Zoom(ZOOM)
{
    cameraPosition = position;
    WorldTop = up;
    Yaw = yaw;
    Pitch = pitch;

    /* set the initial condition of the world */
    cameraPosition0 = cameraPosition;
    lookDirection0 = lookDirection;

    updateCameraVectors();

}

//--------------------------------------------------------------------------------------

oglu::camera::camera(float posX, float posY, float posZ, float upX, float upY, float upZ, float yaw, float pitch) : lookDirection(glm::vec3(0.0f, 0.0f, -1.0f)), MovementSpeed(SPEED), MouseSensitivity(SENSITIVITY), Zoom(ZOOM)
{

    cameraPosition = glm::vec3(posX, posY, posZ);
    cameraPosition0 = cameraPosition;
    WorldTop = glm::vec3(upX, upY, upZ);
    Yaw = yaw;
    Pitch = pitch;

    /* set the initial condition of the world */
    cameraPosition0 = cameraPosition;
    lookDirection0 = lookDirection;

    updateCameraVectors();
}

//--------------------------------------------------------------------------------------

// Returns the view matrix calculated using Euler Angles and the LookAt Matrix
glm::mat4 oglu::camera::GetViewMatrix()
{

    return glm::lookAt(cameraPosition, cameraPosition + lookDirection, currentTop);
    //                 cam pos., where to look, up dir.
}

//--------------------------------------------------------------------------------------

void oglu::camera::ProcessKeyboard(Camera_Movement direction, float deltaTime)
{
    float velocity = MovementSpeed * deltaTime;

     if (direction == CENTER){
        cameraPosition = cameraPosition0;
        lookDirection = lookDirection0;
        currentTop = WorldTop;
        Yaw = -90.0f;
        Pitch = 0.0f;
        updateCameraVectors();
    }   
    if (direction == FORWARD)
        cameraPosition += lookDirection * velocity;
    if (direction == BACKWARD)
        cameraPosition -= lookDirection * velocity;
    if (direction == LEFT)
        cameraPosition -= Right * velocity;
    if (direction == RIGHT)
        cameraPosition += Right * velocity;
    if (direction == UP)
        cameraPosition += currentTop * velocity;
    if (direction == DOWN)
        cameraPosition -= currentTop * velocity;
    if (direction == VIEWIN){
        if (Zoom >= 1.0f && Zoom <= 90.0f)
            Zoom -= velocity;
        if (Zoom <= 1.0f)
            Zoom = 1.0f;
        if (Zoom >= 90.0f)
            Zoom = 90.0f;
    }
    if (direction == VIEWOUT){
        if (Zoom >= 1.0f && Zoom <= 90.0f)
            Zoom += velocity;
        if (Zoom <= 1.0f)
            Zoom = 1.0f;
        if (Zoom >= 90.0f)
            Zoom = 90.0f;
    }

}

//--------------------------------------------------------------------------------------

void oglu::camera::ProcessMouseMovement(float xoffset, float yoffset, GLboolean constrainPitch)
{
    xoffset *= MouseSensitivity;
    yoffset *= MouseSensitivity;

    Yaw   += xoffset;
    Pitch += yoffset;

    float pitchLimit = 89.0;
    // Make sure that when pitch is out of bounds, screen doesn't get flipped
    if (constrainPitch)
    {
        if ( Pitch > pitchLimit )
            Pitch = pitchLimit;
        if ( Pitch < -pitchLimit )
            Pitch = -pitchLimit;
    }

    // Update lookDirection, Right and currentTop Vectors using the updated Euler angles
    updateCameraVectors();
}

//--------------------------------------------------------------------------------------

void oglu::camera::ProcessMouseMovementLeftBtn(float xoffset, float yoffset)
{
    GLboolean constrainPitch = true;

    xoffset *= MouseSensitivity;
    yoffset *= MouseSensitivity;

    Yaw   += xoffset;
    Pitch += yoffset;

    float pitchLimit = 89.0;
    // Make sure that when pitch is out of bounds, screen doesn't get flipped
    if (constrainPitch)
    {
        if ( Pitch > pitchLimit )
            Pitch = pitchLimit;
        if ( Pitch < -pitchLimit )
            Pitch = -pitchLimit;
    }

    // Update lookDirection, Right and currentTop Vectors using the updated Euler angles
    updateCameraVectors();
}

//--------------------------------------------------------------------------------------

void oglu::camera::ProcessMouseMovementRightBtn(float xoffset, float yoffset, float deltaTime)
{
    float velocity = /*MovementSpeed * */deltaTime * MouseSensitivity * 4.0f;

    cameraPosition += Right * xoffset* velocity;
    cameraPosition += currentTop * yoffset* velocity;
}

//--------------------------------------------------------------------------------------

void oglu::camera::ProcessMouseScroll(float yoffset)
{
    float velocity = MouseSensitivity/*MovementSpeed*/ * yoffset;
    cameraPosition += lookDirection * velocity;
}

//--------------------------------------------------------------------------------------

void oglu::camera::updateCameraVectors()
{
    glm::vec3 front;
    front.x = cos(glm::radians(Yaw)) * cos(glm::radians(Pitch));
    front.y = sin(glm::radians(Pitch));
    front.z = sin(glm::radians(Yaw)) * cos(glm::radians(Pitch));
    lookDirection = glm::normalize(front);
    Right = glm::normalize(glm::cross(lookDirection, WorldTop));
    currentTop = glm::normalize(glm::cross(Right, lookDirection));
}

//--------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------
