#include "glwindow.hpp"

//--------------------------------------------------------------------------------------

float oglu::glwindow::lastX;
float oglu::glwindow::lastY;
bool oglu::glwindow::firstMouse;
oglu::camera *oglu::glwindow::main_camera;
float oglu::glwindow::lastFrame;

//--------------------------------------------------------------------------------------

oglu::glwindow::glwindow(std::string _title)
{
    title = _title;
    lastFrame = 0.0f;
    this->construct();
}

//--------------------------------------------------------------------------------------

oglu::glwindow::~glwindow()
{
}

//--------------------------------------------------------------------------------------

void oglu::glwindow::clear()
{
    delete main_camera;
}

//--------------------------------------------------------------------------------------

void oglu::glwindow::construct()
{
    //main_camera = new camera(glm::vec3(0.05f, 0.05f, 1.0f));
    main_camera = new camera(glm::vec3(0.50f, 1.0f, 5.0f));

    firstMouse = true;
}


//--------------------------------------------------------------------------------------

GLFWwindow *oglu::glwindow::load()
{

    if ( !glfwInit() ){
        std::cout << stderr << "Failed to initialise GLFW" << std::endl;
        exit(1);
    }

    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWmonitor* monitor;
    monitor = glfwGetPrimaryMonitor();
    if (monitor == NULL){
        std::cout<<"Primary monitor is not detected !"<<std::endl;
        exit(1);
    }
    int scr_width, scr_height;
    const GLFWvidmode * mode = glfwGetVideoMode(monitor);
    scr_width = mode->width;
    scr_height = mode ->height;

    lastX = scr_width / 2.0f;
    lastY = scr_height / 2.0f;

    GLFWwindow *window;
    window = glfwCreateWindow( (int)scr_width/2.0, (int)scr_height/2.0,title.c_str(), NULL, NULL );
    if ( window == NULL ){
        std::cout << stderr << "Failed to open OpenGL window" <<std::endl;
        glfwTerminate();
        exit(1);
    }

    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);


    glfwSetKeyCallback(window, key_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    //glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);
    //glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

    /* glad: load all OpenGL function pointers */
    if ( !gladLoadGLLoader((GLADloadproc)glfwGetProcAddress )){
        std::cout << "Failed to initialize GLAD" << std::endl;
        exit(1);
    }

    return window;
}

//--------------------------------------------------------------------------------------

void oglu::glwindow::framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    /* glfw: whenever the window size changed (by OS or user resize) this callback function executes */
    glViewport(0, 0, width, height);

}  

//--------------------------------------------------------------------------------------

void oglu::glwindow::handle_input(GLFWwindow *window/*, float secondsElapsed*/)
{
    
    float currentFrame = glfwGetTime();
    float deltaTime = currentFrame - lastFrame;
    lastFrame = currentFrame;
 
    /* process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly */

    if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
    
    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);


    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        main_camera->ProcessKeyboard(FORWARD, /*secondsElapsed*/deltaTime);

    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        main_camera->ProcessKeyboard(BACKWARD, /*secondsElapsed*/deltaTime);

    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        main_camera->ProcessKeyboard(LEFT, /*secondsElapsed*/deltaTime);

    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        main_camera->ProcessKeyboard(RIGHT, /*secondsElapsed*/deltaTime);

}

//--------------------------------------------------------------------------------------

void oglu::glwindow::key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{

    float currentFrame = glfwGetTime();
    float deltaTime = currentFrame - lastFrame;
    lastFrame = currentFrame;

    if (action == GLFW_RELEASE)
        return;
    else
	{   
        if (key == GLFW_KEY_ESCAPE)
            glfwSetWindowShouldClose(window, true);
        if (key == GLFW_KEY_1)
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        if (key == GLFW_KEY_2)
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		if (key == GLFW_KEY_W)
			main_camera->ProcessKeyboard(FORWARD, deltaTime);
		if (key == GLFW_KEY_LEFT)
			main_camera->ProcessKeyboard(LEFT, deltaTime);
		if (key == GLFW_KEY_S)
			main_camera->ProcessKeyboard(BACKWARD, deltaTime);
		if (key == GLFW_KEY_RIGHT)
			main_camera->ProcessKeyboard(RIGHT, deltaTime);
		if (key == GLFW_KEY_UP)
			main_camera->ProcessKeyboard(UP, deltaTime);
		if (key == GLFW_KEY_DOWN)
			main_camera->ProcessKeyboard(DOWN, deltaTime);
		if (key == GLFW_KEY_C)
            main_camera->ProcessKeyboard(CENTER, deltaTime);
		if (key == GLFW_KEY_KP_ADD)
            main_camera->ProcessKeyboard(VIEWIN, deltaTime);
		if (key == GLFW_KEY_KP_SUBTRACT)
            main_camera->ProcessKeyboard(VIEWOUT, deltaTime);
		if (key == GLFW_KEY_F){
            glfwMaximizeWindow(window);
        }
        if (key == GLFW_KEY_N){
            glfwRestoreWindow(window);
        } 
    }
       
}

//--------------------------------------------------------------------------------------

void oglu::glwindow::mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
     if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE) 
    {
        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_RELEASE)
            return;
    }

    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos;

    lastX = xpos;
    lastY = ypos;

    if ( glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS && glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_RELEASE ){
        
        main_camera->ProcessMouseMovementLeftBtn(xoffset, yoffset);
    }
    else if ( glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE && glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS ){

        float currentFrame = glfwGetTime();
        float deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        main_camera->ProcessMouseMovementRightBtn(xoffset, yoffset, deltaTime);
    }
    
}

//--------------------------------------------------------------------------------------

void oglu::glwindow::scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    /* glfw: whenever the mouse scroll wheel scrolls, this callback is called */
    main_camera->ProcessMouseScroll(yoffset);
}

//--------------------------------------------------------------------------------------
