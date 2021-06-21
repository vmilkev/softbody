#ifndef GLWINDOW_HPP
#define GLWINDOW_HPP

#include "glad/glad.h"
#include <GLFW/glfw3.h>
#include <iostream>
#include <string>
#include "camera.hpp"

namespace oglu
{
    class glwindow
    {
        public:

            static camera *main_camera;

            /* CONSTRUCTOR */
            glwindow(std::string _title);
            ~glwindow();            
            GLFWwindow *load();            
            void handle_input(GLFWwindow *window);
            void clear();

        private:

            static float lastX;
            static float lastY;
            static bool firstMouse;
            static float lastFrame;
            //unsigned int scr_width;
            //unsigned int scr_height;
            std::string title;

            static void construct();

            /* CALLBACKs */
            static void framebuffer_size_callback(GLFWwindow* window, int width, int height);
            static void mouse_callback(GLFWwindow* window, double xpos, double ypos);            
            static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
            static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);

    };
}

#endif
