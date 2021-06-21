
#include <iostream>
#include <cmath>

#include "glad/glad.h"
#include <GLFW/glfw3.h> // let GLFW handle the window and the keyboard

#include "glwindow.hpp"
#include "shader.hpp"
#include "texture.hpp"
#include "glmath.hpp"

// settings
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;
std::string wnd_title("my graphics");


int main (){

    //------------------- Create window ----------------------
    
    GLFWwindow *window;
    window = oglu::LoadWindow(SCR_WIDTH, SCR_HEIGHT, wnd_title);

    // glad: load all OpenGL function pointers
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)){
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    //------- Create a shader object and a shader program ------

    oglu::shader firstShader("shaders/CubeVertShader.vshdr", "shaders/CubeFragmShader.fshdr");
    oglu::shader secondShader("shaders/CubeVertShader.vshdr", "shaders/CubeFragmShader2.fshdr");

    // set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------
    // add a new set of vertices to form a second triangle (a total of 6 vertices); the vertex attribute configuration remains the same (still one 3-float position vector per vertex)
    
    float vertices1[] = {
        // triangle_1 coord. // texture coord. 
        0.5f,  0.5f, 0.0f, 0.0f, 1.0f,
        0.5f, -0.4f, 0.0f, 1.0f, 1.0f,
        -0.4f,  0.5f, 0.0f, 1.0f, 0.0f,
    };
    float vertices2[] = {
        // first triangle
        -0.5f,  -0.5f, 0.0f, 0.0f, 1.0f,
        -0.5f, 0.4f, 0.0f, 0.0f, 0.0f,
        0.4f,  -0.5f, 0.0f, 1.0f, 0.0f,
    };

    /*
    float vertices[] = {
        0.5f,  0.5f, 0.0f,  // top right
        0.5f, -0.5f, 0.0f,  // bottom right
        -0.5f, -0.5f, 0.0f,  // bottom left
        -0.5f,  0.5f, 0.0f   // top left 
    };
    unsigned int indices[] = {  // note that we start from 0!
        0, 1, 3,   // first triangle
        1, 2, 3    // second triangle
    }; 
    */

    // vertex and indeces buffer objects
    unsigned int VBO[2], VAO[2];

    glGenVertexArrays(2, VAO);
    glGenBuffers(2, VBO);
    //glGenBuffers(1, &EBO);

    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
    glBindVertexArray(VAO[0]);
    // 0. copy our vertices array in a buffer for OpenGL to use
    glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices1), vertices1, GL_STATIC_DRAW);
    // 1. then set the vertex attributes pointers
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // texture coord attribute
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindVertexArray(VAO[1]);
    glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices2), vertices2, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);  
    // texture coord attribute
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    //glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);


    // note that this is allowed, the call to glVertexAttribPointer registered VBO as the vertex attribute's bound vertex buffer object so afterwards we can safely unbind
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // You can unbind the VAO afterwards so other VAO calls won't accidentally modify this VAO, but this rarely happens. Modifying other
    // VAOs requires a call to glBindVertexArray anyways so we generally don't unbind VAOs (nor VBOs) when it's not directly necessary.
    glBindVertexArray(0); 

    // textures
    oglu::texture firstTexture("textures/wood.jpg");
    oglu::texture secondTexture("textures/sample.jpg");

    // tell opengl for each sampler to which texture unit it belongs to (only has to be done once)

    firstShader.apply(); // don't forget to activate/use the shader before setting uniforms!
    firstShader.setUniformInt("texture_cube1", 0);

    secondShader.apply(); // don't forget to activate/use the shader before setting uniforms!
    secondShader.setUniformInt("texture2", 0);



    // render loop
    while(!glfwWindowShouldClose(window)){

        // change wireframe polygons modes
        if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        // input
        oglu::processInput(window);

        // rendering

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // apply some simple transformations
        glm::mat4 trans1 = glm::mat4(1.0f);
        glm::mat4 trans2 = glm::mat4(1.0f);
        glm::mat4 model1 = glm::mat4(1.0f);
        glm::mat4 view1 = glm::mat4(1.0f);
        glm::mat4 projection1;

        // first triangle
        
        firstTexture.apply();

        // 2. use our shader program when we want to render an object
        firstShader.apply();

        // update the uniform color
        float timeValue = glfwGetTime();
        float greenValue = sin(timeValue) / 2.0f + 0.5f; // we are going to chang the "green" part of vec4 colour vector
        firstShader.setUniformFloat("modifiedGreenColour", greenValue);
        firstShader.setUniformFloat("offsetX", greenValue);
        firstShader.setUniformFloat("offsetY", -greenValue);

        trans1 = glm::rotate(trans1, glm::radians(90.0f*greenValue), glm::vec3(0.5, 0.5, 1.0));
        trans1 = glm::scale(trans1, glm::vec3(0.8+greenValue, 0.8+greenValue, 0.8+greenValue));
        
        model1 = glm::rotate(model1, glm::radians(-55.0f), glm::vec3(1.0f, 0.0f, 0.0f));         
        // note that we're translating the scene in the reverse direction of where we want to move
        view1 = glm::translate(view1, glm::vec3(0.0f, 0.0f, -3.0f));        
        projection1 = glm::perspective(glm::radians(45.0f), (float) (SCR_WIDTH / SCR_HEIGHT), 0.1f, 100.0f);

        firstShader.setUniformMatrFloat("transform", trans1);
        firstShader.setUniformMatrFloat("model", model1);
        firstShader.setUniformMatrFloat("view", view1);
        firstShader.setUniformMatrFloat("projection", projection1);

        // 3. now draw the object 
        glBindVertexArray(VAO[0]);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        

        // second triangle

        secondTexture.apply();
        secondShader.apply();

        secondShader.setUniformMatrFloat("transform", trans2);
        secondShader.setUniformMatrFloat("transform", trans1);
        secondShader.setUniformMatrFloat("model", model1);
        secondShader.setUniformMatrFloat("view", view1);
        secondShader.setUniformMatrFloat("projection", projection1);

        glBindVertexArray(VAO[1]);
        glDrawArrays(GL_TRIANGLES, 0, 3);

        //glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();    
    }

    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------
    glDeleteVertexArrays(1, VAO);
    glDeleteBuffers(1, VBO);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

    return 0;
}
