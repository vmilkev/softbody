#include "gprimitive.hpp"
#include "glmath.hpp"

namespace oglu
{
    //--------------------------------------------------------------

    GraphicsPrimitive::GraphicsPrimitive()
    {
    }

    //--------------------------------------------------------------

    GraphicsPrimitive::~GraphicsPrimitive()
    {
        if(vertices != NULL)
            delete vertices;
        if(_shader != NULL)
            delete _shader;
        if(_texture != NULL)
            delete _texture;
    }

    //--------------------------------------------------------------

    void GraphicsPrimitive::render(
                                    glwindow window,
                                    int width,
                                    int height,
                                    std::vector <glm::vec3> &positions,
                                    int numOfObj,
                                    std::vector <glm::vec3> &scaling,
                                    std::vector <float> &angle,
                                    std::vector <glm::vec3> &rotation,
                                    std::vector <glm::vec3> &colour,
                                    glm::vec3 light
                                    )
    {
        use_shader();

        glBindVertexArray(VAO[id]);

        /* loop to draw multiple boxs in different positions */
        for(unsigned int i = 0; i < numOfObj; i++){

            /* calculate projection matrix (note that in this case it could change every frame) */
            glm::mat4 projection = glm::perspective(glm::radians(window.main_camera->Zoom), (float)width / (float)height, 0.01f, 100.0f);
            
            /* calculate camera/view transformation matrix */
            glm::mat4 view = window.main_camera->GetViewMatrix();
                        
            /* calculate model matrix (for each object) */
            glm::mat4 model = glm::mat4(1.0f);
            model = glm::translate(model, positions[i]);

            //std::cout<<"colour = "<<colour[i].x<<", "<<colour[i].y<<", "<<colour[i].z<<std::endl;
            
            if (angle[2*i+1])
                model = glm::rotate(model, (float)glfwGetTime() * glm::radians(angle[2*i]), rotation[i]);
            else
                model = glm::rotate(model, glm::radians(angle[2*i]), rotation[i]);
            
            /* calculate transformation matrix (for each object) */
            glm::mat4 transform = glm::mat4(1.0f);

            transform = glm::scale(transform, scaling[i]);

            /* pass matrix variables (projection, view, model and transformation) to the shader before drawing */
            _shader->setUniformMatrFloat("Projection", projection);
            _shader->setUniformMatrFloat("View", view);
            _shader->setUniformMatrFloat("Model", model);
            _shader->setUniformMatrFloat("Transform", transform);

            glm::vec3 currColor = glm::vec3(colour[i].x, colour[i].y, colour[i].z);
            //glm::vec3 currLight = light;
            glm::vec3 currLight = glm::vec3(0.8f);

            glm::vec3 currLPos = window.main_camera->cameraPosition;
            glm::vec3 currCamPos = window.main_camera->cameraPosition;
            
            _shader->setUniformVecFloat("viewPos", currCamPos);
            _shader->setUniformVecFloat("objectColor", currColor);
            _shader->setUniformVecFloat("lightColor", currLight);            
            _shader->setUniformVecFloat("lightPos", currLPos);

            /* drawing */
            glDrawArrays(GL_TRIANGLES, 0, drawVert);

        }

    }

    //--------------------------------------------------------------

    void GraphicsPrimitive::apply()
    {
        // bind the Vertex Array Object
        glBindVertexArray(VAO[id]);
        glBindBuffer(GL_ARRAY_BUFFER, VBO[id]);

        if (istextured){
            glBufferData(GL_ARRAY_BUFFER, vertexNum * sizeof(float), vertices, GL_STATIC_DRAW);
            // set the vertex attributes pointers
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
            glEnableVertexAttribArray(0);

            // set the texture attribute
            glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
            glEnableVertexAttribArray(1);

            // set the normals attribute
            glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(5 * sizeof(float)));
            glEnableVertexAttribArray(2);

            _shader->apply();
            _shader->setUniformInt("texture_1", 0);
        }
        else {
            glBufferData(GL_ARRAY_BUFFER, vertexNum * sizeof(float), vertices, GL_STATIC_DRAW);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
            glEnableVertexAttribArray(0);

            // set the normals attribute
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
            glEnableVertexAttribArray(1);

            _shader->apply();
        }
    }

    //--------------------------------------------------------------

    void GraphicsPrimitive::use_shader()
    {
        _shader->apply();

        if (istextured)
            _texture->apply();
    }

    //--------------------------------------------------------------

}
