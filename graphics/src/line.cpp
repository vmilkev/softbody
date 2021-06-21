#include "line.hpp"
#include "glmath.hpp"

namespace oglu
{
    //--------------------------------------------------------------

    void line::use(unsigned int *_vbo, unsigned int *_vao, unsigned int objectId, bool isTexture)
    {
        VBO = _vbo;
        VAO = _vao;
        id = objectId;
        istextured = isTexture;

        drawVert = 2;

        vertexNum = 6;
        vertices = new float [6];
        _shader = new shader(VertexShaderSource, colourFragmentShaderSource);
        
    }
    
    //--------------------------------------------------------------

    void line::apply()
    {
        // bind the Vertex Array Object
        glBindVertexArray(VAO[id]);
        glBindBuffer(GL_ARRAY_BUFFER, VBO[id]);
        glBufferData(GL_ARRAY_BUFFER, vertexNum * sizeof(float), vertices, GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        _shader->apply();

    }
    //--------------------------------------------------------------

    void line::use_shader()
    {
        _shader->apply();

    }

    //--------------------------------------------------------------

    void line::render(
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
        unsigned int j = 0;

        for(unsigned int i = 0; i < 2*numOfObj; i = i + 2){

            vertices[0] = positions[i].x;
            vertices[1] = positions[i].y;
            vertices[2] = positions[i].z;
            vertices[3] = positions[i+1].x;
            vertices[4] = positions[i+1].y;
            vertices[5] = positions[i+1].z;
            
            apply();
 
            use_shader();

            glBindVertexArray(VAO[id]);

            /* calculate projection matrix (note that in this case it could change every frame) */
            glm::mat4 projection = glm::perspective(glm::radians(window.main_camera->Zoom), (float)width / (float)height, 0.1f, 100.0f);
            
            /* calculate camera/view transformation matrix */
            glm::mat4 view = window.main_camera->GetViewMatrix();
                        
            /* calculate model matrix (for each object) */
            glm::mat4 model = glm::mat4(1.0f);

            /* calculate transformation matrix (for each object) */
            glm::mat4 transform = glm::mat4(1.0f);

            /* pass matrix variables (projection, view, model and transformation) to the shader before drawing */
            _shader->setUniformMatrFloat("Projection", projection);
            _shader->setUniformMatrFloat("View", view);
            _shader->setUniformMatrFloat("Model", model);
            _shader->setUniformMatrFloat("Transform", transform);

            glm::vec3 currColor = glm::vec3(colour[j].x, colour[j].y, colour[j].z);
                        
            glm::vec3 currLight = light;
            glm::vec3 currLPos = window.main_camera->cameraPosition;
            glm::vec3 currCamPos = window.main_camera->cameraPosition;
            
            _shader->setUniformVecFloat("viewPos", currCamPos);
            _shader->setUniformVecFloat("objectColor", currColor);
            _shader->setUniformVecFloat("lightColor", currLight);            
            _shader->setUniformVecFloat("lightPos", currLPos);

            glLineWidth(scaling[j].x);

            /* drawing */
            glDrawArrays(GL_LINES, 0, drawVert);

            j++;

        }

    }

    //--------------------------------------------------------------
}