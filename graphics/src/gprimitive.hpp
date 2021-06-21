#ifndef GPRIMITIVE_HPP
#define GPRIMITIVE_HPP

#include "glad/glad.h"
#include <GLFW/glfw3.h>
#include "glwindow.hpp"
#include "shader.hpp"
#include "texture.hpp"
#include "shadervars.hpp"

namespace oglu
{
    class GraphicsPrimitive
    {
        public:

            GraphicsPrimitive();

            virtual ~GraphicsPrimitive();

            virtual void use(
                            unsigned int *_vbo,
                            unsigned int *_vao,
                            unsigned int objectId,
                            bool isTexture
                            ) = 0;

            virtual void render(
                                glwindow window,
                                int width,
                                int height,
                                std::vector <glm::vec3> &positions,
                                int numOfObj,
                                std::vector <glm::vec3> &scaling,
                                std::vector <float> &angle,
                                std::vector <glm::vec3> &rotation,
                                std::vector <glm::vec3> &colour,
                                glm::vec3 light = glm::vec3(1.0f)
                                );

        protected:

            /* the unique ID of object primitive in VAO & VBO */
            unsigned int id;
            /* vertex buffer */
            unsigned int *VBO;
            /* index buffer */
            unsigned int *VAO;
            /* use texture */
            bool istextured;
            /* array of vertices for object primitive */
            float *vertices = NULL;
            /* shader object */
            shader *_shader = NULL;
            /* texture object */
            texture *_texture = NULL;
            /* number of vertices in object primitive */
            uint vertexNum;
            /* number of drawing vertices */
            uint drawVert;
            /* some helper method to call shader and texture objects */
            virtual void use_shader();
            virtual void apply();

    };
}
#endif