#ifndef LINE_HPP
#define LINE_HPP

#include "gprimitive.hpp"

namespace oglu
{
    class line: public GraphicsPrimitive
    {
        public:

            void use(unsigned int *_vbo, unsigned int *_vao, unsigned int objectId, bool isTexture);

            void render(
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
                        ) override;
            
        private:
            void apply() override;
            void use_shader() override;
    };
}

#endif