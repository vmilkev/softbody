#ifndef SPHERE_HPP
#define SPHERE_HPP

#include "gprimitive.hpp"

namespace oglu
{
    class sphere: public GraphicsPrimitive
    {
        public:
            sphere(uint prec = 5, bool light = false)
            {
                precision = prec;
                isLight = light;
            }

            void use(unsigned int *_vbo, unsigned int *_vao, unsigned int objectId, bool isTexture);            

        private:
            void apply() override;
            void use_shader() override;
            void get_verteces();
            uint precision = 0;
            bool isLight;

    };
}

#endif