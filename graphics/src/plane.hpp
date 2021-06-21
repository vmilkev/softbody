#ifndef PLANE_HPP
#define PLANE_HPP

#include "gprimitive.hpp"

namespace oglu
{
    class plane: public GraphicsPrimitive
    {
        public:

            void use(unsigned int *_vbo, unsigned int *_vao, unsigned int objectId, bool isTexture);            

    };
}

#endif