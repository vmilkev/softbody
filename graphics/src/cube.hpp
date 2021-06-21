#ifndef CUBE_HPP
#define CUBE_HPP

#include "gprimitive.hpp"

namespace oglu
{
    class cube: public GraphicsPrimitive
    {
        public:

            void use(unsigned int *_vbo, unsigned int *_vao, unsigned int objectId, bool isTexture);            

    };
}

#endif