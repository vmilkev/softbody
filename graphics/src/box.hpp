#ifndef BOX_HPP
#define BOX_HPP

#include "gprimitive.hpp"

namespace oglu
{
    class box: public GraphicsPrimitive
    {
        public:

            void use(unsigned int *_vbo, unsigned int *_vao, unsigned int objectId, bool isTexture);            

    };
}

#endif