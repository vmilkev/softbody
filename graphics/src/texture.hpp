#ifndef TEXTURE_HPP
#define TEXTURE_HPP

#include "glad/glad.h"
#include <iostream>


namespace oglu
{
    class texture
    {
        public:
            texture(const char * path);
            void apply();

        private:
            unsigned int textureID;
            int imWidth;
            int imHeight;
            int imNrChannels;
    };

}

#endif
