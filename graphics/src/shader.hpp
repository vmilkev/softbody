#ifndef SHADER_HPP
#define SHADER_HPP

#include "glad/glad.h"
#include "glmath.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

namespace oglu
{
    class shader
    {
        public:
            shader( std::vector<const char *> &file_path );
            shader(const char vertex[], const char fragment[]);
            void apply();
            void setUniformBool(const std::string &name, bool value) const;  
            void setUniformInt(const std::string &name, int value) const;   
            void setUniformFloat(const std::string &name, float value) const;
            void setUniformMatrFloat(const std::string &name, glm::mat4 value) const;
            void setUniformVecFloat(const std::string &name, glm::vec3 value) const;

        private:
            unsigned int programID;

    };

}

#endif
