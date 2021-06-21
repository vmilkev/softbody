#include "plane.hpp"
#include "glmath.hpp"

namespace oglu
{

    //--------------------------------------------------------------

    void plane::use(unsigned int *_vbo, unsigned int *_vao, unsigned int objectId, bool isTexture)
    {
        VBO = _vbo;
        VAO = _vao;
        id = objectId;
        istextured = isTexture;

        drawVert = 6;

        if (istextured) {
            vertexNum = 48;
            vertices = new float [48] {
                // top wall
                -0.5f,  0.5f, -0.5f,  0.0f, 1.0f,  0.0f,  1.0f,  0.0f,
                0.5f,  0.5f, -0.5f,  1.0f, 1.0f,  0.0f,  1.0f,  0.0f,
                0.5f,  0.5f,  0.5f,  1.0f, 0.0f,  0.0f,  1.0f,  0.0f,
                0.5f,  0.5f,  0.5f,  1.0f, 0.0f,  0.0f,  1.0f,  0.0f,
                -0.5f,  0.5f,  0.5f,  0.0f, 0.0f,  0.0f,  1.0f,  0.0f,
                -0.5f,  0.5f, -0.5f,  0.0f, 1.0f,  0.0f,  1.0f,  0.0f
            };

            _shader = new shader(VertexShaderSource, FragmentShaderSource);
            _texture = new texture("textures/wood.jpg");

            apply();
        }
        else{
            vertexNum = 36;
            vertices = new float [36] {
                // top wall
                -0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,
                0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,
                0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,
                0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,
                -0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,
                -0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f
            };

            _shader = new shader(colorVertexShaderSource, colourFragmentShaderSource);

            apply();
        }
        
    }
    
    //--------------------------------------------------------------
}