#include "sphere.hpp"
#include "glmath.hpp"
#include <vector>

namespace oglu
{
    //--------------------------------------------------------------

    void sphere::use(unsigned int *_vbo, unsigned int *_vao, unsigned int objectId, bool isTexture)
    {
        VBO = _vbo;
        VAO = _vao;
        id = objectId;
        istextured = isTexture;

        /* the size of vertex array calculated as follows:
           size = (2*precision+1)*(precision-1)*18
        */
        switch (precision)
        {
            case 6:
                vertexNum = 1170;
                vertices = new float [1170];
                break;
            case 7:
                vertexNum = 1620;
                vertices = new float [1620];
                break;
            case 8:
                vertexNum = 2142;
                vertices = new float [2142];
                break;
            case 9:
                vertexNum = 2726;
                vertices = new float [2736];
                break;
            case 10:
                vertexNum = 3402;
                vertices = new float [3402];
                break;
            case 11:
                vertexNum = 4140;
                vertices = new float [4140];
                break;
            case 12:
                vertexNum = 4950;
                vertices = new float [4950];
                break;
            case 13:
                vertexNum = 5832;
                vertices = new float [5832];
                break;
            case 14:
                vertexNum = 6786;
                vertices = new float [6786];
                break;
            case 15:
                vertexNum = 7812;
                vertices = new float [7812];
                break;
            case 16:
                vertexNum = 8910;
                vertices = new float [8910];
                break;
            
            default:
                if(precision < 6){
                    precision = 5;
                    vertexNum = 792;
                    vertices = new float [792];
                }
                else if(precision > 16){
                    precision = 16;
                    vertexNum = 8910;
                    vertices = new float [8910];
                }
                break;
        }

        get_verteces();

        if(isLight)
            _shader = new shader(VertexShaderSource, lightFragmentShaderSource);
        else
            _shader = new shader(colorVertexShaderSource, colourFragmentShaderSource);
        
        apply();
        
    }

    //--------------------------------------------------------------

    void sphere::get_verteces()
    {
        int vertThau = precision + 1;
        std::vector <glm::vec3> vertList;
        glm::vec3 tmpVect;

        drawVert = precision*(precision-1)*12*2; // here I multiply by 2 becaus I reuse vertices as a normals

        double phi = 0.0;
        while(phi <= 2*3.141592654){
            double theta = 0.0;
            while(theta <= 3.141592654){
                tmpVect.x = 0.5*sin(theta)*cos(phi);
                tmpVect.y = 0.5*cos(theta);
                tmpVect.z = 0.5*sin(theta)*sin(phi);
                vertList.push_back(tmpVect);
                theta = theta + 3.14159265/precision; 
            }
            phi = phi + 3.14159265/precision;
        }

        uint v = 0;
        uint ii = 0;
        for(auto i = 0; i < vertList.size(); i = i + vertThau){
            uint j = 0;ii++;
            while(j < precision-1){

                vertices[v] = vertList[i+j].x; v++;
                vertices[v] = vertList[i+j].y; v++;
                vertices[v] = vertList[i+j].z; v++;

                vertices[v] = vertList[i+j+1].x; v++;
                vertices[v] = vertList[i+j+1].y; v++;
                vertices[v] = vertList[i+j+1].z; v++;

                vertices[v] = vertList[i+j+2+precision].x; v++;
                vertices[v] = vertList[i+j+2+precision].y; v++;
                vertices[v] = vertList[i+j+2+precision].z; v++;
                
                j++;
            }
            j = 0;
            while(j < precision-1){
                
                vertices[v] = vertList[i+j+1].x; v++;
                vertices[v] = vertList[i+j+1].y; v++;
                vertices[v] = vertList[i+j+1].z; v++;
                
                vertices[v] = vertList[i+j+2+precision].x; v++;
                vertices[v] = vertList[i+j+2+precision].y; v++;
                vertices[v] = vertList[i+j+2+precision].z; v++;

                vertices[v] = vertList[i+j+3+precision].x; v++;
                vertices[v] = vertList[i+j+3+precision].y; v++;
                vertices[v] = vertList[i+j+3+precision].z; v++;
                
                j++;
            }
        }
        vertList.clear();
        vertList.shrink_to_fit();
    }

    //--------------------------------------------------------------

    void sphere::apply()
    {
        // bind the Vertex Array Object
        glBindVertexArray(VAO[id]);
        glBindBuffer(GL_ARRAY_BUFFER, VBO[id]);

        glBufferData(GL_ARRAY_BUFFER, vertexNum * sizeof(float), vertices, GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(1);

        _shader->apply();

    }

    //--------------------------------------------------------------

    void sphere::use_shader()
    {
        _shader->apply();
    }

    //--------------------------------------------------------------
}