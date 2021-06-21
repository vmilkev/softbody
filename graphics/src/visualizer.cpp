#include "visualizer.hpp"

//tmp includes
#include "cubeposition.hpp"
#include "boxposition.hpp"
#include "planeposition.hpp"
#include "lineposition.hpp"
#include "sphereposition.hpp"
//end of tmp includes

namespace oglu
{
    //--------------------------------------------------------------

    void visualizer::pushIn(
                            std::vector <VisData> &into,
                            std::vector <glm::vec3> &pos,
                            std::vector <glm::vec3> &scal,
                            std::vector <glm::vec3> &rot,
                            std::vector <glm::vec3> &col,
                            std::vector <float> &ang,
                            int num,
                            uint id
                            )
    {
        VisData tstr;
        tstr.positions = pos;
        tstr.objects = num;
        tstr.scaling = scal;
        tstr.angle = ang;
        tstr.rotation = rot;
        tstr.color = col;
        tstr.identifier = id;

        into.push_back(tstr);

    }

    //--------------------------------------------------------------

    void visualizer::pushInGrain(std::vector <VisData> &into, std::vector < std::vector <float> > &gdata)
    {
        //std::cout<<"HERE 1"<<std::endl;

        // debugging param
        float scaleNorm = 1.0;

        std::vector <glm::vec3> pos;
        std::vector <glm::vec3> scal;
        std::vector <glm::vec3> rot;
        std::vector <glm::vec3> col;
        std::vector <float> ang;
        uint id = 4;

        size_t grainNum = gdata.size();

        //std::cout<<"grainNum = "<<grainNum<<std::endl;

        for (auto i = 0; i < grainNum; i++){
            glm::vec3 position(gdata[i][0], gdata[i][1], gdata[i][2]);

            //std::cout<<"position = "<<position.x<<", "<<position.y<<", "<<position.z<<std::endl;
            //std::cout<<"position (init) size = "<<gdata[i].size()<<std::endl;
            //std::cout<<"position (init) = "<<gdata[i][0]<<", "<<gdata[i][1]<<", "<<gdata[i][2]<<std::endl;

            glm::vec3 scaling(gdata[i][3], gdata[i][3], gdata[i][3]);
            glm::vec3 rotation(1.0f, 1.0f, 1.0f);
            glm::vec3 color(1.0f, 1.0f-gdata[i][4], 1.0f-gdata[i][4]);
            float angle = gdata[i][5];

            pos.push_back(position);
            scal.push_back(scaling);
            rot.push_back(rotation);
            col.push_back(color);
            ang.push_back(angle);
        }

        VisData tstr;

        tstr.positions = pos;
        tstr.objects = grainNum;
        tstr.scaling = scal;
        tstr.angle = ang;
        tstr.rotation = rot;
        tstr.color = col;
        tstr.identifier = id;

        into.push_back(tstr);

/*                 std::cout<<"into[indSphere].positions.size() = "<<into[0].positions.size()<<std::endl;
                for (auto i = 0; i < into[0].positions.size(); i++)
                    std::cout<<"positions = "<<into[0].positions[i].x<<", "<<into[0].positions[i].y<<", "<<into[0].positions[i].z<<std::endl;         
 */
    }

    //--------------------------------------------------------------

    void visualizer::getCurrUpdate(std::vector < std::vector <float> > &gdata)
    {
        //---------------------------------------------
        std::vector <VisData> tdata;

        pushInGrain(tdata, gdata);

        std::swap(tdata, data);

        //---------------------------------------------

    }
    //--------------------------------------------------------------

    void visualizer::update_data(size_t curr_step, size_t upd_step, std::vector < std::vector <float> > &gdata)
    {
        if(curr_step % upd_step == 0){

            isCube = false;
            isBox = false;
            isPlane = false;
            isLine = false;
            isSphere = false;

            getCurrUpdate(gdata);

            size_t currPrimitiveNum = data.size();

            if (currPrimitiveNum == 0){
                std::cout<<"There are no data have been given to a visualizer! Exit."<<std::endl;
                exit(1);
            }

            for (auto i = 0; i < currPrimitiveNum; i++){
                uint which_obj  = data[i].identifier;
                switch (which_obj)
                {
                case 0:
                    isCube = true;
                    indCube = i;
                    break;
                case 1:
                    isBox = true;
                    indBox = i;
                    break;
                case 2:
                    isPlane = true;
                    indPlane = i;
                    break;
                case 3:
                    isLine = true;
                    indLine = i;
                    break;
                case 4:
                    isSphere = true;
                    indSphere = i;
                    break;
                
                default:
                    break;
                }
            }

/*                 std::cout<<"data[indSphere].positions.size() = "<<data[indSphere].positions.size()<<"; indSphere = "<<indSphere<<std::endl;
                for (auto i = 0; i < data[indSphere].positions.size(); i++)
                    std::cout<<"positions = "<<data[indSphere].positions[i].x<<", "<<data[indSphere].positions[i].y<<", "<<data[indSphere].positions[i].z<<std::endl;         
 */

        }

    }

    //--------------------------------------------------------------

    void visualizer::show(bool &signal, bool &access, bool &render, std::vector < std::vector <float> > &gdata)
    {

        size_t primitivesNum = 5;

        oglu::glwindow main_window("DEM simulation");

        GLFWwindow *whandle = main_window.load();    

        /* vertex and indeces buffer objects */
        unsigned int VBO[primitivesNum], VAO[primitivesNum];

        glGenVertexArrays(primitivesNum, VAO);
        glGenBuffers(primitivesNum, VBO);
            
        //--------------------------------------------------------------
        oglu::GraphicsPrimitive *texturedCube = new oglu::cube();
        texturedCube->use(VBO, VAO, 0, false);
        //--------------------------------------------------------------

        //--------------------------------------------------------------   
        oglu::GraphicsPrimitive *texturedBox = new oglu::box();
        texturedBox->use(VBO, VAO, 1, false);
        //--------------------------------------------------------------

        //--------------------------------------------------------------
        oglu::GraphicsPrimitive *coloredPlane = new oglu::plane();
        coloredPlane->use(VBO, VAO, 2, true);
        //--------------------------------------------------------------

        //--------------------------------------------------------------
        oglu::GraphicsPrimitive *coloredLine = new oglu::line();
        coloredLine->use(VBO, VAO, 3, false);
        //--------------------------------------------------------------

        //--------------------------------------------------------------
        oglu::GraphicsPrimitive *coloredSphere = new oglu::sphere(16);
        coloredSphere->use(VBO, VAO, 4, false);
        //--------------------------------------------------------------

        glEnable(GL_DEPTH_TEST);

        size_t step = 0;
        size_t upd = 1;

/*         if (access){
            access = false;
            update_data(step, upd);
            access = true;
            render = true;
            step++;
        }
 */
        /* rendering loop */
        while(/*!glfwWindowShouldClose(whandle) ||*/ !signal){

            if(glfwWindowShouldClose(whandle))
                signal = true;

            if (access){
                access = false;
                update_data(step, upd, gdata);
                access = true;
                render = true;
            }

            int _width, _height;
            glfwGetWindowSize(whandle, &_width, &_height);

            main_window.handle_input(whandle);

            glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            /* ------------------- START drawing the object: CUBE ------------------- */
            if(isCube){
                texturedCube->render(
                                    main_window,
                                    _width,
                                    _height,
                                    data[indCube].positions,
                                    data[indCube].objects,
                                    data[indCube].scaling,
                                    data[indCube].angle,
                                    data[indCube].rotation,
                                    data[indCube].color
                                    );
            }
            /* ------------------- END drawing the first object: CUBE --------------- */

            /* ------------------- START drawing the object: BOX -------------------- */
            if(isBox){             
                texturedBox->render(
                                    main_window,
                                    _width,
                                    _height,
                                    data[indBox].positions,
                                    data[indBox].objects,
                                    data[indBox].scaling,
                                    data[indBox].angle,
                                    data[indBox].rotation,
                                    data[indBox].color
                                    );
            }
            /* ------------------- END drawing the first object: BOX ---------------- */

            /* ------------------- START drawing the object: PLANE ------------------ */
            if(isPlane){             
                coloredPlane->render(
                                    main_window,
                                    _width,
                                    _height,
                                    data[indPlane].positions,
                                    data[indPlane].objects,
                                    data[indPlane].scaling,
                                    data[indPlane].angle,
                                    data[indPlane].rotation,
                                    data[indPlane].color
                                    );
            }
            /* ------------------- END drawing the first object: BOX ---------------- */

            /* ------------------- START drawing the object: LINE ------------------- */
            if(isLine){
                coloredLine->render(
                                    main_window,
                                    _width,
                                    _height,
                                    data[indLine].positions,
                                    data[indLine].objects,
                                    data[indLine].scaling,
                                    data[indLine].angle,
                                    data[indLine].rotation,
                                    data[indLine].color
                                    );
            }
            /* ------------------- END drawing the first object: LINE --------------- */

            /* ------------------- START drawing the object: SPHERE ----------------- */
            if(isSphere){
/*                 std::cout<<"data[indSphere].positions.size() = "<<data[indSphere].positions.size()<<"; indSphere = "<<indSphere<<std::endl;
                for (auto i = 0; i < data[indSphere].positions.size(); i++)
                    std::cout<<"positions = "<<data[indSphere].positions[i].x<<", "<<data[indSphere].positions[i].y<<", "<<data[indSphere].positions[i].z<<std::endl;         
 */                
                coloredSphere->render(
                                    main_window,
                                    _width,
                                    _height,
                                    data[indSphere].positions,
                                    data[indSphere].objects,
                                    data[indSphere].scaling,
                                    data[indSphere].angle,
                                    data[indSphere].rotation,
                                    data[indSphere].color
                                    );
            }
            /* ------------------- END drawing the first object: SPHERE ------------- */

            glfwSwapBuffers(whandle);
            glfwPollEvents();

            step++;
        }

        glDeleteVertexArrays(primitivesNum, VAO);
        glDeleteBuffers(primitivesNum, VBO);
        
        /* Close OpenGL window, clear memory and terminate GLFW */

        main_window.clear();

        delete texturedCube;
        delete texturedBox;
        delete coloredLine;
        delete coloredPlane;

        glfwTerminate();

    }

}