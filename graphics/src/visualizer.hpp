#ifndef VISUALIZER_HPP
#define VISUALIZER_HPP

#include "glad/glad.h"
#include <GLFW/glfw3.h>

#include <vector>
#include <bits/stdc++.h> 
#include "glmath.hpp"

#include "glwindow.hpp"
#include "cube.hpp"
#include "box.hpp"
#include "plane.hpp"
#include "line.hpp"
#include "sphere.hpp"

namespace oglu
{
    class visualizer
    {
        public:
            visualizer(/*std::vector < std::vector <float> > &gdata*/)
            {
                isCube = false;
                isBox = false;
                isPlane = false;
                isLine = false;
                isSphere = false;

                //grain_data = gdata;
            }
           
            void show(bool &signal, bool &access, bool &render, std::vector < std::vector <float> > &gdata);

        private:
            struct VisData
            {
                std::vector <glm::vec3> positions;
                int objects;
                std::vector <glm::vec3> scaling;
                std::vector <float> angle;
                std::vector <glm::vec3> rotation;
                std::vector <glm::vec3> color;
                uint identifier;
                                // identifier == 0 for cube
                                // identifier == 1 for box
                                // identifier == 2 for plane
                                // identifier == 3 for line
                                // identifier == 4 for sphere
            };

            void update_data(size_t curr_step, size_t upd_step, std::vector < std::vector <float> > &gdata);

            void getCurrUpdate(std::vector < std::vector <float> > &gdata);

            void pushIn(
                        std::vector <VisData> &into,
                        std::vector <glm::vec3> &pos,
                        std::vector <glm::vec3> &scal,
                        std::vector <glm::vec3> &rot,
                        std::vector <glm::vec3> &col,
                        std::vector <float> &ang,
                        int num,
                        uint id
                        );

            void pushInGrain(std::vector <VisData> &into, std::vector < std::vector <float> > &gdata);

            std::vector <VisData> data;

            bool isCube, isBox, isPlane, isLine, isSphere;
            size_t indCube, indBox, indPlane, indLine, indSphere;

            /* container for grain data: [x y z R P] */
            //std::vector < std::vector <float> > grain_data;

    };
}

#endif